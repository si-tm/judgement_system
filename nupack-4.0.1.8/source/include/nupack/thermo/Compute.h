#pragma once
#include "Schedule.h"
#include "Rigs.h"
#include "Engine.h"
#include "Constants.h"
#include "Subopt.h"
#include "Sample.h"
#include <nupack/common/Costs.h>

namespace nupack::thermo {

/******************************************************************************************/

template <class Lengths>
std::size_t subblock_memory(Lengths const &l) {
    return 10 * front(l) * back(l) * (len(l) == 1 ? 1 : 2);
}

/******************************************************************************************/

struct JobRequests {
    Sink sink;
    std::size_t n_structures = 0;
    real energy_gap = 0;
    Sparsity sparsity;
    bool pairs = false;

    NUPACK_REFLECT(JobRequests, sink, n_structures, energy_gap, sparsity, pairs);
};

/******************************************************************************************/

struct PairWork {
    Mat<real> matrix; // matrix to write into
    Sink sink; // where to send the results after
    Complex complex;
    Sparsity sparsity;
    std::atomic<std::size_t> n_remaining; // number of remaining subblocks

    PairWork(Sink f, Complex c, Sparsity s)
        : matrix(sum(c, len), sum(c, len)), sink(std::move(f)), complex(std::move(c)),
          sparsity(std::move(s)), n_remaining(sq(len(complex)) + len(complex)) {}
};

/******************************************************************************************/

struct PairUpdate {
    std::shared_ptr<PairWork> work;
    std::size_t row, col; // starting row and column that should be written to
};

/******************************************************************************************/

template <class O, class T, class Converter>
O &transfer_item(Converter const &c, Item &item) {
    if (auto ptr = item.get<O>()) return *ptr;
    auto old = item.template take<T>();
    auto &p = item.template emplace<O>();
    for_each_zip(members_of(p), members_of(old), [&](auto &x, auto &y) {
        if constexpr(std::is_integral_v<std::decay_t<decltype(x)>>) x = y;
        else x.assign_and_clear(std::move(y), c);
    });
    return p;
}

/******************************************************************************************/

// Try to extract an item in the requested format
template <class O, class X>
auto try_extract(X &x) {
    Guard::ReadLock lock(x.guard);
    return x.data.template get<O>();
}

/******************************************************************************************/

// Extract an item in the requested format. Convert it if it is in the wrong format.
template <class O, class T, class Converter, class X>
auto extract_item(Converter const &c, X &x) {
    O *p = try_extract<O>(x);
    if (!p) {
        // if upgrade guard can be acquired, translate item to new strategy
        Guard::UpgradeLock guard(x.guard);
        if (guard) p = &transfer_item<O, T>(c, x.data);
        else {guard.wait(); p = &x.data.template ref<O>();}
    }
    return p;
}

/******************************************************************************************/

template <class Model>
auto reserve_backup_model(ModelData &mod) {
    NUPACK_QUICK_REQUIRE(mod.capacity, >, 0);
    std::shared_lock model_lock(mod.mutex);
    auto &model = mod.backup.ref<Model>();
    if (model.capacity() < mod.capacity) {
        model_lock.unlock();
        {std::unique_lock unique(model_lock); model.reserve(mod.capacity);}
        model_lock.lock();
    }
    return std::make_pair(&model, std::move(model_lock));
}

/******************************************************************************************/

// Try gathering dependencies with the non-backup strategy
template <class S, bool Forward, class D, class SequenceList, class F>
bool try_with_constants(BlockData &block, D const &deps, ModelData &mod, SequenceList const &v, Action const &action, F &&f) {
    auto subblock = try_extract<typename S::block>(block);
    auto strands = vmap(deps.strands, [&](auto const &s) {return std::make_pair(try_extract<typename S::strand>(*s), &s->guard);});
    auto blocks = vmap(deps.blocks, [&](auto const &b) {return std::make_pair(try_extract<typename S::block>(*b), &b->guard);});
    if (!subblock || !all_of(strands, first_of) || !all_of(blocks, first_of)) return false;

    auto const &model = mod.data.ref<typename S::model>();
    NUPACK_QUICK_REQUIRE(model.capacity(), >=, 2 * *max_element(indirect_view(v, len)), "model reservation failed");
    using Algebra = if_t<Forward, typename S::forward, typename S::backward>;
    if (n_strands(v) == 1) {
        using C = SingleConstants<S::complexity::value, typename S::block, typename S::model, Algebra>;
        f(C(*subblock, view(front(v)), model, action, *front(strands).first));
    } else {
        using C = MultiConstants<S::complexity::value, typename S::block, typename S::model, Algebra>;
        f(C(*subblock, v, model, action, std::move(blocks), std::move(strands)));
    }
    return true;
}

template <class S, bool Forward, class D, class SequenceList, class F>
void with_constants(BlockData &block, D const &deps, ModelData &mod, SequenceList const &v, Action const &action, F &&f) {
    bool computed = try_with_constants<S, Forward>(block, deps, mod, v, action, f);
    if (!computed) if constexpr(S::backup::value) {
        using B = typename S::backup;
        auto subblock = extract_item<typename B::block, typename S::block>(typename S::converter(), block);
        auto strands = vmap(deps.strands, [&](auto const &s) {
            return std::make_pair(extract_item<typename B::strand, typename S::strand>(typename S::converter(), *s), &s->guard);
        });
        auto blocks = vmap(deps.blocks, [&](auto const &b) {
            return std::make_pair(extract_item<typename B::block, typename S::block>(typename S::converter(), *b), &b->guard);
        });
        NUPACK_QUICK_ASSERT(subblock && all_of(blocks, first_of) && all_of(strands, first_of), subblock, blocks, strands);
        auto [model, model_lock] = reserve_backup_model<typename B::model>(mod);

        NUPACK_QUICK_REQUIRE(model->capacity(), >=, 2 * *max_element(indirect_view(v, len)), "model reservation failed");

        using Algebra = if_t<Forward, typename B::forward, typename B::backward>;
        if (n_strands(v) == 1) {
            using C = SingleConstants<S::complexity::value, typename B::block, typename B::model, Algebra>;
            f(C(*subblock, view(front(v)), *model, action, *front(strands).first));
        } else {
            using C = MultiConstants<S::complexity::value, typename B::block, typename B::model, Algebra>;
            f(C(*subblock, v, *model, action, std::move(blocks), std::move(strands)));
        }
        computed = true;
    }
    NUPACK_ASSERT(computed, v, deps.strands, deps.blocks);
}

/******************************************************************************************/

template <class S, class Blocks, class Strands, NUPACK_IF(!S::is_partition_function::value)>
Result::Subopt run_backtrack(Job::Subopt const &job, Complex const &k, ModelData &model, Blocks blocks, Strands strands) {
    if (job.gap < 0) return {};
    std::optional<real> mfe;
    // Initialize queues with Q element
    SuboptQueue q(k);
    q.initialize(Segment(0, back(k).size() - 1, Q.priority()), PairList(sum(k, len)));
    // Proceed largest to smallest blocks
    for (auto o : ~range(k.size())) for (auto i : range(k.size() - o)) {
        auto const j = i + o;
        Dependencies deps;
        deps.strands = view(strands, i, j+1);
        for (auto n : range(o)) { // descending in size
            deps.blocks.emplace_back(blocks(i, j-n-1));
            deps.blocks.emplace_back(blocks(i+n+1, j));
        }
        NUPACK_ASSERT(all_of(deps.blocks) && all_of(deps.strands), "missing dependency");
        bool exit_now = false;
        with_constants<S, false>(*blocks(i, j), deps, model, view(k, i, j+1), blocks(i, j)->action, [&](auto const &c) {
            if (!mfe) exit_now = mfe.emplace(c.model.complex_result(c.block.Q.corner(0, 1), k)) == inf<real>();
            if (exit_now) return;
            bool finished = q.consume(job.gap + 1e-4, job.max_number, c, i, j);
            NUPACK_ASSERT(finished, "Too many suboptimal structures found", job, i, j);
        });
        if (exit_now) return {};
    }

    auto const &energy_model = model.data.ref<typename S::model>().energy_model;
    auto strucs = vmap(q.structures, [&](auto &s){
        real stack_energy = s.gap + *mfe;
        real energy = has_subensemble(energy_model.ensemble) ? energy_model.structure_energy(k, s.structure) : stack_energy;
        return StructureEnergy{std::move(s.structure), energy, stack_energy};
    });
    for (auto const &s : strucs) NUPACK_QUICK_ASSERT(std::isfinite(s.stack_energy) && std::isfinite(s.energy), s, k, job, *mfe);

    spreadsort_float_map(strucs, [](auto const &x) {return x.energy;});
    return {std::move(strucs)};
}

/******************************************************************************************/

template <class S, class Blocks, class Strands, NUPACK_IF(S::is_partition_function::value)>
Result::Sample run_backtrack(Job::Sample const &job, Complex const &k, ModelData &model, Blocks blocks, Strands strands) {
    NUPACK_REQUIRE(job.number, >, 0);
    vec<PairList> samples(job.number, PairList(sum(k, len)));
    DefaultRNG rng(job.seed);
    SampleQueue q(k);
    // Initialize queues with Q element
    q.queues(0, k.size()-1).emplace(Segment(0, back(k).size() - 1, Q.priority()), range(job.number));
    // Proceed largest to smallest blocks
    bool finite = true;
    for (auto o : ~range(k.size())) for (auto i : range(k.size() - o)) {
        if (!finite) continue;
        auto const j = i + o;
        Dependencies deps;
        deps.strands = view(strands, i, j+1);
        for (auto n : range(o)) { // descending in size
            deps.blocks.emplace_back(blocks(i, j-n-1));
            deps.blocks.emplace_back(blocks(i+n+1, j));
        }
        NUPACK_ASSERT(all_of(deps.blocks) && all_of(deps.strands), "missing dependency");
        with_constants<S, false>(*blocks(i, j), deps, model, view(k, i, j+1), blocks(i, j)->action, [&](auto const &c) {
            if (i == 0 && j+1 == k.size() && !std::isfinite(S::forward::rig_type::as_logarithm(c.block.Q.corner(0, 1)))) finite = false;
            else q.consume(samples, c, i, j, rng);
        });
    }
    random_shuffle(samples, rng);
    return {std::move(samples)};
}

/******************************************************************************************/

Result::Pairs finalize_pairs(Mat<real> P, Job::Pairs const &p);
Result::Costs finalize_pairs(Mat<real> P, Job::Costs const &p);

/******************************************************************************************/

template <class S, class Model>
Computer make_computer(Model model) {
    auto const alphabet = model.alphabet();
    Computer com;
    com.cache = Cache::create();
    // Prepare models for computation
    com.prepare = [model=std::move(model)](tf::Taskflow &flow, std::shared_ptr<ModelData> mod, SharedError err) {
        return flow.emplace([err=std::move(err), model=model, mod=std::move(mod)]() mutable {
            err.invoke_noexcept([&] {
                using M = typename S::model;

                if (!mod->data) mod->data = Item::make<M>(typename M::model_type(model));
                mod->data.template ref<M>().reserve(mod->capacity); // reserve done ahead of time

                if constexpr(S::backup::value) {
                    using M = typename S::backup::model;
                    if (!mod->backup) mod->backup = Item::make<M>(typename M::model_type(model)); // reserve is done on demand
                }
            });
        });
    };
    // Compute strand data
    com.compute_strand = [](StrandData &strand, auto const &mod) {
        if (strand.data) return;
        auto const &model = mod.data.template ref<typename S::model>();
        auto &data = strand.data.template emplace<typename S::strand>();
        data.calculate(strand.sequence, model);
    };
    // Compute block data
    com.compute_block = [](tf::Subflow &flow, BlockData &block, SharedError &err, ModelData &mod, auto const &deps, bool full_requested) {
        int const end = len(block.complex) == 1 || !full_requested ? len(block.complex[0]) : len(block.complex[0]) + len(block.complex.back()) - 1;

        bool fresh = !block.data;
        if (fresh) block.data.emplace<typename S::block>();
        else if (block.stat.value >= end) return;

        // another layer of subflow needed because forward_pass joins its given subflow
        auto task = flow.emplace([&block, fresh, err=err, end, &deps, &mod](tf::Subflow &flow) mutable {
            if (!err.is_set()) err.invoke_noexcept([&] {
                with_constants<S, true>(block, deps, mod, block.complex, block.action, [&](auto const &c) {
                    if (fresh) allocate_block(c.block, c);
                    // Use this for profiling memory
                    // print(block.complex, c.model.energy_model.ensemble, real(memory::measure(c.block)) / S::complex_memory(block.complex, c.model.energy_model.alphabet()));
                    block.stat = c.forward_pass(flow, err, block.stat.start_diagonal(), end, false);
                    err.rethrow_if_set();
                });
            });
        });

        if constexpr(S::backup::value) task = flow.emplace([&block, &mod, &deps, end, err=err](tf::Subflow &flow) mutable {
            if (!err.is_set()) err.invoke_noexcept([&] {
                if (block.stat.value >= end) return;
                block.guard.set_stage(true);
                transfer_item<typename S::backup::block, typename S::block>(typename S::converter(), block.data); // no lock necessary

                with_constants<S, true>(block, deps, mod, block.complex, block.action, [&](auto const &c) {
                    block.stat = c.forward_pass(flow, err, block.stat.start_diagonal(), end, true);
                    err.rethrow_if_set();
                });
            });
        }).succeed(task);

        flow.join();
        NUPACK_REQUIRE(block.stat.value, ==, end, "Dynamic programming block computation failed", block.complex, TypeName<S>());
    };
    // Compute block memory
    com.block_memory = [alphabet](Complex const &x) {return S::complex_memory(x, alphabet);};
    // Compute strand memory
    com.strand_memory = [alphabet](Sequence const &s) {return S::strand_memory(s, alphabet);};
    // Compute partition function
    com.extensions.emplace_back([](CallbackList &callbacks, Item input, auto function) {
        auto key = input.take<Complex>();
        callbacks.emplace_back(Root{.blocks={{std::move(key), true}}}, [f=std::move(function)](Ignore, auto deps, auto const &model) {
            real raw;
            {
                if (auto const p = try_extract<typename S::block>(*deps.blocks[0]); p) {
                    raw = S::forward::rig_type::as_logarithm(p->Q.corner(0, 1));
                } else {
                    if constexpr(S::backup::value) {
                        using B = typename S::backup;
                        auto const p = try_extract<typename B::block>(*deps.blocks[0]);
                        NUPACK_ASSERT(p);
                        raw = B::forward::rig_type::as_logarithm(p->Q.corner(0, 1));
                    } else NUPACK_ERROR("missing dependency");
                }
            }
            using T = if_t<S::is_partition_function::value, Result::PF, Result::MFE>;
            T res{model.data.template ref<typename S::model>().complex_result(raw, deps.blocks[0]->complex), raw};
            f(Item::from(std::make_pair(deps.blocks[0]->complex, std::move(res))));
        });
    });
    // Compute backtrack
    com.extensions.emplace_back([](CallbackList &callbacks, Item input, auto function) {
        if constexpr(decltype(S::ensemble::can_backtrack())::value) {
            using T = if_t<S::is_partition_function::value, Job::Sample, Job::Subopt>;
            auto [job, key] = input.take<std::pair<T, Complex>>();
            Root root;
            for (auto j : indices(key)) for (auto i : range(j+1))
                root.blocks.emplace_back(view(key, i, j+1), true);
            root.strands = key;
            callbacks.emplace_back(std::move(root), [f=std::move(function), key=key, job=job](Ignore, auto deps, auto &model) {
                Triangle<std::shared_ptr<BlockData>> tri(len(key));
                tri.values = std::move(deps.blocks);
                auto res = run_backtrack<S>(std::move(job), key, model, std::move(tri), std::move(deps.strands));
                f(Item::from(std::make_pair(key, std::move(res))));

            });
        } else NUPACK_ERROR("backtracking not supported right now");
    });
    // Compute pair probabilities
    com.extensions.emplace_back([](CallbackList &callbacks, Item input, auto function) {
        using T = if_t<S::is_partition_function::value, Job::Pairs, Job::Costs>;
        auto [job, key] = input.take<std::pair<T, Complex>>();
        auto const n = sum(key, len);

        auto const p = std::make_shared<std::tuple<Mat<real>, T, decltype(function), Complex, std::atomic<std::size_t>>>(
            Mat<real>(n, n), std::move(job), std::move(function), key, len(key) * (len(key) + 1));

        auto const dup = duplicate(key);
        auto const pre = prefixes(true, indirect_view(key, len));

        for (auto I : indices(key)) for (auto J : range(I, I+len(key)+1)) {
            Root root{.blocks={{view(dup, I, J+1), J != I + len(key)}}};

            callbacks.emplace_back(std::move(root), [p, is=pre[I], js=pre[J % len(key)], I, J](Ignore, auto deps, Ignore) {
                auto &[P, job, f, key, count] = *p;
                bool ok;
                {
                    auto Q = try_extract<typename S::block>(*deps.blocks.at(0));
                    if ((ok = bool(Q))) {
                        for (auto i : range(Q->B.shape()[0]))
                            for (auto j : range(I == J ? i+1 : 0, I+len(key) == J ? i : Q->B.shape()[1]))
                                P(is + i, js + j) = S::forward::rig_type::as_logarithm(Q->B(i, j));

                        if (I == 0 && J == len(key)-1) P(0, 0) = S::forward::rig_type::as_logarithm(Q->Q.corner(0, 1));
                    }
                }
                if (!ok) if constexpr(S::backup::value) { // sorry for the copy paste ...
                    auto Q = try_extract<typename S::backup::block>(*deps.blocks.at(0));
                    if ((ok = bool(Q))) {
                        for (auto i : range(Q->B.shape()[0]))
                            for (auto j : range(I == J ? i+1 : 0, I+len(key) == J ? i : Q->B.shape()[1]))
                                P(is + i, js + j) = S::backup::forward::rig_type::as_logarithm(Q->B(i, j));

                        if (I == 0 && J == len(key)-1) P(0, 0) = S::backup::forward::rig_type::as_logarithm(Q->Q.corner(0, 1));
                    }
                }
                NUPACK_ASSERT(ok, "missing dependency");

                if (count.fetch_sub(1, std::memory_order_relaxed) == 1) {
                    f(Item::from(std::make_pair(std::move(key), finalize_pairs(std::move(P), std::move(job)))));
                }
            });
        }
    });

    return com;
}

/******************************************************************************************/

}
