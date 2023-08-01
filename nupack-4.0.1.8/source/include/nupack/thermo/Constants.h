#pragma once
#include "Action.h"
#include "../types/Sequence.h"
#include <nupack/model/Model.h>

namespace nupack::thermo {

/******************************************************************************************/

template <class F>
Stat iterate_diagonal(tf::Subflow &flow, SharedError &err, uint start, uint end, F const &f) noexcept {
    // Proceed through tasks as long as ok is true.
    AtomicStat stat;
    tf::Task last;
    bool stop = false;

    for (auto d : range(start, end)) {
        tf::Task next = flow.emplace([&f, &err, &stat, &stop, d] (tf::Subflow &flow) noexcept {
            if (stop) return; // Conditional task approach didn't work so just short-circuit here.
            f(flow, stat, d);
            flow.join();
            if (d % 16 == 0) err.invoke_noexcept([&] {throw_if_signal();});
            if (err.is_set()) {
                stat.store_failure(Stat::interrupt());
                stop = true;
            } else {
                stop = stat.should_stop();
            }
        });//.name("diag-" + std::to_string(d));

        if (!last.empty()) last.precede(next);
        last = std::move(next);
    }

    flow.join();
    return stat.load(end);
}

/******************************************************************************************/

uint diagonal_chunk_size(uint d) noexcept;

/******************************************************************************************/

template <class F>
void call_noexcept(SharedError &err, AtomicStat &stat, F &&run, uint i, uint j, int d) noexcept {
    try {
        if (unlikely(!run(i, j))) stat.store_failure(Stat(d));
    } catch (...) {
        stat.store_failure(Stat::bug());
        err.set_to_current_exception();
    }
}

// Call for_each_index but chunk the input range [b, e) into blocks of size n
template <class F>
void for_each_chunked(tf::Subflow &flow, uint const b, uint const e, uint const n, F &&f) noexcept {
    flow.for_each_index(b, e, n, [f=std::move(f), n, e](uint i) noexcept {
        for (auto j = i; j != min(i+n, e); ++j) f(j);
    });
}

/// Single strand top-level partition function iteration
// Takes index of first diagonal to do
// Returns Stat::finished() or the first diagonal that failed
template <class F, class G>
Stat iterate_single_strand(tf::Subflow &flow, SharedError &err, uint n, uint start_diag, F &&run, G &&next) noexcept {
    // Arrange diagonal by diagonal tasks
    return iterate_diagonal(flow, err, start_diag, n, [&next, &run, &err, n](auto &flow, auto &stat, uint d) noexcept {
        for_each_chunked(flow, 0, n-d, diagonal_chunk_size(d), [&err, &run, &stat, d](uint i) noexcept {
            call_noexcept(err, stat, run, i, i+d, d);
        });
        if (!stat.should_stop()) next();
    });
}

/// Multiple strand top-level partition function iteration
template <class F, class G>
Stat iterate_multi_strand(tf::Subflow &flow, SharedError &err, uint l, uint r, uint mid, uint const start_diag, uint const end_diag, F &&run, G &&next) noexcept {
    return iterate_diagonal(flow, err, start_diag, end_diag, [&next, &run, &err, l, r, mid](auto &flow, auto &stat, uint d) noexcept {
        for_each_chunked(flow, max(l-1, d) - l+1, min(r, d+1), diagonal_chunk_size(mid + d), [&err, &run, &stat, d, l](uint j) noexcept {
            call_noexcept(err, stat, run, l-1 + j - d, j, d);
        });
        if (!stat.should_stop()) if (!next()) stat.store_failure(Stat(d));
    });
}

/******************************************************************************************/

template <class Block, class Algebra>
struct BlockRef : Algebra {
    using value_type = typename Block::value_type;
    using block_type = Block;

    Block &block;

    BlockRef(Block &b, Ignore={}) : block(b) {}

    template <class R>
    static auto find(R recursion) {
        auto pred = [](auto x) {return std::is_base_of<R, std::decay_t<decltype(*x)>>();};
        return find_c(as_pack<decltype(members_of(std::declval<Block &>()))>(), pred);
    }

    template <class R>
    static constexpr bool contains(R recursion) {return !std::is_same_v<decltype(find(recursion)), not_found>;}

    template <class R, class ...Ts>
    decltype(auto) operator()(R recursion, Ts const &...ts) const {
        using N = decltype(find(recursion));
        static_assert(!std::is_same_v<N, not_found>, "Member of this type is not present");
        return Algebra::element(std::get<N::value>(members_of(static_cast<Block const &>(block))), ts...);
    }

    // template <class R, class I, class J, class K, class ...Ts>
    // decltype(auto) operator()(R recursion, I arg1, J arg2, K arg3, Ts const &...ts) const {
    //     using N = decltype(find(recursion));
    //     static_assert(!std::is_same_v<N, not_found>, "Member of this type is not present");
    //     auto blk = std::get<N::value>(members_of(static_cast<Block const &>(block)));
    //     if(typeid(blk).name().find("nupack::thermo::Table") != std::string::npos){
    //         return Algebra::element(blk.table[arg3], arg1, arg2, ts...);
    //     }
    //     else{
    //         return Algebra::element(blk, arg1, arg2, arg3, ts...);
    //     }
    // }
};


/******************************************************************************************/

template <class Block, class Algebra>
struct TagBlockRef : BlockRef<Block, Algebra> {
    using base_type = BlockRef<Block, Algebra>;
    int strand;

    TagBlockRef(Block &b, int s) : base_type(b), strand(s) {}

    template <class R, class ...Ts>
    decltype(auto) operator()(R recursion, Ts const &...ts) const {return base_type::operator()(recursion, strand, ts...);}
};

/******************************************************************************************/

template <int N, class Block, class Model, class Algebra>
struct SingleConstants : BlockRef<Block, Algebra> {
    using BlockRef<Block, Algebra>::block;
    using complexity = size_constant<N>;
    using model_type = Model;
    using is_multi = False;

    Sequence const sequence;
    Model const &model;
    Action action;
    typename Block::StrandData const &linear;

    SingleConstants(Block &b, Sequence s, Model const &mod, Action a, typename Block::StrandData const &l)
        : BlockRef<Block, Algebra>(b), sequence(std::move(s)), model(mod), action(std::move(a)), linear(l) {}

    SingleStrand multi() const {return {};}

    ComplexView strands() const {return {&sequence, &sequence + 1};}
    constexpr std::size_t n_strands() const {return 1;}

    template <class E>
    auto pairing(int i, int j, E &&expression) const noexcept {
        bool can_pair = model.can_pair(false, next(sequence, i), next(sequence, j));
        return action(i, j, can_pair, Algebra(), block, strands(), model, static_cast<E &&>(expression));
    }

    template <class Matrix>
    auto left_plus(Matrix const &b, Ignore, Index i, Index j) const {NUPACK_QUICK_REQUIRE(i, <, j); return b(i+1, j);}

    template <class Matrix>
    auto right_minus(Matrix const &b, Ignore, Index i, Index j) const {NUPACK_QUICK_REQUIRE(i, <, j); return b(i, j-1);}

    auto left(Index i) const {NUPACK_QUICK_REQUIRE(i, >=, 0); NUPACK_QUICK_REQUIRE(i, <, lsize()); return sequence[i];}
    auto left(Index i, Index j) const {return view(sequence, i, j);}
    auto lsize() const {return sequence.size();}

    constexpr bool on_left(Index i, Index x) const {return x > 0 ? i + x < lsize() : i + x >= 0;}
    constexpr bool on_right(Index i, Index x) const {return x > 0 ? i + x < rsize() : i + x >= 0;}

    auto right(Index i) const {NUPACK_QUICK_REQUIRE(i, >=, 0); NUPACK_QUICK_REQUIRE(i, <, lsize()); return sequence[i];}
    auto right(Index i, Index j) const {return view(sequence, i, j);}
    auto rsize() const {return sequence.size();}

    Stat forward_pass(tf::Subflow &flow, SharedError &err, uint start, Ignore, Ignore) const noexcept {
        return iterate_single_strand(flow, err, lsize(), start,
            [&](auto i, auto j) {
                bool ok = true;
                for_each(members_of(block), [&](auto &m) {
                    ok = ok && m.calculate(SingleStrand(), i, j, *this);
                });
                return ok;
            },
            [&]() noexcept {
                for_each(members_of(block), [&](auto &m) noexcept {m.next_diagonal();});
            }
        );
    }

    template <class F>
    [[nodiscard]] auto sum_bases(F &&f) const {return this->total(model.energy_model.alphabet().all(), std::forward<F>(f));}
};

/******************************************************************************************/

template <class Block, class Model, class Algebra>
struct OpConstants : BlockRef<Block const, Algebra> {
    Sequence left, right;
    Model const &model;
    uint n_strands;
    // typename Block::StrandData const &linear;

    NUPACK_REFLECT(OpConstants, left, right, n_strands);

    OpConstants(Block const &b, Model const &mod, Sequence s1, Sequence s2, uint n_strands, Algebra={})
        : BlockRef<Block const, Algebra>(b), left(std::move(s1)), right(std::move(s2)), model(mod), n_strands(n_strands) {}

    // bool is_first(Index i) const;
    bool included(Index j) const {return j < len(right);}

    // auto const n_strands() const {return n_strands;}
    bool single() const {return n_strands == 1;}
    bool adjacent(Index i, Index j) const {
        return (n_strands == 1 && i+1 == j)
            || (n_strands == 2 && i+1 == len(left) && j == 0);
    }
};

template <class Block, class Model, class Algebra>
OpConstants(Block const &, Model const &, Sequence, Sequence, uint, Algebra) -> OpConstants<Block, Model, Algebra>;

/******************************************************************************************/

// Example of the contents of "middle"
    // 01 -> (0,0)  NULL
    //        NULL (1,1)
    // 012 -> (0,0) (0,1) NULL
    //        NULL  (1,2) (2,2)
    // 0123 -> (0,0) (0,1) (0,2) NULL
    //         NULL  (1,3) (2,3) (3,3)
template <int N, class Block, class Model, class Algebra>
struct MultiConstants : BlockRef<Block, Algebra> {
    using complexity = size_constant<N>;
    using base_type = BlockRef<Block, Algebra>;
    using model_type = Model;
    using Linear = typename Block::StrandData const;
    using Ref = if_t<Algebra::is_forward::value, BlockRef<Block const, Algebra>, TagBlockRef<Block const, Algebra>>;
    using is_multi = True;

    // There are N middle blocks for N sequences
    // The first one has null left, the last one has null right, all others non-null
    struct MiddleBlock {
        Block const *left;
        Block const *right;
        Linear const *linear;
        int index;
        uint last;
        NUPACK_REFLECT(MiddleBlock, left, right, linear, index, last);
    };

    vec<Guard *> locks;

    // as many of these as there are strands
    vec<MiddleBlock> middle;
    Sequence const left_sequence, right_sequence;
    Model const &model;

    SequenceList strand_list;
    Action action;
    using base_type::block;

    ComplexView strands() const {return {strand_list.data(), strand_list.data() + strand_list.size()};}
    auto n_strands() const {return middle.size();}

    template <class ...Ts>
    decltype(auto) left_linear(Ts &&...ts) const {return (*front(middle).linear)(static_cast<Ts &&>(ts)...);}

    decltype(auto) left_linear() const {return *front(middle).linear;}

    template <class ...Ts>
    decltype(auto) right_linear(Ts &&...ts) const {return (*back(middle).linear)(static_cast<Ts &&>(ts)...);}

    decltype(auto) right_linear() const {return *back(middle).linear;}

    template <class V1, class V2>
    MultiConstants(Block &b, ConstView<SequenceList> const &v, Model const &mod, Action a, V1 const &blocks, V2 const &strands)
        : base_type(b), left_sequence(front(v)), right_sequence(back(v)), model(mod), strand_list(v.begin(), v.end()), action(std::move(a)) {
            auto const n = std::size(v);
            NUPACK_REQUIRE(n, >=, 2, v);
            NUPACK_REQUIRE(len(blocks), ==, 2 * (n-1), v);
            locks.reserve(len(blocks) + len(strands));
            extend(locks, item_view(strands));
            extend(locks, item_view(blocks));

            for (auto const &x : blocks) NUPACK_QUICK_ASSERT(x.first);
            // The following order is assumed to match that in the block map construction
            middle.push_back({blocks.end()[-2].first, nullptr, strands[0].first, 0, uint(front(v).size() - 1u)});
            for (auto i : range<int>(1, n-1)) {
                NUPACK_ASSERT(!v[i].empty(), "Empty sequences are not allowed");
                middle.push_back({blocks.end()[-2*i-2].first, blocks[2*i-1].first, strands[i].first, int(i), uint(v[i].size() - 1u)});
            }
            middle.push_back({nullptr, back(blocks).first, strands[n-1].first, int(n-1), uint(back(v).size() - 1u)});
            NUPACK_QUICK_REQUIRE(middle.size(), ==, n);
        }

    MultiStrand multi() const {return {};}

    template <class E>
    auto pairing(Index i, Index j, E &&expression) const noexcept {
        bool can_pair = model.can_pair(true, next(left_sequence, i), next(right_sequence, j));
        return action(i, j, can_pair, Algebra(), *this, strands(), model, static_cast<E &&>(expression));
    }

    template <bool Right>
    static Ref ref(MiddleBlock const &m) {
        if constexpr(Right) {NUPACK_QUICK_ASSERT(m.right); return {*m.right, -m.index};}
        else {NUPACK_QUICK_ASSERT(m.left); return {*m.left, m.index};} // 0 is unambiguous: it must be on the left
    }

    template <class F>
    auto first_middle(F &&f) const {auto const &m = middle[1]; return f(ref<1>(m), *m.linear, uint(0));}

    template <class F>
    auto last_middle(F &&f) const {auto const &m = middle.end()[-2]; return f(ref<0>(m), *m.linear, m.last);}

    auto middle_at(Index i) const {
        auto const &m = middle[i];
        return std::make_tuple(Ref{m.left ? *m.left : this->block},
                               Ref{m.right ? *m.right : this->block}, m.last);
    }

    auto linear_at(Index i) const {
        auto const &m = middle[i];
        return *m.linear;
    }

    auto middle_at_first() const{
        return ref<0>(front(middle));
    }

    auto middle_at_last() const{
        return ref<1>(back(middle));
    }

    auto left(Index i) const {NUPACK_QUICK_REQUIRE(i, >=, 0); NUPACK_QUICK_REQUIRE(i, <, lsize()); return left_sequence[i];}
    auto right(Index i) const {NUPACK_QUICK_REQUIRE(i, >=, 0); NUPACK_QUICK_REQUIRE(i, <, rsize()); return right_sequence[i];}
    auto left(Index i, Index j) const {return view(left_sequence, i, j);}
    auto right(Index i, Index j) const {return view(right_sequence, i, j);}
    bool adjacent(Index i, Index j) const {return j == 0 && n_strands() == 2 && i + 1 == lsize();}

    constexpr bool on_left(Index i, Index x) const {return x > 0 ? i + x < lsize() : i + x >= 0;}
    constexpr bool on_right(Index i, Index x) const {return x > 0 ? i + x < rsize() : i + x >= 0;}

    template <class Matrix, class R>
    decltype(auto) left_plus(Matrix const &b, R, Index i, Index j) const {
        return i+1 == lsize() ? ref<1>(middle[1])(R(), 0, j) : b(i+1, j);
    }
    template <class Matrix, class R>
    decltype(auto) right_minus(Matrix const &b, R, Index i, Index j) const {
        return j ? b(i, j-1) : ref<0>(middle.end()[-2])(R(), i, middle.end()[-2].last);
    }

    // Go over all k, k+1 bases which are on the same strand, k ∈ [i, j-1]
    template <class F>
    auto sandwich(Index i, Index j, F &&f) const {
        NUPACK_QUICK_REQUIRE(i, >=, 0);
        NUPACK_QUICK_REQUIRE(j, >=, 0);
        NUPACK_QUICK_REQUIRE(i, <, lsize());
        NUPACK_QUICK_REQUIRE(j, <, rsize());
        base_type const &blk = *this;
        return Algebra::sum(
            // left strand
            on_left(i, +1) ? Algebra::nullable(f(ref<0>(front(middle)), blk, span(i, lsize() - 1), *front(middle).linear)) : Algebra::zero(),
            // right strand
            on_right(j, -1) ? Algebra::nullable(f(blk, ref<1>(back(middle)), span(0, j), *back(middle).linear)) : Algebra::zero(),
            // middle strands
            Algebra::total(view(middle).offset(1, -1), [&, f](auto const &m) {
                // second axis of left should match first axis of right (checked in constructor)
                return m.last ? Algebra::nullable(f(ref<0>(m), ref<1>(m), span(0, m.last), *m.linear)) : Algebra::zero();
            })
        );
    }

    template <class F>
    auto nick_total(F &&f) const {
        base_type const &blk = *this;
        return Algebra::total(range(n_strands() - 1),
            [&, f](int i) {return f(ref<0>(middle[i]), ref<1>(middle[i+1]), *middle[i].linear, *middle[i+1].linear, middle[i].last);}
        );
    }

    auto lsize() const {return left_sequence.size();}
    auto rsize() const {return right_sequence.size();}

    /******************************************************************************************/

    Stat forward_pass(tf::Subflow &flow, SharedError &err, uint start, uint end, bool stage) const noexcept {
        return iterate_multi_strand(flow, err,
            lsize(), rsize(), sum(view(strand_list).offset(1, -1), len), start, end,
            [&](auto i, auto j) {
                // everything inside here invoked in serial!
                bool ok = true;
                Guard::TryReadLock lock(locks, ok, stage);
                for_each(members_of(block), [&](auto &m) {
                    ok = ok && m.calculate(MultiStrand(), i, j, *this);
                });
                return ok;
            },
            [&]() noexcept {
                bool ok = true;
                Guard::TryReadLock lock(locks, ok, stage);
                if (ok) for_each(members_of(block), [&](auto &m) noexcept {m.next_diagonal();});
                return ok;
            }
        );
    }

    template <class F>
    [[nodiscard]] auto sum_bases(F &&f) const {return this->total(model.energy_model.alphabet().all(), std::forward<F>(f));}
};

/******************************************************************************************/

}
