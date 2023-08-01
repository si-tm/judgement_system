#include <nupack/thermo/Schedule.h>
#include <nupack/thermo/Cache.h>
#include <nupack/standard/Variant.h>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>

namespace nupack {

/******************************************************************************************/

void SharedError::set_to_current_exception() noexcept {
    if (!ptr->set.exchange(true, std::memory_order_relaxed)) {
        ptr->exception = std::current_exception();
        ptr->stored.store(true, std::memory_order_relaxed);
    }
}

bool SharedError::clear() noexcept {
    bool out = ptr->stored.exchange(false, std::memory_order_relaxed);
    ptr->exception = nullptr;
    ptr->set.exchange(false, std::memory_order_relaxed);
    return out;
}

void SharedError::rethrow_if_set() const {
    if (is_set()) {
        // this should practically never matter ... but theoretically it could
        while (!ptr->stored.load(std::memory_order_relaxed)) {}
        std::rethrow_exception(ptr->exception);
    }
}

bool SharedError::is_set() const noexcept {
    return ptr->set.load(std::memory_order_relaxed);
}

}

namespace nupack::thermo {

/******************************************************************************************/

std::atomic<std::size_t> AllocatedBytes{0};

/******************************************************************************************/

std::size_t get_hash(BlockData const &d) {return std::hash<Complex>()(d.complex);}

std::size_t get_hash(StrandData const &d) {return std::hash<Sequence>()(d.sequence);}

struct CacheImpl final : Cache {
    struct Value {
        std::variant<std::shared_ptr<BlockData>, std::shared_ptr<StrandData>> data;
        std::size_t memory = 0, time;

        NUPACK_REFLECT(Value, data, memory, time);

        Value(std::size_t t, std::size_t mem, Sequence const &s) : data(std::make_shared<StrandData>(s)), memory(mem), time(t) {}
        Value(std::size_t t, std::size_t mem, Complex const &x, Action const &a) : data(std::make_shared<BlockData>(x, a)), memory(mem), time(t) {}
    };

    struct GetKey {
        using result_type = Value const &;
        result_type operator()(Value const &v) const {return v;}
    };

    struct GetTime {
        using result_type = std::size_t;
        result_type operator()(Value const &v) const {return v.time;}
    };

    struct Hash {
        std::size_t operator()(Value const &v) const {
            std::size_t h = v.data.index();
            boost::hash_combine(h, fork(v.data, [](auto const &p) {return get_hash(*p);}));
            return h;
        }
    };

    struct Equal {
        bool operator()(Value const &v, Value const &o) const {
            return v.data == o.data;
        }
    };

    using Map = boost::multi_index_container<
        Value,
        boost::multi_index::indexed_by<
            boost::multi_index::hashed_unique<GetKey, Hash, Equal>,
            boost::multi_index::ordered_non_unique<GetTime>
        >
    >;

    std::mutex mut;
    std::shared_ptr<ModelData> model;
    Map map;
    std::size_t time = 0, mem = 0;

    explicit CacheImpl() : model(std::make_shared<ModelData>()) {}

    std::size_t used() const override {return mem;}
    void clear() override {map.clear();}

    auto emplace_strand(std::size_t mem, Sequence const &s) {
        auto pair = map.emplace(time, mem, s);
        if (!pair.second) update_time(pair.first);
        return pair;
    } 

    auto emplace_complex(std::size_t mem, Complex const &s, Action const &action) {
        auto pair = map.emplace(time, mem, s, action);
        if (!pair.second) update_time(pair.first);
        return pair;
    }

    template <class Iter>
    void update_time(Iter &it) {NUPACK_ASSERT(map.modify(it, [&](Value &v) {v.time = time;}));}
    
    std::size_t shrink_to_limit(std::size_t limit) override {
        if (mem > limit) {
            std::lock_guard<std::mutex> lk(mut);
            while (mem > limit) {
                NUPACK_ASSERT(!map.empty());
                auto it = map.get<1>().begin();
                mem -= it->memory;
                map.get<1>().erase(it);
            }
        }
        return mem;
    }
};

std::shared_ptr<Cache> Cache::create() {return std::make_shared<CacheImpl>();}

/******************************************************************************************/

template <class F>
void for_each_dependency(Complex const &key, F &&f) {
    for (auto n : range(len(key)-1)) {
        f(view(key, 0, len(key)-n-1));
        f(view(key, n+1, len(key)));
    }
}

struct WorkMap {
    struct Block {
        std::shared_ptr<bool> full_requested;
        std::optional<tf::Task> task;
        decltype(std::declval<CacheImpl>().map.emplace(std::shared_ptr<BlockData>(), 0).first) iterator;

        Block() : full_requested(std::make_shared<bool>(false)) {}
    };

    struct Strand {
        std::optional<tf::Task> task;
        decltype(std::declval<CacheImpl>().map.emplace(std::shared_ptr<StrandData>(), 0).first) iterator;
    };

    std::unordered_map<Sequence, Strand> strands;
    std::unordered_map<Complex, Block> blocks;
    std::size_t const start_time;
    std::size_t memory;
    vec<tf::Task> root_tasks;

    WorkMap(std::size_t t, std::size_t mem) : start_time(t), memory(mem) {}

    /**************************************************************************************/

    tf::Task request_root(tf::Taskflow &flow, CacheImpl const &cache, Root const &root, Callback c, SharedError const &error) {
        vec<tf::Task> tasks;
        Dependencies deps;
        for (auto const &[k, f] : root.blocks) {
            auto &a = blocks.at(k);
            deps.blocks.emplace_back(std::get<0>(a.iterator->data));
            if (a.task) tasks.emplace_back(*a.task);
        }
        for (auto const &s : root.strands) {
            auto &a = strands.at(s);
            deps.strands.emplace_back(std::get<1>(a.iterator->data));
            if (a.task) tasks.emplace_back(*a.task);
        }

        auto task = flow.emplace([err=error, c=std::move(c), deps=std::move(deps), model=cache.model](tf::Subflow &sub) mutable {
            // static std::atomic<int> count{0};
            // print("compute root", ++count, AllocatedBytes.load(), tf::node_pool.capacity());
            if (!err.is_set()) err.invoke_noexcept(c, sub, std::move(deps), *model);
            deps.reset();
            sub.reset();
        });
        for (auto &t : tasks) t.precede(task);
        return task;
    }

    void shrink_to_limit_no_lock(CacheImpl &cache, vec<tf::Task> &new_tasks, std::size_t max) {
        if (max != 0) while (memory > max) {
            auto iter = cache.map.get<1>().begin();
            if (iter->time == cache.time) break;
            if (iter->time >= start_time) // make sure the old root is done
                for (auto &t : new_tasks) root_tasks.at(iter->time - start_time).precede(t);

            if (iter->data.index()) strands.erase(std::get<1>(iter->data)->sequence);
            else blocks.erase(std::get<0>(iter->data)->complex);

            memory -= iter->memory;
            cache.mem -= iter->memory;
            cache.map.get<1>().erase(iter);
        }
    }
};

/******************************************************************************************/

// Compare complexes by their number of strands
// This is used so that blocks may be iterated in order of their logical dependency
struct CompareByLength {
    bool operator()(Complex const &a, Complex const &b) const {
        return std::forward_as_tuple(len(a), a) > std::forward_as_tuple(len(b), b);
    }
};

/******************************************************************************************/

Schedule schedule_all(tf::Taskflow &flow, Computer const &computer, ComputeOptions const &ops, CallbackList callbacks) {
    // print("computing", callbacks, bool(computer.cache), bool(ops.cache_complexes));
    Schedule out;
    if (callbacks.empty()) {
        out.start = flow.emplace([]{});
        out.finish = flow.emplace([]{}).succeed(out.start);
        return out;
    }
    NUPACK_ASSERT(computer.cache);

    // Lock while computation of the cache elements is in progress and while this schedule is being created
    out.start = flow.emplace([cache=computer.cache]{
        dynamic_cast<CacheImpl &>(*cache).mut.lock();
    });
    out.finish = flow.emplace([cache=computer.cache]{
        dynamic_cast<CacheImpl &>(*cache).mut.unlock();
    });
    auto &cache = dynamic_cast<CacheImpl &>(*computer.cache);
    std::lock_guard<std::mutex> lock(cache.mut);

    std::optional<CacheImpl> fake_cache;
    CacheImpl &complex_cache = ops.cache_complexes ? cache : fake_cache.emplace();
    std::size_t max_memory = 0;

    for (auto const &[r, c] : callbacks) for (auto const &[k, f] : r.blocks)
        for (auto const &s : k) max_eq(cache.model->capacity, 2 * len(s));

    auto prepare = computer.prepare(flow, cache.model, ops.error).succeed(out.start);

    WorkMap map(cache.time, cache.mem);
    map.root_tasks.reserve(len(callbacks));

    for (auto const &[root, c] : callbacks) {        
        vec<tf::Task> new_tasks;

        auto put_strand = [&](Sequence const &s) {
            auto [it, inserted] = map.strands.try_emplace(s);
            if (!inserted) {cache.update_time(it->second.iterator); return;}
            std::size_t const mem = computer.strand_memory(s);
            auto [iter, uncached] = cache.emplace_strand(mem, s);
            it->second.iterator = iter;
            if (uncached) {
                it->second.task.emplace(flow.emplace([err=ops.error, c=computer.compute_strand, model=cache.model, p=std::get<1>(iter->data)]() mutable {
                    // print("compute strand", p->sequence, err.is_set());
                    if (!err.is_set()) err.invoke_noexcept(c, *p, *model);
                    p.reset();
                })).succeed(prepare);
                
                map.memory += mem;
                cache.mem += mem;
                new_tasks.emplace_back(*it->second.task);
            } 
        };
        for (auto const &s : root.strands) put_strand(s);

        std::map<Complex, bool, CompareByLength> blocks;

        auto reserve_complex = [&](Complex const &key, bool const full) {
            auto [it, inserted] = map.blocks.try_emplace(key);
            *it->second.full_requested |= full;

            if (inserted) {
                std::size_t const mem = computer.block_memory(key);
                auto [iter, uncached] = complex_cache.emplace_complex(mem, key, ops.action);
                it->second.iterator = iter;
            } else {
                complex_cache.update_time(it->second.iterator);
            }

            std::size_t required_diagonal = (full && len(key) > 1) ? len(key.front()) + len(key.back()) - 1 : len(key[0]);
            if (it->second.task || (std::get<0>(it->second.iterator->data)->stat.value >= required_diagonal)) {
                return false;
            } else {
                for_each_dependency(key, [&](auto const &k) {
                    blocks.try_emplace(k).first->second = true;
                });
                return true;
            }
        };

        for (auto const &[x, full] : root.blocks) // iterate from larger to smaller blocks
            blocks.try_emplace(x, false).first->second |= full;

        for (auto &[x, full] : blocks) {
            full = reserve_complex(x, full);
            for (auto const &s : x) put_strand(s); 
        }

        auto put_complex_task = [&](Complex const &key) {
            auto &x = map.blocks.at(key);
            vec<tf::Task> tasks;
            Dependencies deps;

            for_each_dependency(key, [&](auto const &k) {
                auto &a = map.blocks.at(k);
                deps.blocks.emplace_back(std::get<0>(a.iterator->data));
                if (a.task) tasks.emplace_back(*a.task);
            });

            for (auto const &s : key) {
                auto &a = map.strands.at(s);
                deps.strands.emplace_back(std::get<1>(a.iterator->data));
                if (a.task) tasks.emplace_back(*a.task);
            }

            x.task.emplace(flow.emplace([err=ops.error, c=computer.compute_block, model=cache.model, p=std::get<0>(x.iterator->data), deps=std::move(deps), full=x.full_requested](tf::Subflow &flow) mutable {               
                // print("compute block", p->complex, err.is_set());
                if (!err.is_set()) err.invoke_noexcept(c, flow, *p, err, *model, std::move(deps), *full);
                deps.reset();
                p.reset();
                flow.reset();
            }).succeed(prepare));

            map.memory += x.iterator->memory;
            complex_cache.mem += x.iterator->memory;
            new_tasks.emplace_back(*x.task);

            for (auto &t : tasks) t.precede(*x.task);
        };

        for (auto const &[k, needs_compute] : view(blocks.rbegin(), blocks.rend())) // iterate from smaller to larger blocks
            if (needs_compute) put_complex_task(k);

        map.root_tasks.emplace_back(map.request_root(flow, cache, root, c, ops.error).succeed(prepare));
        new_tasks.emplace_back(map.root_tasks.back());

        map.shrink_to_limit_no_lock(cache, new_tasks, ops.max_bytes);
        max_eq(max_memory, map.memory);

        ++cache.time;
    }

    for (auto &t : map.root_tasks) out.finish.succeed(t);
    out.root_tasks = std::move(map.root_tasks);
    return out;
}

/******************************************************************************************/

}
