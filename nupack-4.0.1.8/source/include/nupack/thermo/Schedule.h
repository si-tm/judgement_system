#pragma once
#include <atomic>
#include <nupack/standard/Map.h>
#include <nupack/types/Sequence.h>
#include <nupack/thermo/Action.h>
#include <nupack/thermo/Common.h>
#include <nupack/thermo/Cache.h>
#include <nupack/execution/Executor.h>
#include "../common/Threading.h"

namespace nupack::thermo {

/******************************************************************************************/

struct Resource {
    real cost=0;
    std::size_t memory=0;
    NUPACK_REFLECT(Resource, cost, memory);
};

/******************************************************************************************/

struct Root {
    vec<std::pair<Complex, bool>> blocks;
    vec<Sequence> strands;

    NUPACK_REFLECT(Root, blocks, strands);
};

/******************************************************************************************/

struct Item {
    struct Concept {
        virtual std::unique_ptr<Concept> clone() const = 0;
        virtual std::string describe() const = 0;
        virtual ~Concept() {}
    };

    template <class T>
    struct Impl final : Concept {
        T value;
        template <class ...Ts>
        Impl(Ts &&...ts) : value(fw<Ts>(ts)...) {}
        std::unique_ptr<Concept> clone() const override {return std::make_unique<Impl>(value);}

        std::string describe() const override {
            std::stringstream ss;
            dump_os(ss, TypeName<T>(), " ", value);
            return ss.str();
        }
    };

    std::unique_ptr<Concept> ptr;
    NUPACK_REFLECT(Item, ptr);

    Item() = default;
    Item(Item &&) noexcept = default;
    Item& operator=(Item &&) noexcept = default;

    Item(Item const &i) : ptr(i.ptr ? i.ptr->clone() : std::unique_ptr<Concept>()) {}
    Item& operator=(Item const &i) {if (i.ptr) ptr = i.ptr->clone(); else ptr.reset(); return *this;}

    std::string describe() const {return ptr ? ptr->describe() : "null";}

    template <class T, class ...Ts>
    static Item make(Ts &&...ts) {
        Item i;
        i.ptr = std::make_unique<Impl<T>>(fw<Ts>(ts)...);
        return i;
    }

    template <class T>
    static Item from(T &&t) {return make<std::decay_t<T>>(fw<T>(t));}

    template <class T, class ...Ts>
    T & emplace(Ts &&...ts) {
        auto p = std::make_unique<Impl<T>>(fw<Ts>(ts)...);
        T &o = p->value;
        ptr = std::move(p);
        return o;
    }

    template <class T>
    T *get() {
        static_assert(std::is_same_v<T, std::decay_t<T>>);
        auto p = dynamic_cast<Impl<T> *>(ptr.get());
        return p ? &p->value : nullptr;
    }
    
    template <class T>
    T const *get() const {return const_cast<Item &>(*this).get<T>();}

    template <class T>
    T const &ref() const {auto t = get<T>(); NUPACK_QUICK_ASSERT(t, "unmatched data type", bool(ptr), describe(), TypeName<T>()); return *t;}
    
    template <class T>
    T &ref() {auto t = get<T>(); NUPACK_QUICK_ASSERT(t, "unmatched data type", bool(ptr), describe(), TypeName<T>()); return *t;}
    
    template <class T>
    T take() {return std::move(ref<T>());}

    explicit operator bool() const {return bool(ptr);}
};

/******************************************************************************************/

// The way taskflow interacts with read/write locking is really annoying
// You would think for translating blocks that a shared_mutex would work
// but actually it deadlocks commonly because of how weirdly the control flow goes inside taskflow
// The only solution is to only read/write lock while actively touching the block, i.e. at a level beneath which there are no subflows.
// I did this with an atomic variable as below. Spinning only rarely occurs, so it should be pretty efficient.
// Odd/even is used in the atomic to signal what stage the block is in (normal or backup)
// Multiples of 2 are added for each active reader of the block
// A large negative number is used for the write lock to make sure other threads don't read while it is being written

struct Guard {
    std::atomic<int> val{0};

    // return whether a read lock for the given stage could be acquired
    [[nodiscard]] bool read_lock(bool stage) noexcept {
        auto v = val.fetch_add(2, std::memory_order_acquire);
        return v >= 0 && (v % 2 == stage);
    }

    int value() const {return val.load();}

    void set_stage(bool s) {val.store(int(s));}

    void read_unlock() noexcept {
        val.fetch_sub(2, std::memory_order_release);
    }
    
    // Read lock which spins until the block is in a stable state (not being upgraded)
    class ReadLock {
        Guard &guard;
    public:
        ReadLock(Guard &g) : guard(g) {
            auto v = guard.val.fetch_add(2, std::memory_order_acquire);
            while (v < 0) v = guard.val.load(std::memory_order_relaxed); // Wait until it's not being written
        }
        
        ~ReadLock() {guard.val.fetch_sub(2, std::memory_order_release);}
    };

    // Upgrade lock which spins until the block can be put safely into the upgraded form
    class UpgradeLock {
        Guard &guard;
        bool acquired;
        static constexpr int big = 100000; // this should be safely bigger than the number of threads
    public:
        UpgradeLock(Guard &g) : guard(g) {
            // Spin until either upgrade has been performed by another thread or this thread can do the upgrade
            // This is the only operation that can move val to an odd number, and it can only occur once
            while (true) {
                int v = 0;
                acquired = guard.val.compare_exchange_weak(v, 1-big, std::memory_order_acquire);
                if (acquired || (v >= 0 && (v % 2 == 1))) return;
            }
        }
        // Whether the write lock has been acquired
        explicit operator bool() const noexcept {return acquired;}
        // Wait until another thread finishes the upgrade
        void wait() const {while (guard.val.load(std::memory_order_relaxed) < 0) {}}

        ~UpgradeLock() {
            if (acquired) guard.val.fetch_add(big, std::memory_order_release);
        }
    };

    // This is the read lock type used at the lowest level of the dynamic programs
    // Basically it tries to get a read lock for the given stage, and if one isn't available, the control flow
    // should exit because the running program will not be able to run safely
    struct TryReadLock {
        vec<Guard *> const &locks;

        TryReadLock(vec<Guard *> const &v, bool &ok, bool stage) noexcept : locks(v) {
            for (auto l : locks) ok &= l->read_lock(stage);
        }

        ~TryReadLock() {for (auto l : locks) l->read_unlock();}
    };

};


/******************************************************************************************/

struct BlockData {
    Guard guard;
    Item data;
    Complex complex;
    Action action;
    Stat stat = Stat(0);

    BlockData(Complex x, Action a) : complex(std::move(x)), action(std::move(a)) {}

    NUPACK_REFLECT(BlockData, guard, data, complex, action, stat);
};

/******************************************************************************************/

struct StrandData {
    Guard guard;
    Item data;
    Sequence sequence;

    StrandData(Sequence s) : sequence(std::move(s)) {}

    NUPACK_REFLECT(StrandData, guard, data, sequence);
};

/******************************************************************************************/

struct Dependencies {
    vec<std::shared_ptr<BlockData>> blocks;
    vec<std::shared_ptr<StrandData>> strands;

    void reset() {blocks.clear(); strands.clear();}

    NUPACK_REFLECT(Dependencies, blocks, strands);
};

/******************************************************************************************/

// Here the normal and backup models are just kept separately
// Pretty sure a normal shared_mutex is fine because the calculation happens once and at the top level
struct ModelData {
    Item data, backup;
    std::shared_mutex mutex;
    std::size_t capacity = 0;
};

/******************************************************************************************/

using Callback = std::function<void(tf::Subflow &, Dependencies, ModelData &)>;
using CallbackList = vec<std::pair<Root, Callback>>;

/******************************************************************************************/

struct Computer {
    std::shared_ptr<Cache> cache;
    
    std::function<tf::Task(tf::Taskflow &, std::shared_ptr<ModelData>, SharedError)> prepare;
    std::function<void(StrandData &, ModelData &)> compute_strand;
    std::function<void(tf::Subflow &, BlockData &, SharedError &, ModelData &, Dependencies const &, bool)> compute_block;
    std::function<std::size_t(Complex const &)> block_memory;
    std::function<std::size_t(Sequence const &)> strand_memory;
    
    using Extension = std::function<void(CallbackList &, Item, std::function<void(Item)>)>;
    vec<Extension> extensions;    
};

/******************************************************************************************/

struct ComputeOptions {
    Action action;
    SharedError error;
    SharedExecutor executor;
    std::size_t max_bytes = 0;
    bool cache_complexes = true;

    NUPACK_REFLECT(ComputeOptions, action, error, executor, max_bytes, cache_complexes);
};

/******************************************************************************************/

struct Schedule {
    vec<tf::Task> root_tasks;
    tf::Task start, finish;
    std::size_t max_bytes = 0;
};

Schedule schedule_all(tf::Taskflow &flow, Computer const &computer, ComputeOptions const &ops, CallbackList callbacks);

/******************************************************************************************/

template <class Ensemble, class T>
struct MemoryFactor;

}
