#pragma once
#include "Operations.h"
#include "Executor.h"
#include "../standard/Vec.h"
#include "../iteration/Range.h"

#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/partitioner.h>
#include <oneapi/tbb/task_arena.h>
#include <oneapi/tbb.h>
#include <tbb/cache_aligned_allocator.h>
#include <thread>
#include <atomic>

namespace nupack {

/******************************************************************************************/

inline auto this_thread_id() {return std::this_thread::get_id();}

/******************************************************************************************/

inline auto default_thread_number() {return std::min<int>(TotalCPU, std::thread::hardware_concurrency());}

/******************************************************************************************/

struct SharedImpl {
    struct State {
        using simple_type = True;
        State() = default;
        State(usize n) : arena(std::make_unique<oneapi::tbb::task_arena>(n)), max(n) {}

        std::unique_ptr<oneapi::tbb::task_arena> arena;
        std::mutex mut;
        usize max;

        auto save_repr() const {return make_members(max);}
        void load_repr(usize m) {max = m; arena = std::make_unique<oneapi::tbb::task_arena>(max);}
    };

    SharedImpl(usize n=0) : state(std::make_shared<State>(n ? n : default_thread_number())) {}

    template <class V, class F>
    void implement_map(V const &v, GrainSize g, F const &f, OrderedSplit) const {
        std::atomic<usize> count{0u};
        SharedError err;
        state->arena->execute([&] {tbb::parallel_for(usize(0), usize(len(v)), g.value, [&] (usize i) {
            if (!err.is_set()) err.invoke_noexcept([&] {f(v[count++]);});
        }, tbb::auto_partitioner());});
        err.rethrow_if_set();
    }

    /// Parallel for with manually specified granularity
    template <class V, class F>
    void implement_map(V const &v, GrainSize g, F const &f, AffinitySplit) const {
        SharedError err;
        state->arena->execute([&] {tbb::parallel_for(tbb::blocked_range<usize>(0u, len(v), g.value), [&](auto const &b) {
            if (!err.is_set()) err.invoke_noexcept([&] {for (auto i : iterators(b)) f(v[i]);});
        }, tbb::auto_partitioner());});
        err.rethrow_if_set();
    }

    // template <class R, class V>
    // auto reduce(V const &v, R const &r) const {
    //     auto range = tbb::blocked_range<const_iterator_of<V>>(begin_of(v), end_of(v));
    //     auto const acc = DefaultAccumulator();
    //     return tbb::parallel_reduce(range, value_type_of<V>(),
    //         [&r, &acc](auto const &a, auto const &o) {return acc(begin_of(a), end_of(a), o, r);}, r);
    // }

    auto n_workers() const {return state->max;}

    NUPACK_REFLECT(SharedImpl, state);

private:
    std::shared_ptr<State> state;
};

/******************************************************************************************/

}
