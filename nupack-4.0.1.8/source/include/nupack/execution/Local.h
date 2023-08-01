/**
 * @brief Definition of Local executor for shared-memory parallelism
 *
 * @file Local.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "Operations.h"

#include "Shared.h"
#include "Serial.h"
#include "../standard/Vec.h"
#include "../standard/Variant.h"

namespace nupack {

/******************************************************************************************/

/// Default constructor Local() makes a serial executor
struct Local {
    Variant<SerialImpl, SharedImpl> executor;
    Local(uint n=1) : executor(n == 1 ? SerialImpl() : Variant<SerialImpl, SharedImpl>(SharedImpl(n))) {}

    NUPACK_REFLECT(Local, executor);

    auto n_workers() const {return fork(executor, [](auto const &ex) -> uint {return ex.n_workers();});}

    /**
     * @brief Parallelize a functor across v which overwrites v
     * @param v: a mutable container, should generally be random access
     * @param g: minimum number of elements to run in a task
     * @param f: functor returning new_element from (*this, element, index)
     * @param tag: execution tag
     */
    template <class V, class F, class T=OrderedSplit>
    auto map(V const &range, F const &fun, GrainSize g=1, T tag={}) const {
        using R = std::decay_t<decltype(fun(*std::begin(range)))>;
        if constexpr(std::is_same_v<void, R>) {
            fork(executor, [&](auto const &ex) {ex.implement_map(range, g, fun, tag);});
        } else {
            vec<R> out(len(range));
            map(indices(range), [&](auto i) {out[i] = fun(range[i]);}, g, tag);
            return out;
        }
    }

    // template <class R=DefaultReducer, class V>
    // auto reduce(V const &v, R const &r=DefaultReducer()) const {
    //     return fork(executor, [&](auto const &ex) {return ex.reduce(v, r);});
    // }

    template <class R=DefaultReducer, class V, class I, class F, class T=OrderedSplit>
    I map_reduce(V &&range, I init, F const &fun, R const &op={}, GrainSize g=1, T tag={}) const {
        auto results = map(fw<V>(range), fun, g, tag);
        return std::accumulate(move_begin(results), move_end(results), std::move(init), op);
    }

};


/******************************************************************************************/

}
