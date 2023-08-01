#pragma once
#include "Operations.h"
#include "../iteration/Patterns.h"
#include "../iteration/Range.h"

namespace nupack {

/******************************************************************************************/

struct SerialImpl {
    // template <class R, class V>
    // auto reduce(V const &v, R const &r) const {
    //     return std::accumulate(begin_of(v), end_of(v), value_type_of<V>(), r);
    // }

    template <class V, class F>
    void implement_map(V &&range, GrainSize, F const &fun, Ignore) const {
        for (auto &&x : range) fun(std::forward<decltype(x)>(x));
    }

    constexpr auto n_workers() const {return 1u;}
};

/******************************************************************************************/

}
