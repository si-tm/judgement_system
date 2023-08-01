/**
 * @brief Backtracking through dynamic programs
 */
#pragma once
#include "CachedModel.h"
#include "../types/Sequence.h"
#include "../types/PairList.h"
#include "../types/Triangle.h"
#include "../common/Random.h"
#include "../iteration/Search.h"

#include <map>
#include <set>
#include <utility>

namespace nupack::thermo {

/******************************************************************************************/

/** Represents an element of a recursion matrix that has been sampled. Includes
  matrix type (e.g. "Q", "MS", etc.) priority is determined automatically
  during lookup from the structure of the Block object, and thus Segment works
  without modification for both basic and coaxial matrix sets.
*/
struct Segment : MemberComparable, WeaklyOrdered {
    uint i, j;
    Priority priority;
    NUPACK_REFLECT(Segment, i, j, priority)

    constexpr Segment(uint i, uint j, Priority p) : i(i), j(j), priority(p) {}

    friend std::ostream & operator<<(std::ostream &os, Segment const &t) {
        return os << recursion_name(t.priority) << '(' << t.i << ", " << t.j << ')';
    }

    // Bigger elements have a longer range [i:j] and a higher priority
    // Include i as a tie-breaker in comparison
    bool operator<(Segment const &b) const noexcept {
        return std::forward_as_tuple(b.i + j, priority, b.i) < std::forward_as_tuple(i + b.j, b.priority, i);
    }
};

/******************************************************************************************/

struct StrandSegment : Segment, MemberOrdered {
    NUPACK_EXTEND_REFLECT(StrandSegment, Segment, left, right);
    uint left, right;

    constexpr StrandSegment(uint l, uint r, Segment const &s) : Segment(s), left(l), right(r) {}

    template <class T>
    StrandSegment(T const &t, uint l, uint r) : Segment(t.i, t.j, t.recursion.priority()), left(l), right(r) {
        NUPACK_QUICK_REQUIRE(t.i, >=, 0);
        NUPACK_QUICK_REQUIRE(t.j, >=, 0);
        if constexpr(traits::has_recursion_strand<T>) {
            if (t.recursion.strand >= 0) {
                NUPACK_QUICK_REQUIRE(l + t.recursion.strand, <=, r, l, r, t.recursion.strand);
                right = l + t.recursion.strand;
            } else {
                NUPACK_QUICK_REQUIRE(l + -t.recursion.strand, <=, r, l, r, t.recursion.strand);
                left = l + -t.recursion.strand;
            }
        }
    }

    friend std::ostream & operator<<(std::ostream &os, StrandSegment const &t) {
        return os << recursion_name(t.priority) << '[' << t.left << ", " << t.right << "](" << t.i << ", " << t.j << ')';
    }

    // Compare this way so strands with a larger span and lower start are considered bigger
    bool operator<(StrandSegment const &b) const noexcept {
        return std::forward_as_tuple(right + b.left, b.left, static_cast<Segment const &>(*this))
             < std::forward_as_tuple(b.right + left, left, static_cast<Segment const &>(b));
    }
};

/******************************************************************************************/

template <class B, class F>
auto backtrack(B const &block, Segment const &s, F &&f) {
    bool done = false;
    for_each(members_of(block), [&](auto const &m) {
        if constexpr(can_backtrack<decltype(m)>) {
            if (!done && m.priority() == s.priority) {f(m, m(s.i, s.j)); done = true;}
        }
    });
    NUPACK_ASSERT(done, "Could not find matching recursion", s);
}

/******************************************************************************************/

}
