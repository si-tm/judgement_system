#pragma once
#include "../iteration/Spreadsort.h"
#include "../types/Complex.h"
#include "Algebras.h"
#include "Backtrack.h"
#include "Action.h"

namespace nupack { namespace thermo {

/******************************************************************************************/

using mark_t = vec<uint>;

/** Takes a list of marks and then associates each mark with a random value in
  the interval [0, value] and returns the list sorted by the random values.
*/
vec<std::pair<real, uint>> compute_weights(mark_t const & marks, DefaultRNG &);

/******************************************************************************************/

/** Determines which recursion matrix in block contains element t and then adds
  the corresponding Segment to each of the samples in used.
  This overload is for subelements from within the same block.
*/
template <class Q, class T, class U, NUPACK_IF(!traits::has_recursion_strand<T>)>
void push_segment(Ignore, Ignore, Q &queue, Ignore, T const &t, U used) {
    NUPACK_QUICK_REQUIRE(t.i, >=, 0);
    NUPACK_QUICK_REQUIRE(t.j, >=, 0);
    extend(queue.try_emplace(Segment(t.i, t.j, t.recursion.priority())).first->second, item_view(used));
}

// This overload is for subelements from different blocks.
template <class M, class T, class U, NUPACK_IF(traits::has_recursion_strand<T>)>
void push_segment(int i, int j, Ignore, M &queues, T const &t, U used) {
    NUPACK_QUICK_REQUIRE(t.i, >=, 0);
    NUPACK_QUICK_REQUIRE(t.j, >=, 0);
    auto &q = t.recursion.strand >= 0 ? queues(i, i + t.recursion.strand) : queues(i + -t.recursion.strand, j);
    extend(q.try_emplace(Segment(t.i, t.j, t.recursion.priority())).first->second, item_view(used));
}

/******************************************************************************************/

/** Core of the sampling algorithm. Finds the matrix element in the block
  refered to by the segment and replays the N4 computation of this element,
  adding new segments to marked samples as the partial sum crosses their
  corresponding weights.
*/
template <class Constants, class Q, class M>
void sample_element(uint i, uint j, Constants const &c, Q &queue, M &queues, Segment const &s, mark_t const & marks, DefaultRNG &rng) {
    backtrack(c.block, s, [&](auto const &m, auto value) {
        // compute weights in [0:1]
        auto w = compute_weights(marks, rng);
        auto weights = view(w);
        using Rig = typename Constants::model_type::rig_type;

        // If accum ends up slightly less than 1 and there are weights remaining, rescale so that those weights will be hit
        for (uint time = 0; time != 3; ++time) {
            auto const inverse_value = Rig::invert()(value);
            real accum = 0.0;

            c.recurse(m.expression(c.multi(), s.i, s.j, c), c.one(), [&](auto const s, auto const &...ns) {
                NUPACK_QUICK_REQUIRE(len(weights), >, 0);

                accum += Rig::as_exponential(Rig::times()(s, inverse_value));

                if (front(weights).first < accum) {
                    auto const it = std::upper_bound(weights.begin()+1, weights.end(), accum, [](auto const &a, auto const &p) {return a < p.first;});
                    auto const used = view(begin_of(weights), it);
                    (push_segment(i, j, queue, queues, ns, used), ...);
                    weights.set_begin(it);
                    return weights.empty(); // short circuit if no more weights
                } else return false;
            });

            if (weights.empty()) break;
            else if (time == 0) Rig::adjust_exponential(value, accum); // shrinks value by a factor of accum to deal with numerical imprecision.
            else if (time == 1) Rig::adjust_exponential(value, 0.999); // back up case when the numerical precision is still an issue
            else NUPACK_QUICK_ASSERT(weights.empty(), "contribution not found", i, j, s, weights, value, accum, inverse_value);
        }
    });
}

/******************************************************************************************/

struct SampleQueue {
    using Queue = std::map<Segment, mark_t>;
    Triangle<Queue> queues;
    SequenceList sequences;

    SampleQueue(SequenceList k) : queues(k.size()), sequences(std::move(k)) {}

    template <class V, class Constants>
    void consume(V &samples, Constants const &c, uint i, uint j, DefaultRNG &rng, bool print_segments=false) {
        auto queue = std::move(queues(i, j));
        auto const i_start = sum(view(sequences, 0, i), len),
                   j_start = sum(view(sequences, 0, j), len);

        while (!queue.empty()) {
            auto const begin = queue.begin();
            auto const cur = std::move(*begin);
            queue.erase(begin);
            if (print_segments) print("popping: ", cur.first, cur.second);

            /* samples are updated on pop of "B" matrix element. If coaxial stacking and
            dangle states are eventually captured, which will need to be done at
            the time Segments are pushed, this logic should be moved to the same
            point for consistency. */
            if (cur.first.priority == B.priority())
                for (auto i : cur.second) samples[i].add_pair(i_start + cur.first.i, j_start + cur.first.j);
            sample_element(i, j, c, queue, queues, cur.first, cur.second, rng);
        }
    }
};


/******************************************************************************************/

}}
