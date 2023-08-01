#pragma once
#include "Backtrack.h"
#include "Action.h"
#include "../types/Complex.h"
#include "../standard/Deque.h"
#include "../standard/Set.h"
#include "../iteration/Spreadsort.h"

// mfe information manipulated externally when pushing segments onto stack
namespace nupack::thermo {

/******************************************************************************************/

struct SuboptQueue {
    /// Ordered map from segment to which structures contains it (needs to be ordered)
    using StructureIndex = std::uint32_t;
    using SegmentMap = std::map<Segment, vec<StructureIndex>>;
    // Set of strand segments, compared by address
    using SegmentIter = typename SegmentMap::iterator;
    struct Hash {
        std::size_t operator()(SegmentIter const &it) const {return hash_of(&*it);}
    };
    struct Equal {
        bool operator()(SegmentIter const &a, SegmentIter const &b) const {return &*a == &*b;}
    };
    using SegmentSet = std::unordered_set<typename SegmentMap::iterator, Hash, Equal>;
    
    struct Structure {
        PairList structure;
        SegmentSet segments; // remaining segments that are unrecursed for this structure
        real gap; // minimum energy of this structure - MFE
    };

    std::deque<Structure> structures; //< list of structures and their current segments (deque is good for large number of structures and not invalidating references)
    Triangle<SegmentMap> segment_maps; //< segment map for each strand block (segment -> list of containing structures)
    SequenceList sequences; //< list of all the sequences (a constant)

    SuboptQueue(SequenceList k) : segment_maps(k.size()), sequences(std::move(k)) {}

    // put a new structure from its root structure, its strand segments, and its energy difference
    template <class V>
    void add_structure(StructureIndex t, V const &s_iters, real diff) {
        // Put in new structure with copies of segments from root structure and an updated gap
        auto const new_index = structures.size();
        auto &new_structure = structures.emplace_back(structures[t]);
        new_structure.gap += diff;

        // Add new segments to the created strucutre
        new_structure.segments.insert(s_iters.begin(), s_iters.end());

        // Put new structure in each of segment's sets
        for (auto const &iter : new_structure.segments)
            iter->second.emplace_back(new_index);
    }

    auto length() const {return segment_maps.length();}

    void initialize(Segment s, PairList p) {
        auto iter = segment_maps(0, length() - 1).try_emplace(s).first;
        iter->second.emplace_back(0);
        structures.emplace_back(Structure{std::move(p), SegmentSet{iter}, 0});
    }

    SegmentIter emplace(StrandSegment const &s) {
        NUPACK_QUICK_REQUIRE(s.i, <, len(sequences[s.left]));
        NUPACK_QUICK_REQUIRE(s.j, <, len(sequences[s.right]));
        return segment_maps(s.left, s.right).try_emplace(s).first;
    }

    template <class Constants>
    void pop(SegmentMap &q, StrandSegment const s, Constants const &c, real const gap) {
        auto const iter = std::prev(q.end());
        NUPACK_QUICK_REQUIRE(iter->first, ==, static_cast<Segment const &>(s));
        vec<StructureIndex> v = std::move(iter->second); // list of structures containing the segment
        NUPACK_ASSERT(!std::empty(v));
        spreadsort_float_map(v, [&](auto i) {return structures[i].gap;}); // sort structures by mfe differentials

        // The first-encountered MFE structure is handled separately
        small_vec<SegmentIter, 2> best;

        // go through each structure and pop this segment
        for (auto const i : v) structures[i].segments.erase(iter);

        q.erase(iter); // remove segment

        backtrack(c.block, s, [&](auto const &m, auto const value) {
            NUPACK_ASSERT(std::isfinite(value), value, s);

            // Numerical precision issues probably only happens in fast interior loops.
            // But easier to just redo the recursion so we don't have this issue
            auto min_observed_value = inf<decltype(value)>();
            c.recurse(m.expression(c.multi(), s.i, s.j, c), c.one(), [&](auto const alt, auto const &...ns) {
                min_eq(min_observed_value, alt);
                return False();
            });
            NUPACK_QUICK_REQUIRE(abs(min_observed_value - value), <, 1e-2, min_observed_value, value, s);

            bool best_done = false;
            c.recurse(m.expression(c.multi(), s.i, s.j, c), c.one(), [&](auto const alt, auto const &...ns) {
                // keep track of the actually observed minimum value in rare cases of precision issues
                // difference between the alternate substructure and the minimum one
                // it should be non-negative but might drift a tiny bit negative due to floating point precision
                auto const diff = alt - min_observed_value;
                bool const same = diff <= 0;

                // if the most stable structure can't hold this within gap, then short-circuit
                if (!same && gap < structures[front(v)].gap + diff) return False();
                NUPACK_ASSERT(std::isfinite(diff), diff, alt, value, min_observed_value, s);

                if (same && !best_done) {
                    // every structure is compatible because no energy cost is being introduced. Hold these segments for later.
                    best.reserve(sizeof...(ns));
                    (best.emplace_back(emplace(StrandSegment(ns, s.left, s.right))), ...);
                    best_done = true;
                } else {
                    // get iterators for each segment of ns...
                    std::array<SegmentIter, sizeof...(ns)> const iters{emplace(StrandSegment(ns, s.left, s.right))...};

                    auto n = same ? v.size() : lower_bound(view(v).offset(1), gap - diff, [&](auto i) {return structures[i].gap;}) - v.begin();

                    // add new structures as copies of v, modified energy, new segments
                    for (auto i : view(v, 0, n)) add_structure(i, iters, diff);
                }
                return False();
            });
            NUPACK_ASSERT(best_done, "missing optimal structure", s, value, min_observed_value, gap);
        });

        if (!best.empty()) { // all structures compatible with the list of segments in first
            for (auto const &s : best) extend(s->second, v); // put all structures into each segment
            for (auto t : v) structures[t].segments.insert(best.begin(), best.end()); // put each segment into all structures
        }
    }

    auto prefix(uint n) const {return sum(view(sequences, 0, n), len);}

    template <class Constants>
    [[nodiscard]] bool consume(real gap, std::size_t max, Constants const &c, uint l, uint r) {
        gap = std::max(gap, 0.0);
        auto const i_start = prefix(l), j_start = prefix(r);
        SegmentMap &q = segment_maps(l, r);

        while (!q.empty() && structures.size() <= max) {
            StrandSegment cur(l, r, back(q).first);

            /* samples are updated on pop of "B" matrix element. If coaxial stacking and
            dangle states are eventually captured, which will need to be done at
            the time Segments are pushed, this logic should be moved to the same
            point for consistency.
            */
            if (cur.priority == B.priority())
                for (auto t : back(q).second)
                    structures[t].structure.add_pair(i_start + cur.i, j_start + cur.j);

            pop(q, std::move(cur), c, gap);
        }

        return q.empty();
    }

};

}
