/** \file System.h
 * @brief Contains System, representing all strands present in a State or set of States
 */
#pragma once
#include "../types/Complex.h"
#include "../types/Structure.h"
#include "../types/Alphabet.h"
#include "../iteration/Transform.h"
#include "../standard/Ptr.h"
#include "../standard/Set.h"
#include "StateComplexes.h"

namespace nupack {

/**************************************************************************************/

struct SystemData {
    /// Contiguous container for concatenated strands
    Sequence total_sequence;
    /// Views into total_sequence representing each strand
    SubsequenceList strands;
    /// Start/end of each strand, strand of each base
    vec<iseq> nicks, strand_map;
    Alphabet alphabet = DNA;

    NUPACK_REFLECT(SystemData, total_sequence, strands, nicks, strand_map, alphabet);

    // void make_strands() {
    //     strands.clear();
    //     if (total_sequence.empty()) return;
    //     for (auto i : range(len(nicks) - 1))
    //         strands.emplace_back(iterator_at(nicks[i]), iterator_at(nicks[i+1]));
    // }
};

/**************************************************************************************/

/// The system is stored as a contiguous, concatenated array of bases
/// for all the strands. The strands are kept as views into this sequence
/// This is easier for lookup, but worse for dynamic modification
/// of the system contents. However, the system is so small in memory
/// that it is probably easiest to reserve a large chunk of sequence space
/// so that iterators will not be invalidated
struct System : public ConstIndexable<System>, TotallyOrdered {
    std::shared_ptr<SystemData const> data;
    using StrandIter = SubsequenceList::const_iterator;

    NUPACK_REFLECT(System, data);

    System() = default;

    /// Make a System from a container of strands
    System(SequenceList const &v);
    
    /**************************************************************************************/

    Sequence const &sequence() const {return data->total_sequence;}
    SubsequenceList const & iter() const {return data->strands;}
    auto const & strands() const {return data->strands;}

    SequenceList to_complex() const {return indirect_view(strands(), [](auto const &v) {return view(v).offset(1, -1);});}

    bool operator<(System const &o) const {return data != o.data && sequence() < o.sequence();}
    bool operator==(System const &o) const {return data == o.data || sequence() == o.sequence();}
    auto hash() const {return hash_of(sequence());}

    /// Return index corresponding to iterator in the total sequence
    auto index(BaseIter it) const {return it - sequence().data();}
    /// Return starting index of a strand from the strand iterator
    auto begin_of_strand(StrandIter it) const {return at(data->nicks, it - begin_of(data->strands));}
    /// Return past-the-end index of a strand from the strand iterator
    auto end_of_strand(StrandIter it) const {return at(data->nicks, it - begin_of(data->strands) + 1);}
    /// Return next strand within a loop structure recursion
    template <class V>
    StrandIter next_strand_it(value_type_of<V> j, V const &pairs) const {
        while (j != begin_of_strand(strand_it_of(j))) j = pairs[--j]; 
        return strand_it_of(j);
    }
    /// Return strand index of a sequence iterator
    auto strand_of(std::size_t i) const {return data->strand_map[i];}
    /// Return strand index of a sequence index
    auto strand_of(BaseIter it) const {return strand_of(it - sequence().data());}
    /// Return strand iterator of a sequence index
    StrandIter strand_it_of(int loc) const {return next(data->strands, data->strand_map[loc]);}
    /// Return whether a position i is the past-the-end index of
    bool is_strand_end(iseq i) const {
        if (i == len(sequence())) return true;
        return (sequence()[i] == Base::null() && sequence()[i - 1] == Base::null());
    }
    /// Return iterator in sequence from an index
    BaseIter iterator_at(iseq i) const {return sequence().data() + i;}
    Base base_at(iseq i) const {return sequence()[i];}
    BaseIter total_begin() const {return sequence().data();}
    BaseIter total_end() const {return sequence().data() + sequence().size();}
    /// Number of nucleotides in the system
    auto n_bases() const {return len(sequence()) - 2 * len(data->strands);}

    auto save_repr() const {return vmap<SequenceList>(data->strands, view);}

    void load_repr(SequenceList const &seqs) {if (!seqs.empty()) *this = System(seqs);}
};

/******************************************************************************************/

template <class V, class ...Ts>
void build_complex(V &loops, System const &s, PairList const &pairs, vec<System::StrandIter> &strands, Ts const &...ts) {
    std::vector<std::tuple<usize, iseq, iseq>> queue;
    if constexpr(has_capacity<V>) queue.reserve(loops.capacity() - loops.size());
    queue.emplace_back(len(loops), s.begin_of_strand(strands.front()), s.end_of_strand(strands.front()) - 1);
    loops.emplace_back(len(loops), Ether, strands.front()->begin(), ts...);

    while (!queue.empty()) {
        auto &q = queue.back(); // LIFO so that a depth first search is used
        auto const index = std::get<0>(q);

        if (loops[index].next_pair(s, len(loops), std::get<1>(q), std::get<2>(q), pairs, strands)) {
            auto d = std::get<1>(q), e = pairs[d];
            std::get<1>(q) = e;
            NUPACK_QUICK_REQUIRE(d, !=, e);
            queue.emplace_back(len(loops), d, e);
            loops.emplace_back(len(loops), index, s.iterator_at(d), ts...);
        } else {
            queue.pop_back();
        }
    }
}

/******************************************************************************************/

/// Recurse through a complex of loops calling a callback on each new loop, takes starting index and returns ending index
// Edge build_complex(vec<StrandIter> &its, Edge, PartialLoop, PairList const &pairs, Callback const &o) const;
/// Recurse through each complex in a PairList and return a StateComplexes
template <class V, class ...Ts>
StateComplexes build_complex_set(V &loops, System const &s, PairList const &pairs, Ts const &...ts) {
    if (len(pairs) != sum(s.strands(), len))
        NUPACK_ERROR("number of nucleotides doesn't match length of pair list", len(pairs), sum(s.strands(), len));

    StateComplexes out{len(s.strands())};

    auto pool = make_set(iterators(s.strands()));
    while(!pool.empty()) {
        auto p = *begin_of(pool); // take the first strand not incorporated already
        vec<System::StrandIter> strands = {p};
        build_complex(loops, s, pairs, strands, ts...);
        // reconstruct the order of strands in the complex
        out.emplace_back(vmap<StateComplexes::Indices>(strands, [&](auto it) {
            pool.erase(it);
            return it - begin_of(s.strands());
        }));
    }
    return out;
}

/******************************************************************************************/

}
