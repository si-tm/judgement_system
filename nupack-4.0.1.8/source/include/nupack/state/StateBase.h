/** \file StateBase.h
 * @brief Contains non-templated base classes for StaticState and JumpState objects
 */
#pragma once
#include "StateComplexes.h"
#include "System.h"
#include "../types/PairList.h"
#include "../standard/Variant.h"

namespace nupack {

/******************************************************************************************/

struct StateBase {
    using is_member_ordered = True;

    StateBase() = default;
    explicit StateBase(System, PairList p={});

    NUPACK_REFLECT(StateBase, sys, pairs, complexes);

    /// Pointer to system for access to strands
    System sys;
    /// Vector of which bases are paired
    PairList pairs;
    /// Complex data for the state
    StateComplexes complexes;
    /// Unpseudoknotted dot-parens
    string dp() const;
    /// Base sequences ordered to match dot-parens
    string sequence() const;
    /// Number of bases
    auto n_bases() const {return sys.n_bases();};
    /// Number of strands
    auto n_strands() const {return len(sys.strands());};
    /// Number of strands
    auto n_complexes() const {return len(complexes);};
    /// Align another State's pairs to agree with the strand ordering of this one
    PairList aligned_pairs(StateBase const &) const;
    /// All possible alignments of another State's pairs to agree with the strand ordering of this one
    vec<PairList> all_aligned_pairs(StateBase const &, bool only_one) const;

    std::size_t symmetry() const;

    /// Print the state
    friend std::ostream & operator<<(std::ostream &os, StateBase const &w) {
        dump_os(os, "State('", w.sequence(), "', '", w.dp(), "')");
        return os;
    }

    /**************************************************************************************/

    // auto hash() const {return sys ? hash_of(*sys, pairs) : hash_of(pairs);}

    auto operator^(StateBase const &w) const {return pairs ^ w.pairs;}

    static constexpr auto repr_names() {return make_names("pairs", "sys");}
};

/******************************************************************************************/

NUPACK_DEFINE_TYPE(is_state_base, StateBase);

}
