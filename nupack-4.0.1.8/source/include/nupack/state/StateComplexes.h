/** \file StateComplexes.h
 * @brief Contains StateComplexes, a class representing all strands/complexes present in a State
 */
#pragma once

#include "../model/Move.h"

#include "../iteration/Search.h"
#include "../algorithms/Operators.h"
#include "../algorithms/Utility.h"

#include "../types/Fenwick.h"
#include "../types/IO.h"
#include "../standard/Variant.h"

namespace nupack {

/******************************************************************************************/

struct StateComplexes : ConstIndexable<StateComplexes>, CompareByKey<StateComplexes> {
    using Index = iseq;
    using Indices = vec<iseq>;

    struct StrandData {
        Index x, pos, loop;
        NUPACK_REFLECT(StrandData, x, pos, loop);
    };

    /**************************************************************************************/

    /// Strand # -> (Complex #, index within complex, loop #)
    vec<StrandData> strand_map;
    /// Complex # -> Ordered strand indices in that complex
    vec<Indices> complex_indices;
    
    NUPACK_REFLECT(StateComplexes, strand_map, complex_indices);

    /**************************************************************************************/

    /// Check that indices of strand_map and complex_indices are consistent
    bool check() const;

    /**************************************************************************************/

    StateComplexes() = default;
    
    StateComplexes(std::size_t n) : strand_map(n) {}

    /// Reserve space for "s" number of strands
    //void reserve(Index s) {strand_map.resize(s); join_rates.resize(s, );}
    /// Add complex by giving its strand indices
    template <class... Args> void emplace_back(Args&&...);
    /// Reorder indices so strand "s" is first in its complex (unused)
    void rotate(Index s);
    /// Update loop index of strand "i" to "o"
    void set_loop_index(Index i, Index o) {strand_map[i].loop = o;}
    /// Register that strands i and j have joined
    void register_join(Index i, Index j);
    /// Register that strands i and j have split. Return list of new strands in formed complex and the old complex index.
    std::pair<Indices, Index> register_split(Index i, Index j);

    auto const & iter() const {return complex_indices;}

    auto const & compare_key() const {return complex_indices;}

    /**************************************************************************************/

    friend std::ostream & operator<<(std::ostream &os, StateComplexes const &c) {dump_os(os, c.complex_indices); return os;}
};

/******************************************************************************************/

template <class... Args>
void StateComplexes::emplace_back(Args&&... args) {
    complex_indices.emplace_back(fw<Args>(args)...);
    for (auto i : indices<iseq>(complex_indices.back()))
        strand_map[complex_indices.back()[i]] = {static_cast<iseq>(complex_indices.size() - 1), i, static_cast<iseq>(-1)};
}

/******************************************************************************************/

}
