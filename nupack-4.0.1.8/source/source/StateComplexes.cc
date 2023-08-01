#include <nupack/state/StateComplexes.h>
#include <nupack/state/StateComplexes.h>
#include <nupack/iteration/Transform.h>

namespace nupack {

/******************************************************************************************/

void StateComplexes::rotate(Index s) {
    auto const x = strand_map[s].x, pos = strand_map[s].pos; // complex and location within
    std::rotate(complex_indices[x].begin(), complex_indices[x].begin() + pos, complex_indices[x].end());
    for (auto i : indices<Index>(complex_indices[x])) strand_map[complex_indices[x][i]].pos = i;
};

/******************************************************************************************/

bool StateComplexes::check() const {
    // BEEP(*this);
    for (auto p : indices<Index>(strand_map))
        if (p != complex_indices[strand_map[p].x][strand_map[p].pos]) return false;
    return true;
}

/******************************************************************************************/

std::pair<StateComplexes::Indices, StateComplexes::Index> StateComplexes::register_split(Index p, Index k) {
    // Loops must be from same complex
    NUPACK_ASSERT(check());
    NUPACK_REQUIRE(strand_map[p].x, ==, strand_map[k].x, "tried to split two strands from different complexes", strand_map, *this);
    NUPACK_REQUIRE(complex_indices[strand_map[p].x][strand_map[p].pos], ==, p);
    NUPACK_REQUIRE(complex_indices[strand_map[k].x][strand_map[k].pos], ==, k);

    Index const old_ix = strand_map[p].x, new_ix = complex_indices.size();

    Index i1, i2; std::tie(i1, i2) = std::minmax(strand_map[k].pos, strand_map[p].pos);

    complex_indices.emplace_back();
    auto &oldx = complex_indices[old_ix], &newx = complex_indices[new_ix];
    NUPACK_REQUIRE(i2, <, len(oldx));

    newx.assign(oldx.begin() + i1, oldx.begin() + i2);
    oldx.erase(oldx.begin() + i1, oldx.begin() + i2);
    
    for (auto j : indices<Index>(oldx)) {
        strand_map[oldx[j]].x = old_ix;
        strand_map[oldx[j]].pos = j;
    }
    for (auto j : indices<Index>(newx)) {
        strand_map[newx[j]].x = new_ix;
        strand_map[newx[j]].pos = j;
    }
    return {newx, old_ix};
}

/******************************************************************************************/

void StateComplexes::register_join(Index s1, Index s2) {
    auto x1 = strand_map[s1].x, x2 = strand_map[s2].x;
    auto pos1 = strand_map[s1].pos, pos2 = strand_map[s2].pos;
    auto & xi1 = complex_indices[x1], & xi2 = complex_indices[x2];

    xi1 = catted<vec<Index>>(xi1.begin(), xi1.begin() + pos1, xi2.begin() + pos2, xi2.end(),
                             xi2.begin(), xi2.begin() + pos2, xi1.begin() + pos1, xi1.end());

    for (auto i : indices<Index>(xi1)) {
        strand_map[complex_indices[x1][i]].x = x1;
        strand_map[complex_indices[x1][i]].pos = i;
    }

    if (x2 + 1 != complex_indices.size()) for (auto i : complex_indices.back()) strand_map[i].x = x2;
    swap_erase(complex_indices, x2);
}

/******************************************************************************************/

}
