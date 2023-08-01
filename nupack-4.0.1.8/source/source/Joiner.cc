#include <nupack/state/StateComplexes.h>
#include <nupack/state/Joiner.h>
#include <nupack/iteration/Transform.h>

namespace nupack::kmc {

#pragma message("should go to fixed point arithmetic for Joiner?")

/******************************************************************************************/

void ProductJoiner::update_join_rates_matrix(StateComplexes const &complexes, Index i, Matrix const &new_join_rates) {
    auto const x = complexes.strand_map[i].x;
    complex_rates.update(x, complex_rates[x] + new_join_rates - join_rates[i]);
    complex_rates_squared.update(x, calculate_rate_product(complex_rates[x]));
    join_rates[i] = new_join_rates;
}

/******************************************************************************************/

void ProductJoiner::register_join(StateComplexes const &complexes, Index s1, Index s2) {
    auto x1 = complexes.strand_map[s1].x, x2 = complexes.strand_map[s2].x;

    complex_rates.update(x1, complex_rates[x1] + complex_rates[x2]);
    complex_rates_squared.update(x1, calculate_rate_product(complex_rates[x1]));
    complex_rates.swap_erase(x2); complex_rates_squared.swap_erase(x2);
}

/******************************************************************************************/

// Register that two complexes have joined. Zero all of their join interactions.
void QuadraticJoiner::register_join(StateComplexes const &complexes, Index s1, Index s2) {
    for (auto i : complexes.complex_indices[complexes.strand_map[s1].x])
        for (auto j : complexes.complex_indices[complexes.strand_map[s2].x])
            join_rates.zero(join_index(i, j));
}

/******************************************************************************************/

/// Search for "t" in the cumulative sum function C(f1) + scale * C(f2)
template <class Matrix, class Column>
auto fenwick_find_product(Fenwick<Matrix> const &f1, Fenwick<Column> const &f2, uint index, real t, real const &scale) {
    return fenwick_find<real>(t, 
        0, // prefix
        f1.total()(index, 0) * scale - f2.total()(index), // total
        f1.size(), // length
        [&](auto k) {return f1.tree[k](index, 0) * scale - f2.tree[k](index);}); // map
}

/******************************************************************************************/

JoinMove ProductJoiner::get_join_move_nondimensional(StateComplexes const &complexes, real const rate) const {
    auto const &complex_indices = complexes.complex_indices;
    auto const &strand_map = complexes.strand_map;
    JoinMove ret;
    Index x1, x2;

    // Choose an interaction index from the total array
    Column const total = la::eval(calculate_rate_product(complex_rates.total()) - complex_rates_squared.total());
    auto const [iter, r] = find_cumulative(total, rate);
    ret.index = iter - total.begin();

    // Choose the first complex, knowing what kind of base pair it is
    real margin;
    std::tie(x1, margin) = fenwick_find_product(complex_rates, complex_rates_squared, ret.index, r, complex_rates.total()(ret.index, 1));

    // Choose the second complex, which can't equal the first
    // Scale down the running margin to work in the x2 space
    margin /= complex_rates[x1](ret.index, 0);

    // There's a very minor chance that this could be hit in a multithreaded context, but really minor.
    // Generally the same state will not be hit this way from multiple threads so should be fine.
    auto save = complex_rates[x1];
    const_cast<Fenwick<Matrix> &>(complex_rates).zero(x1); // can't choose the same complex
    std::tie(x2, margin) = complex_rates.find(margin, [&](auto const &m){return m(ret.index, 1);});
    const_cast<Fenwick<Matrix> &>(complex_rates).update(x1, std::move(save));
    NUPACK_REQUIRE(x1, !=, x2, *this);

    // Choose the strand in the second complex
    NUPACK_REQUIRE(margin, <=, complex_rates[x2](ret.index, 1));
    auto const s2 = find_if(complex_indices[x2], [&](auto const &s) {return minus_if(margin, join_rates[s](ret.index, 1));});
    NUPACK_ASSERT(s2 != end_of(complex_indices[x2]));
    ret.o2 = strand_map[*s2].loop;

    // Scale back out
    margin *= complex_rates[x1](ret.index, 0) / join_rates[*s2](ret.index, 1);

    // Choose the strand in the first complex
    NUPACK_REQUIRE(margin, <=, complex_rates[x1](ret.index, 0));
    auto const s1 = find_if(complex_indices[x1], [&](auto const &s) {return minus_if(margin, join_rates[s](ret.index, 0));});
    NUPACK_ASSERT(s1 != end_of(complex_indices[x1]));
    ret.o1 = strand_map[*s1].loop;

    NUPACK_REQUIRE(margin, <=, join_rates[*s1](ret.index, 0));
    ret.margin = margin * join_rates[*s2](ret.index, 1);
    ret.scale = join_rates[*s2](ret.index, 1);

    NUPACK_REQUIRE(ret.o1, !=, ret.o2, *this);
    return ret;
}

/******************************************************************************************/

}
