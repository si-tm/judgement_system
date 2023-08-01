#pragma once
#include "../types/Sequence.h"
#include "../types/Matrix.h"
#include "ParameterStorage.h"

// Loop free energy calculation using a dynamic program-ish method
// Linear time, constant space in the number of base pairs in the loop
// Works by matrix multiplication

namespace nupack {

/******************************************************************************************/

/// Matrix logarithm of 2 x 2 matrix, not tested, simplified in Mathematica, not incorporated yet
/// It could be used to prevent potential overflow on a pathologic multi/exterior loop
// actually not really, forgot the logs don't commute.
template <class A>
A log2x2(A const &a) {
    using T = value_type_of<A>;
    T const rad = std::sqrt(4 * a(0,1) * a(1,0) + sq(a(0,0) - a(1,1)));
    T const logq = std::log(a(0,0) + a(1,1) - rad);
    T const logp = std::log(a(0,0) + a(1,1) + rad);
    T const dif = (logp - logq) / rad;
    T const add = T(0.5) * ((logp + logq) - T(1.3862943611198906)); // log(4)

    return {{T(0.5) * dif * (a(0,0) - a(1,1)) + add, a(0,1) * dif},
            {a(1,0) * dif, T(0.5) * dif * (a(1,1) - a(0,0)) + add}};
}

/******************************************************************************************/

template <class T, class Plus, class Times>
struct StackMatrix {
    using array = std::array<std::array<T, 2>, 2>;
    array m;

    NUPACK_REFLECT(StackMatrix, m);

    StackMatrix() = default;
    StackMatrix(std::pair<T, T> const &a, std::pair<T, T> const &b) :
        m{{{{a.first, a.second}}, {{b.first, b.second}}}} {}

    template <class I, class J, NUPACK_IF(std::is_integral_v<I> && std::is_integral_v<J>)>
    T operator()(I const &i, J const &j) const {return m[i][j];}

    static T plus(T t, T u) {return Plus()(t, u);}

    friend StackMatrix operator*(StackMatrix const &a, StackMatrix const &b) {
        auto dot = [&](auto i, auto j) {
            return Plus()(Times()(a(i, 0), b(0, j)), Times()(a(i, 1), b(1, j)));
        };
        return {{dot(0, 0), dot(0, 1)},
                {dot(1, 0), dot(1, 1)}};
    }

    constexpr T trace() const {return Plus()(m[0][0], m[1][1]);}
};

/******************************************************************************************/

template <class T, class Plus, class Times>
struct StackColumn {
    using array = std::array<T, 2>;
    array m;

    NUPACK_REFLECT(StackColumn, m);

    StackColumn() = default;
    StackColumn(T const &a, T const &b) : m{a, b} {}

    template <class I, NUPACK_IF(std::is_integral_v<I>)>
    T operator()(I const &i) const {return m[i];}

    friend StackColumn operator*(StackMatrix<T, Plus, Times> const &a, StackColumn const &b) {
        return {
            Plus()(Times()(a(0, 0), b(0)), Times()(a(0, 1), b(1))),
            Plus()(Times()(a(1, 0), b(0)), Times()(a(1, 1), b(1))) 
        };
    }

    // b^T A = A^T b
    friend StackColumn operator*(StackColumn const &b, StackMatrix<T, Plus, Times> const &a) {
        return {
            Plus()(Times()(a(0, 0), b(0)), Times()(a(1, 0), b(1))),
            Plus()(Times()(a(0, 1), b(0)), Times()(a(1, 1), b(1))) 
        };
    }
};

/******************************************************************************************/

/// Returns matrix for a given base pair assuming its sequence on the left is l long, on the right is r long
template <class M, class E>
M factor(std::size_t l, std::size_t r, bool disable, E const &e) {
    // row (1, 2): this pair (is not / is) stacking towards the left
    // col (1, 2): next pair (is not / is) stacking towards the left
    // I noticed l==3 is same as l > 3 so I combined them.
    if (disable) { // this pair can't stack. used for the fake base pair in exterior loops.
        return              {{e.one(),            e.one()}, // this "pair" cannot stack to the left
                             {e.zero(),           e.zero()}}; // all configurations with the next pair stacking are allowed
    } else if (l == 2) {
        switch (r) {
            case 2:  return {{e.one(),            e.one()},
                             {e.stack(),          e.zero()}};
            case 3:  return {{M::plus(e.one(), e.right()),  e.one()},
                             {e.stack(),                    e.stack()}};
            default: return {{M::plus(e.one(), e.right()),  M::plus(e.one(), e.right())}, 
                             {e.stack(),                    e.stack()}};
        }
    } else {
        switch (r) {
            case 2:  return {{e.one(),            e.one()},
                             {e.left(),           e.zero()}};
            case 3:  return {{M::plus(e.one(), e.right()),  e.one()}, // this pair not stacking to left. allow stacking to right depending on whether next pair occupies the dangle
                             {M::plus(e.left(), e.both()),  e.left()}};
            default: return {{M::plus(e.one(), e.right()),  M::plus(e.one(), e.right())},  // this pair is not stacking to the left
                             {M::plus(e.left(), e.both()),  M::plus(e.left(), e.both())}}; // this pair is stacking to the left
        }
    }
}

/******************************************************************************************/

template <class M, bool B>
struct StackModel {
    M const &model;
    // With respect to 2 sequences bounding a base pair
    Base ip; //< base paired to front of first
    Base id; //< penultimate of first
    Base i;  //< last of first
    Base j;  //< first base of second
    Base jd; //< second base of second

    using type = if_t<B, real, typename M::value_type>; // always use double precision if exponentiating

    constexpr type boltz(type e) const {return B ? model.boltz(e) : e;}
    type stack()       const  {return ip == Base::null() ? zero() : boltz(model.coaxial_stack_energy(ip, id, i, j));}
    type left()        const  {return boltz(model.dangle3(id, i, j));}
    type right()       const  {return boltz(model.dangle5(i, j, jd));}
    type both()        const  {return boltz(model.terminal_mismatch(id, i, j, jd));}

    type one() const {return B ? 1 : 0;}
    type zero() const {return B ? 0 : std::numeric_limits<type>::infinity();}
};

// Returns a StackModel for a given sequence  v[s0]
template <bool B, class V, class M>
StackModel<M, B> stack_model(V const &v, std::size_t s0, M const &model) {
    auto const s1 = s0 ? s0 - 1 : v.size() - 1;
    auto const s2 = s1 ? s1 - 1 : v.size() - 1;
    return {model, back(v[s2]), back_index(v[s1], 1), back(v[s1]), front(v[s0]), front(v[s0], 1)};
}

// Returns a StackModel for a given sequence  v[s0]
// template <class V, class M>
// StackModel<M> stack_model(V const &v, std::size_t s0, M const &model) {
//     auto const s1 = s0 ? s0 - 1 : v.size() - 1;
//     return {model, Base::null(), back_index(v[s1], 1), back(v[s1]), front(v[s0]), front(v[s0], 1)};
// }

/******************************************************************************************/

template <class T, class V, class M>
T stacking_sum(V const &v, int const nick, M const &model) {
    using A = StackMatrix<T, std::plus<T>, std::multiplies<T>>;
    A q = factor<A>(len(back(v)), len(front(v)), nick == 0, stack_model<true>(v, 0, model));
    for (auto s : range(1, len(v)))
        q = q * factor<A>(len(v[s-1]), len(v[s]), nick == s, stack_model<true>(v, s, model));
    return q.trace();
}

/******************************************************************************************/

template <class T, class V, class M>
T minimum_stacking_energy(V const &v, int const nick, M const &model) {
    using A = StackMatrix<T, min_t, std::plus<T>>;
    A q = factor<A>(len(back(v)), len(front(v)), nick == 0, stack_model<false>(v, 0, model));
    for (auto s : range(1, len(v)))
        q = q * factor<A>(len(v[s-1]), len(v[s]), nick == s, stack_model<false>(v, s, model));
    return q.trace();
}

}
