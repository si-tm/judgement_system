#pragma once
#include <cmath>
#include <algorithm>

namespace nupack { namespace concentration {

template <class T>
constexpr T sq(T t) {return t * t;}

template <class T>
T log_sum_exp(T const &t, T const &u) {
    auto const p = std::minmax(t, u);
    return p.second + std::log(1 + std::exp(p.first - p.second));
}

template <class V, class F>
void izip(V &&v, F &&f) {
    uword i = 0;
    for (auto &&x : v) f(i++, x);
}

/******************************************************************************************/

/// solve quadratic formula with first coefficient=1, return both roots
template <class T>
auto quadratic_solve(T const b, T const c) {
    auto r1 = (-b - std::sqrt(sq(b) - 4 * c)) / 2; // should be optimized away
    auto r2 = (-b + std::sqrt(sq(b) - 4 * c)) / 2;
    return std::make_pair(r1, r2);
}

/// solve quadratic formula, return both roots
template <class T>
auto quadratic_solve(T a, T b, T c) {return quadratic_solve(b / a, c / a);}

/******************************************************************************************/

/// Return the rotational symmetry number of a sequence (1 if there is no symmetry)
template <class V>
std::size_t rotational_symmetry(V const &v) {
    auto const b = std::begin(v), e = std::end(v);
    std::size_t const n = std::distance(b, e);
    for (std::size_t i = 1; i <= n / 2; ++i) if (n % i == 0)
        if (std::equal(b+i, e, b, e-i) && std::equal(b, b+i, e-i, e)) return n / i;
    return 1;
}

/******************************************************************************************/

template <class V, class T, class F>
auto maxmap(V &&v, T const &init, F &&f) {
    std::decay_t<decltype(f(*std::begin(v)))> out = init;
    for (auto &&x : v) {auto tmp = f(x); if (tmp > out) out = tmp;}
    return out;
}

/******************************************************************************************/

struct LexicographicalCompare {
    template <class T, class U>
    bool operator()(T const &t, U const &u) const {
        return std::lexicographical_compare(std::begin(t), std::end(t), std::begin(u), std::end(u));
    }
};

/******************************************************************************************/

}}
