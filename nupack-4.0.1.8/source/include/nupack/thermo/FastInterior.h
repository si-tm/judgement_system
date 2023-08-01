/// Fast interior loops

#pragma once
#include "Common.h"
#include "Tensor.h"
#include "Big.h"

namespace nupack::thermo {

/******************************************************************************************/

// Partition function of [i, j] given that i, j close an extensible interior loop
// At a given diagonal in the algorithm, the X(i, s) is the partition function of loops that stem from within (i:j) of size s. 
// j is implied by the current diagonal
template <class Constants, NUPACK_IF(Constants::complexity::value == 3 && Constants::is_forward::value)>
auto fast_interior_pairs(SingleStrand, Index i, Index j, Constants const &c) {return c.dot(c(X, i, lspan(8, j-i-5)));}

template <class Constants, NUPACK_IF(Constants::complexity::value == 4 && Constants::is_forward::value)>
auto fast_interior_pairs(SingleStrand, Index i, Index j, Constants const &c) {
    return c.total(range(5, min(c.model.int_max-3, j-i-8)), [=, &c](Index s) {
        cspan R(5, min(c.model.int_max+2-s, j-i-s-3));
        return c.dot(c(YA, R+i, j-s), c.model.int_size(R+s-2), c.model.int_asym(s, R));
    });
}

template <class Constants, NUPACK_IF(!Constants::is_forward::value)>
auto fast_interior_pairs(SingleStrand, Index i, Index j, Constants const &c) {
    return c.total(range(10, min(c.model.int_max+2, j-i-3)), [=, &c](auto z) {
        return c.total(range(5, z-4), [=, &c](auto s) {
            return c.product(c(YA, i+s, j+s-z), c.model.int_size(z-2), c.model.int_asym(s, z-s));
        });
    });
}

/******************************************************************************************/

template <class Constants, NUPACK_IF(Constants::complexity::value == 3 && Constants::is_forward::value)>
auto fast_interior_pairs(MultiStrand, Index i, Index j, Constants const &c) {
    return c.dot(c(X, i, lspan(min(8, j + c.lsize() - i - 2), j + c.lsize() - i - 2)));
}

template <class Constants, NUPACK_IF(Constants::complexity::value == 4 && Constants::is_forward::value)>
auto fast_interior_pairs(MultiStrand, Index i, Index j, Constants const &c) {
    return c.total(range(5, min(c.model.int_max-3, j+1)), [=, &c] (Index s) {
        cspan R(5, min(c.model.int_max+2-s, c.lsize() - i));
        return c.dot(c(YA, R+i, j-s), c.model.int_size(R+s-2), c.model.int_asym(s, R));
    });
}

template <class Constants, NUPACK_IF(!Constants::is_forward::value)>
auto fast_interior_pairs(MultiStrand, Index i, Index j, Constants const &c) {
    //  z := j-e + d-i
    // s:= d-i
    return c.total(lrange(10, min(c.model.int_max+2, c.lsize()-i+j)), [=, &c](auto z) {
        return c.total(lrange(max(5, max(j, z)-j), min(z-4, c.lsize()-i)), [=, &c](Index s) {
            NUPACK_QUICK_REQUIRE(i+s, <, c.lsize());
            NUPACK_QUICK_REQUIRE(j+s-z, >=, 0);
            NUPACK_QUICK_REQUIRE(int(j+s)-z, >=, 0);
            NUPACK_QUICK_REQUIRE(s, >=, 5);
            NUPACK_QUICK_REQUIRE(z-s, >=, 5);
            return c.product(c(YA, i+s, j+s-z), c.model.int_size(z-2), c.model.int_asym(s, z-s));
        });
    });
}

/******************************************************************************************/

template <class X0, class X2, class Constants>
bool calculate_fast_interior(SingleStrand, X0 &&x, X2 const &x2, Index i, Index j, Constants const &c) {
    if (i+15 < j) {
        auto const y1 = c(YA, i+5, span(i+9, j-6));
        auto const y2 = c(YA, span(i+7, j-8), j-5);
        cspan K(10, j-i-5);
        if constexpr(std::is_arithmetic_v<typename X2::value_type>) {
            simd::map(x, K, [&](auto const k) {
                return c.plus(
                    c.times(c.model.int_scale(k-2), x2(k-2)),
                    c.times(c.model.int_asym(k), c.plus(y1(j-i-6-k), y2(k-10))));
            });
        } else {
            simd::map(x, K, [&](auto const k) {
                // auto const e = -x2.exponent(k-2);
                auto const e = simd::max_of(x2.exponent(k-2), y1.exponent(j-i-6-k), y2.exponent(k-10));
                return std::make_pair(c.plus(
                    c.ldexp(c.times(c.model.int_scale(k-2), x2.mantissa(k-2)), x2.exponent(k-2) - e), // calculation of mantissa in the exponent of x2
                    c.times(c.model.int_asym(k), c.plus(
                        c.ldexp(y1.mantissa(j-i-6-k), y1.exponent(j-i-6-k) - e),
                        c.ldexp(y2.mantissa(k-10), y2.exponent(k-10) - e)))), e);
            });
        }
    }
    if (i+13 < j)  x.set(8, c.times(c.model.int_asym(8), c(YA, i+5, j-5)));
    if (i+14 == j) x.set(9, c.times(c.model.int_asym(9), c(YA, i+6, j-5)));
    if (i+14 < j)  x.set(9, c.times(c.model.int_asym(9), c.plus(c(YA, i+6, j-5), c(YA, i+5, j-6))));
    return true;
}

/******************************************************************************************/

template <class X0, class X2, class Constants>
bool calculate_fast_interior(MultiStrand, X0 &&x, X2 const &x2, Index i, Index j, Constants const &c) {
    if (c.on_left(i, +6) && c.on_right(j, -6)) {
        auto const y1 = c(YA, i+5, span(0, j));
        auto const y2 = c(YA, span(i+5, c.lsize()), j-5);
        Index const r = j+4, s = c.lsize()-i+3, p = j+c.lsize()-i-2;

        if constexpr(std::is_arithmetic_v<typename X2::value_type>) {
            simd::map(x, 10, min(r, s), [&](auto k) {
                // calculation of mantissa in exponent of x2
                return c.plus(c.times(x2(k-2), c.model.int_scale(k-2)),
                    c.times(c.model.int_asym(k), c.plus(y1(j+3-k), y2(k-8))));
            });

            simd::map(x, r, s, [&](auto k) {
                return c.plus(c.times(x2(k-2), c.model.int_scale(k-2)), c.times(c.model.int_asym(k), y2(k-8)));
            });

            simd::map(x, s, r, [&](auto k) {
                return c.plus(c.times(x2(k-2), c.model.int_scale(k-2)), c.times(c.model.int_asym(k), y1(j+3-k)));
            });

            simd::map(x, max(r, s), p, [&](auto k) {
                return c.times(x2(k-2), c.model.int_scale(k-2));
            });
        } else {
            simd::map(x, 10, min(r, s), [&](auto k) {
                auto const e = simd::max_of(x2.exponent(k-2), y1.exponent(j+3-k), y2.exponent(k-8));
                return std::make_pair(c.plus(c.ldexp(c.times(x2.mantissa(k-2), c.model.int_scale(k-2)), x2.exponent(k-2) - e),
                    c.times(c.model.int_asym(k), c.plus(
                        c.ldexp(y1.mantissa(j+3-k), y1.exponent(j+3-k) - e),
                        c.ldexp(y2.mantissa(k-8), y2.exponent(k-8) - e)))), e);
            });
            simd::map(x, r, s, [&](auto k) {
                auto const e = simd::max_of(y2.exponent(k-8), x2.exponent(k-2));
                return std::make_pair(c.plus(c.ldexp(c.times(x2.mantissa(k-2), c.model.int_scale(k-2)), x2.exponent(k-2) - e),
                              c.ldexp(c.times(c.model.int_asym(k), y2.mantissa(k-8)), y2.exponent(k-8) - e)), e);
            });
            simd::map(x, s, r, [&](auto k) {
                auto const e = simd::max_of(y1.exponent(j+3-k), x2.exponent(k-2));
                return std::make_pair(c.plus(c.ldexp(c.times(x2.mantissa(k-2), c.model.int_scale(k-2)), x2.exponent(k-2) - e),
                              c.ldexp(c.times(c.model.int_asym(k), y1.mantissa(j+3-k)), y1.exponent(j+3-k) - e)), e);
            });
            simd::map(x, max(r, s), p, [&](auto k) {
                return std::make_pair(c.times(x2.mantissa(k-2), c.model.int_scale(k-2)), x2.exponent(k-2));
            });
            // set_aligned_exponents(x, span(10, p), x2, span(8, p-2));
        }

        x.set(9, c.times(c.model.int_asym(9), c.plus(c(YA, i+6, j-5), c(YA, i+5, j-6))));
    }

    if (c.on_left(i, +5) && !c.on_left(i, +6) && c.on_right(j, -6)) {
        auto const y1 = c(YA, i+5, span(0, j));
        cspan R(9, j+4);
        if constexpr(std::is_arithmetic_v<typename X2::value_type>) {
            simd::map(x, R, [&](auto k) {return c.times(c.model.int_asym(k), y1(j+3-k));});
        } else {
            simd::map(x.mantissa, R, [&](auto k) {return c.times(c.model.int_asym(k), y1.mantissa(j+3-k));});
            reverse_copy(y1.exponent(j+4-R), x.exponent(R)); // put x in the same exponent as y1, reversed
        }
    }

    if (c.on_left(i, +6) && c.on_right(j, -5) && !c.on_right(j, -6)) {
        auto const y2 = c(YA, span(i+5, c.lsize()), j-5);
        cspan R(9, c.lsize()-i+3);
        if constexpr(std::is_arithmetic_v<typename X2::value_type>) {
            simd::map(x, R, [&](auto k) {return c.times(c.model.int_asym(k), y2(k-8));});
        } else {
            simd::map(x.mantissa, R, [&](auto k) {return c.times(c.model.int_asym(k), y2.mantissa(k-8));});
            forward_copy(y2.exponent(R-8), x.exponent(R));
        }
    }

    if (c.on_left(i, +5) && c.on_right(j, -5))
        x.set(8, c.times(c.model.int_asym(8), c(YA, i+5, j-5)));
    return true;
}

/******************************************************************************************/

template <class T>
struct FastInterior : X_t {
    static constexpr auto can_backtrack() {return False();}
    std::array<Matrix<T>, 3> matrices;
    NUPACK_REFLECT(FastInterior, matrices);

    void clear() noexcept {for (auto &m : matrices) m.clear();}

    template <class Constants>
    void allocate(Constants const &c, std::size_t i, std::size_t j) {
        if constexpr(Constants::complexity::value == 3 && Constants::is_forward::value)
            for (auto &x : matrices) {x.allocate(i+1, max(i + j + 1, 4) - 4); x.fill(T(c.model.zero()));}
    }

    template <class U, class F>
    void assign_and_clear(FastInterior<U> &&u, F &&f) {
        matrices[0].assign_and_clear(std::move(u.matrices[0]), f);
        matrices[1].assign_and_clear(std::move(u.matrices[1]), f);
        matrices[2].assign_and_clear(std::move(u.matrices[2]), f);
    }

    void next_diagonal() noexcept {
        swap(matrices[2], matrices[1]);
        swap(matrices[1], matrices[0]);
    }

    template <class ...Ts>
    decltype(auto) operator()(Ts const &...ts) const {return matrices[0](ts...);}

    auto & operator[](uint i) {return matrices[i];}
    auto const & operator[](uint i) const {return matrices[i];}

    // Fast interior loops matrices update for single and multiple strands
    template <class S, class Constants>
    bool calculate(S, Index i, Index j, Constants const &c) {
        if constexpr(Constants::complexity::value == 3 && Constants::is_forward::value) {
            if (matrices[0].shape()[0] < 4 || matrices[0].shape()[1] < 4) return true;
            cspan s(0, matrices[0].shape()[1]);
            auto x = matrices[0](i, s);
            auto const x2 = add_const(matrices)[2](i+1, s);
            try {
                // x.fill(0);
                // return true;
                return calculate_fast_interior(S(), x, x2, i, j, c);
            } catch (...) {
                BEEP(i, j, c.strands(), typeid(*this).name());
                throw;
            }
        } else return true;
    }
};

/******************************************************************************************/

}
