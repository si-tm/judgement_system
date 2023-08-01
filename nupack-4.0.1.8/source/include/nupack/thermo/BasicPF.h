/**
 * @brief Dynamic programming recursions for non-coaxial stacking codes
 *
 * @file BasicPF.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#ifndef NUPACK_ONLY_RECURSIONS
#include "Algebras.h"
#include "../types/IO.h"
#include "../iteration/View.h"
#include "../iteration/Patterns.h"
#include "../types/Sequence.h"

#define IFF(cond, expression) ((cond) ? iff_nullable(c, expression) : ZeroHelper<decltype(c)>{c})
#define AND &&

namespace nupack::thermo {

template <class C>
struct ZeroHelper {
    C const &c;

    template <class T, NUPACK_IF(!std::is_reference_v<T> && can_construct<T, decltype(c.zero())>)>
    constexpr operator T() const {return c.zero();}

    template <class T, NUPACK_IF(!std::is_reference_v<T> && !can_construct<T, decltype(c.zero())> && math::Splat<T, decltype(c.model.zero())>::valid::value)>
    constexpr operator T() const {return math::Splat<T, decltype(c.model.zero())>::from(c.model.zero());}
};

template <class T, class C, NUPACK_IF(std::is_scalar_v<std::decay_t<T>>)>
constexpr T iff_nullable(C const &c, T &&t) {return t;}

template <class T, class C, NUPACK_IF(!std::is_scalar_v<std::decay_t<T>>)>
constexpr auto iff_nullable(C const &c, T &&t) {return c.nullable(fw<T>(t));}

/******************************************************************************************/

/// Enum for denoting which part of a block has/had to be calculated
enum class Region : char {upper='U', lower='L', all='A', cached='C'};

/******************************************************************************************/

template <class Recursion, class T, class Layout>
struct RecursionMatrix : Recursion, Matrix<T, Layout> {
    static constexpr auto can_backtrack() {return True();}

    using base_type = Matrix<T, Layout>;
    void next_diagonal() const noexcept {}

    template <class S, class Constants, class ...Ts>
    bool calculate(S, Index i, Index j, Constants const &c, Ts const &...ts) noexcept {
        auto t = Recursion().expression(S(), i, j, c, ts...);
        if (!c.valid(t)) return false;
        for_each(this->layout.set(i, j), [&](auto i) {this->storage.begin()[i] = t;});
        return true;
    }

    void allocate(Ignore, std::size_t i, std::size_t j) {base_type::allocate(i, j);}
};

/******************************************************************************************/

template <class Recursion, class T, class Layout>
struct RecursionMatrix<Recursion, Big<T>, Layout> : Recursion, Matrix<Big<T>, Layout> {
    static constexpr auto can_backtrack() {return True();}

    using base_type = Matrix<Big<T>, Layout>;

    void next_diagonal() const noexcept {}

    template <class S, class Constants, class ...Ts>
    bool calculate(S, Index i, Index j, Constants const &c, Ts const &...ts) {
        Exponent<T> e = (!decltype(c.multi())::value && i == j) ? 0
            : max(c.left_plus(*this, Recursion(), i, j).exponent,
                  c.right_minus(*this, Recursion(), i, j).exponent);
        for (auto shift = 0; shift != 2 << 18; shift += std::numeric_limits<T>::max_exponent) {
            auto t = c.evaluate_mantissa(Recursion().expression(S(), i, j, c, ts...), -(e + shift));
            // if (shift) print("iterating...", e, shift, t, i, j, TypeName<S>());
            if (c.valid(t)) {
                if (t) {
                    auto const p = Constants::rig_type::ifrexp()(t);
                    NUPACK_QUICK_REQUIRE(p.first, >=, 0, p);
                    NUPACK_QUICK_REQUIRE(p.first, <, inf<T>(), p);
                    for_each(this->layout.set(i, j), [&](auto i) {
                        this->storage.mantissa.begin()[i] = p.first;
                        this->storage.exponent.begin()[i] = e + shift + p.second;
                    });
                } else {
                    // Actually this can happen (0 * inf = nan) AND the proper result is 0.
                    // NUPACK_QUICK_ASSERT(!shift, "should not have overflowed if 0", i, j, t, e, shift, c.strands(), std::numeric_limits<T>::max_exponent, typeid(*this).name());
                    for_each(this->layout.set(i, j), [&](auto i) {
                        this->storage.mantissa.begin()[i] = 0;
                        this->storage.exponent.begin()[i] = e;
                    });
                }
                return true;
            }
        }
        print("NUPACK: overflow failure", i, j, e, c.strands(), type_name(*this));
        return false;
    }

    void allocate(Ignore, std::size_t i, std::size_t j) {base_type::allocate(i, j);}
};

/******************************************************************************************/

// On the first or last strand of the multiple strand section
template <class V>
bool on_bread(int i, int j, V const &nicks) noexcept {return (i < front(nicks) && j >= back(nicks));}

// For i, j, go over all d, d+1 bases which are on the same strand, d ∈ [i, j-1]
template <class V, class F, class Algebra>
auto sandwich(int i, int j, V const &nicks, Algebra A, F &&f) noexcept {
    return A.sum(
        NUPACK_WHERE(on_bread(i+1, j, nicks)) & f(span(i, front(nicks)-1)),
        NUPACK_WHERE(on_bread(i, j-1, nicks)) & f(span(back(nicks), j)),
        NUPACK_WHERE(on_bread(i, j, nicks)) & A.total(iterators(nicks, 0, -1), [f, A](auto n) {
            return NUPACK_WHERE(n[1] - n[0] > 1) & f(span(n[0], n[1]-1));
        })
    );
}

/******************************************************************************************/

template <class T>
struct DangleData {
    using value_type = T;
    Matrix<T> data; // 35 by N. 35 = 3 + 4^2 + 4^2. This will become 3 + 2 * M^2
    uint n = 0;
    NUPACK_REFLECT(DangleData, data, n);

    DangleData() = default;

    constexpr uint index_i(Base b, Base c) const noexcept {return (n * +b) + (+c) + 3;} 
    constexpr uint index_j(Base b, Base c) const noexcept {return (n * (n + +b)) + (+c) + 3;} 

    template <class Constants>
    constexpr auto operator()(Constants const &c, Index i, Index j) const {
        static_assert(!decltype(c.multi())::value);
        NUPACK_QUICK_REQUIRE(i, <=, j);
        return c.model.boltz(c.model.energy_model.dangle_switch([&](auto dangle) -> T {
            // Mixup earlier. Model says if one side doesn't exist AND other side is positive, positive number still taken.
            bool const l = c.on_left(i, -1), r = c.on_right(j, +1);
            auto const d5 = l ? dangle.energy5(c.model.complement(c.left(i-1)), c.left(i-1), c.left(i)) : 0;
            auto const d3 = r ? dangle.energy3(c.right(j), c.right(j+1), c.model.complement(c.right(j+1))) : 0;

            if (l) {
                if (r) return dangle.combine(d3, d5, i == j); // special case if j == i
                else return d5;
            } else {
                if (r) return d3;
                else return 0;
            }
        }));
    }

    // Base i dangling on both sides
    template <class Constants>
    constexpr auto operator()(Constants const &c, Index i) const {return data(2, i);}

    // called by right_linear. i-1 always valid
    template <class Constants>
    constexpr auto operator()(Constants const &c, span i, Index j) const {
        return c.on_right(j, +1) ? data(index_j(c.right(j), c.right(j+1)), i) : data(1, i);
    }

    // called by left_linear
    template <class Constants>
    constexpr auto operator()(Constants const &c, Index i, span j) const {
        return c.on_left(i, -1) ? data(index_i(c.left(i-1), c.left(i)), j) : data(0, j);
    }

    template <class Model>
    void calculate(Sequence const &s, Model const &model) {
        n = model.energy_model.alphabet().length();
        data.allocate(3 + 2 * n * n, s.size());

        auto const all = model.energy_model.alphabet().all();
        if (s.size()) model.energy_model.dangle_switch([&](auto dangle) {
            for (auto i : indices(s)) {
                auto const e5 = i > 0 ? dangle.energy5(model.complement(s[i-1]), s[i-1], s[i]) : inf<real>();
                auto const e3 = i+1 < s.size() ? dangle.energy3(s[i], s[i+1], model.complement(s[i+1])) : inf<real>();

                data(0, i) = model.boltz(e3);
                data(1, i) = model.boltz(e5);
                data(2, i) = model.boltz(dangle.combine(e5, e3, true));

                for (auto b : all) for (auto c : all) {
                    auto const bc5 = dangle.energy5(model.complement(b), b, c);
                    auto const bc3 = dangle.energy3(b, c, model.complement(c));
                    data(index_i(b, c), i) = model.boltz(e3 + bc5);
                    data(index_j(b, c), i) = model.boltz(e5 + bc3);
                }
            }
            data(2, 0) = data(0, 0);
            data(2, s.size()-1) = data(1, s.size()-1);
        });
    }
};

template <class T>
struct NoDangleData {
    NoDangleData() = default;

    template <class U>
    NoDangleData(NoDangleData<U> const &u) {}

    auto members() {return make_members();}
    auto names() {return make_names();}

    template <class Constants>
    constexpr auto operator()(Constants const &c, Ignore, Ignore) const noexcept {return c.one();}

    void calculate(Ignore, Ignore) noexcept {}
};

#endif

/******************************************************************************************/

// Simple expressions. These ones don't depend on whether the sequence is multistranded
struct BasicYA : YA_t {
    template <class Constants>
    auto expression(Ignore, Index i, Index j, Constants const &c) const {
        return IFF(c.on_left(i, -1) AND c.on_right(j, +1), c.product(c(B, i, j),
            c.model.mismatch(c.left(i-1), c.left(i), c.right(j), c.right(j+1))));
    }
};

struct BasicYB : YB_t {
    template <class Constants>
    auto expression(Ignore, Index i, Index j, Constants const &c) const {
        return c.product(c(B, i, j), c.model.mismatch(c.left(i), c.right(j)));
    }
};

struct BasicZ : Z_t {
    template <class Constants>
    auto expression(Ignore, Index i, Index j, Constants const &c) const {
        return c.product(c(B, i, j), c.model.terminal(c.left(i), c.right(j)));
    }
};

struct BasicD : D_t {
    template <class Constants>
    auto expression(Ignore, Index i, Index j, Constants const &c) const {
        return IFF(c.model.can_close(c.left(i), c.right(j)), c(Z, i, j));
    }
};

/******************************************************************************************/

// TODO bonuses missing from interior loop contributions

// Partition function of [i, j] given that i, j close a small interior loop (both sides have < 4 unpaired nucleotides)
// Small x Small interior loops
template <class Constants>
auto inextensible_interior(SingleStrand, Index i, Index j, Constants const &c) {
    return c.total(lrange(i+1, i+5), [=, &c](Index d) {
        return c.total(lrange(max(j-4, d+4), j), [=, &c](Index e) {
            return c.product(c(B, d, e), c.model.interior(c.left(i, d+1), c.right(e, j+1)));
        });
    });
}

template <class Constants>
auto inextensible_interior(MultiStrand, Index i, Index j, Constants const &c) {
    return c.total(range(i+1, min(c.lsize(), i+5)), [=, &c](Index d) { // d between i AND first nick
        return c.total(range(max(4, j) - 4, j), [=, &c](Index e) { // e between last nick AND j
            return c.product(c(B, d, e), c.model.interior(c.left(i, d+1), c.right(e, j+1)));
        });
    });
}

/******************************************************************************************/

// Partition function of [i, j] given that i, j close a non-special case interior loop
template <class Constants>
auto extensible_interior(SingleStrand, Index i, Index j, Constants const &c) {
    auto asymmetric = [&] (auto z, auto Y) { // Large asymmetric loops
        span r(i+4+z, j-5);
        return IFF(i+z+9 < j, c.sum(c.dot(c.model.int_rsize(r-i-z-z, z), c(Y, i+1+z, r+1)),
                                    c.dot(c.model.int_size(r-i-z-z, z), c(Y, r+1-z, j-1-z))));
    };
    auto const mm = IFF(j > i, c.model.mismatch(c.right(j-1), c.right(j), c.left(i), c.left(i+1)));// : c.model.zero();
    span r(4, j-i-5);
    return c.sum(
        IFF(i+9 < j, c.product(c.sum(c.dot(c(Z, i+1, r+i+1), c.model.rbulge(r)),
                                     c.dot(c(Z, r+i+1, j-1), c.model.bulge(r))), c.model.terminal(c.right(j), c.left(i)))),
        c.product(asymmetric(1, YB), c.model.mismatch(c.right(j), c.left(i))),
        c.product(c.sum(asymmetric(2, YA), asymmetric(3, YA),
                        IFF(13 < j-i, fast_interior_pairs(SingleStrand(), i, j, c))), mm)
    );
}

template <class Constants>
auto extensible_interior(MultiStrand, Index i, Index j, Constants const &c) {
    span r(4, j), s(4, c.lsize()-i-1);
    return c.sum(
        IFF(c.on_left(i, +1) AND c.on_right(j, -5), c.product(c.dot(c(Z, i+1, r-4), c.model.rbulge(r)), c.model.terminal(c.right(j), c.left(i)))),
        IFF(c.on_left(i, +5) AND c.on_right(j, -1), c.product(c.dot(c(Z, s+i+1, j-1), c.model.bulge(s)), c.model.terminal(c.right(j), c.left(i)))),
        c.product(c.sum( // Large asymmetric loops 1
            IFF(c.on_left(i, +5) AND c.on_right(j, -2), c.dot(c.model.int_size(s-1, 1), c(YB, s+i+1, j-2))),
            IFF(c.on_left(i, +2) AND c.on_right(j, -5), c.dot(c.model.int_rsize(r-1, 1), c(YB, i+2, r-4)))),
            c.model.mismatch(c.right(j), c.left(i))
        ),
        c.product(c.sum( // Large asymmetric loops 2
            IFF(c.on_left(i, +3) AND c.on_right(j, -5), c.dot(c.model.int_rsize(r-2, 2), c(YA, i+3, r-4))),
            IFF(c.on_left(i, +4) AND c.on_right(j, -5), c.dot(c.model.int_rsize(r-3, 3), c(YA, i+4, r-4))),
            IFF(c.on_left(i, +5) AND c.on_right(j, -3), c.dot(c.model.int_size(s-2, 2), c(YA, s+i+1, j-3))),
            IFF(c.on_left(i, +5) AND c.on_right(j, -4), c.dot(c.model.int_size(s-3, 3), c(YA, s+i+1, j-4))),
            IFF(c.on_left(i, +5) AND c.on_right(j, -5), fast_interior_pairs(MultiStrand(), i, j, c))),
            IFF(c.on_left(i, +1) AND c.on_right(j, -1), c.model.mismatch(c.right(j-1), c.right(j), c.left(i), c.left(i+1)))
        )
    );
}

/******************************************************************************************/

#define NUPACK_BONUS c.one()

/// MS recursion for single AND multiple strands. 
// Sum of all configurations where there is a single base pair from i to (i, j]
struct BasicMS : MS_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {
        // span s(i, j-3);
        // return IFF(i + 3 < j, c.product(c.dot(c(D, i, s+4), c.model.multi3rs(s-i)), c.model.multi2)); // only for QB non GU
        span s(i, j-4);
        return c.product(c.sum(
            c(D, i, j), 
            IFF(i + 4 < j, c.product(c.dot(c(D, i, s+4), c.model.multi3rs(s+1-i)), NUPACK_BONUS)) // only for QB non GU //BONUS
        ), c.model.multi2);
    }
    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        span r(0, j); // no nicks between d AND j
        return c.product(c.sum(
            c(D, i, j),
            IFF(j, c.product(c.dot(c(D, i, r), c.model.multi3rs(r+1)), NUPACK_BONUS)) //BONUS
        ), c.model.multi2);
    }
};

/// MS recursion for single AND multiple strands
struct DangleMS : MS_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {
        span r(i, j-5);
        return c.product(c.sum(
            IFF(i+3 < j, c(D, i, j)), // only for QB non GU
            IFF(i+4 < j, c.product(c(D, i, j-1), c.model.multi3, c.linear(c, j, j))), // only for QB non GU //BONUS
            IFF(i+5 < j, c.dot(c(D, i, r+4), c.model.multi3rs(r-i+2), c.linear(c, r+5, j))) // only for QB non GU //BONUS
        ), c.model.multi2);
    }
    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        span r(0, j-1); // no nicks between d AND j
        return c.product(c.sum(
            c(D, i, j),
            IFF(0 < j, c.product(c(D, i, j-1), c.model.multi3, c.right_linear(c, j))), //BONUS
            IFF(1 < j, c.dot(c(D, i, r), c.model.multi3rs(r+2), c.right_linear(c, r+1, j))) //BONUS
        ), c.model.multi2);
    }
};

/******************************************************************************************/

/// M1 recursion for single AND multiple strands - no dangles
// Sum of all multiloop configurations where there is exactly least one base pair in [i, j]
struct BasicM1 : M1_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {
        span r(i+1, j-3);
        return c.sum(
            c(MS, i, j),
            IFF(i+4 < j, c.product(c.dot(c(MS, r, j), c.model.multi3s(r-i)), NUPACK_BONUS)) // SPLIT AND BONUS
        );
    }
    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        span s(i+1, c.lsize()); // no nicks between i AND e
        return c.sum(
            c(MS, i, j),
            IFF(i+1 < c.lsize(), c.product(c.dot(c(MS, s, j), c.model.multi3s(s-i)), NUPACK_BONUS))
        );
    }
};


/// M recursion for single AND multiple strands - no dangles
// Sum of all multiloop configurations where there is at least one base pair in [i, j]
struct BasicM : M_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {
        span r(i+1, j-3), s(i, j-4);
        return c.sum(
            c(MS, i, j),
            IFF(i+4 < j, c.product(c.sum(
                c.dot(c(MS, r, j), c.model.multi3s(r-i)), // SPLIT AND BONUS
                c.dot(c(M, i, s), c(MS, s+1, j))
            ), NUPACK_BONUS))
        );
    }
    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        // no nicks between i AND e
        span s(i+1, c.lsize());
        return c.sum(
            c(MS, i, j),
            c.sum(c.product(c.sum(
                IFF(i+1 < c.lsize(), c.dot(c(MS, s, j), c.model.multi3s(s-i))),
                c.sandwich(i, j, [=, &c](auto const &L, auto const &R, span s, Ignore) {
                    return c.dot(L(M, i, s), R(MS, s+1, j));
                })
            ), NUPACK_BONUS))
        );
    }
};

// struct BasicA : A_t {
//     template <class Constants>
//     auto expression(SingleStrand, Index i, Index j, Constants const &c) const {
//         return c.model.multi3s(j+1-i);
//     }
//     template <class Constants>
//     auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
//         return c.zero();
//     }
// };

/******************************************************************************************/

/// M recursion for single AND multiple strands
struct DangleM : M_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {
        span const r(i, j-4), s(i+1, j-4);
        return c.sum(
            IFF(i+3 < j, c(MS, i, j)),
            IFF(i+4 < j, c.dot(c(M, i, r), c(MS, r+1, j))),
            IFF(i+4 < j, c.product(c(MS, i+1, j), c.linear(c, i, i), c.model.multi3)),
            IFF(i+5 < j, c.dot(c(MS, s+1, j), c.linear(c, i, s), c.model.multi3s(s+1-i)))
        );
    }
    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        span s(i+2, c.lsize());
        return c.sum(
            c.product(c(MS, i, j)),
            IFF(c.on_left(i, +1), c.product(c(MS, i+1, j), c.model.multi3, c.left_linear(c, i))),
            IFF(c.on_left(i, +2), c.dot(c(MS, s, j), c.model.multi3s(s-i), c.left_linear(c, i, s-1))),
            c.sandwich(i, j, [=, &c](auto const &L, auto const &R, span s, Ignore) {
                return c.dot(L(M, i, s), R(MS, s+1, j));
            })
        );
    }
};

/******************************************************************************************/

// Sum of all multiloop configurations where there are 2 or more base pairs in [i, j] (inclusive)
struct BasicM2 : M2_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {
        span s(i+4, j-4);
        return IFF(i+8 < j AND i != 0 AND j+1 < c.rsize(), c.product(c.dot(c(M, i, s), c(MS, s+1, j)), NUPACK_BONUS));
    }

    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        return IFF(i != 0 AND j+1 < c.rsize(), c.product( 
            c.sandwich(i, j, [=, &c](auto const &L, auto const &R, span s, Ignore) {
                return c.dot(L(M, i, s), R(MS, s+1, j));
            }), NUPACK_BONUS
        ));
    }
};

// Sum of all multiloop configurations where there are 3 or more base pairs in [i, j] (inclusive)
struct BasicM3 : M3_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {
        span s(i+8, j-4);
        return IFF(i+13 < j AND i != 0 AND j+1 < c.rsize(), c.product(
            c.dot(c(M2, i, s), c(MS, s+1, j)), NUPACK_BONUS));
    }

    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        return IFF(i != 0 AND j+1 < c.rsize(), c.product(
            c.sandwich(i, j, [=, &c](auto const &L, auto const &R, span s, Ignore) {
                return c.dot(L(M2, i, s), R(MS, s+1, j));
            }), NUPACK_BONUS
        ));
    }
};

/******************************************************************************************/

/// S recursion for single AND multiple strands - no dangles
struct BasicS : S_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {
        if constexpr(Constants::is_forward::value) {
            // sum of D(i, e) where e >= i+4 AND e <= j
            return IFF(i+3 < j, c.sum(c(S, i, j-1), c(D, i, j))); // only for non GU  // SPLIT AND BONUS
        } else {
            return IFF(i+3 < j, c.dot(c(D, i, span(i+4, j+1)))); // only for non GU
        }
    }
    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        if constexpr(Constants::is_forward::value) {
            // sum of D(i, e) where e >= 0 AND e <= j
            return c.sum(IFF(c.on_right(j, -1), c(S, i, j-1)), c(D, i, j));
        } else {
            return c.dot(c(D, i, span(0, j+1))); // higher complexity but more direct for backtracking
        }
    }
};

/// S recursion for single AND multiple strands
struct DangleS : S_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {
        span s(i+4, j-1);
        return c.sum(
            IFF(i+3 < j, c(D, i, j)),
            IFF(i+4 < j, c.product(c(D, i, j-1), c.linear(c, j))),
            IFF(i+5 < j, c.dot(c(D, i, s), c.linear(c, s+1, j)))
        );
    }
    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        span s(0, j-1);
        return c.sum(
            c(D, i, j),
            IFF(0 < j, c.product(c(D, i, j-1), c.right_linear(c, j))),
            IFF(1 < j, c.dot(c(D, i, s), c.right_linear(c, s+1, j)))
        );
    }
};

/******************************************************************************************/

/// c.N (exterior loop) recursions
// N is the pfunc of [i:j] given there is a single nick strictly between i AND j
struct BasicN : N_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {return c.model.zero();}

    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        return c.product(c.nick_total([=, &c](auto const &L, auto const &R, Ignore, Ignore, auto last) { /// (0, 1)(2, 2), (0, 0)(1, 2)
            return c.product(L(Q, i, last), R(Q, 0, j)); // 43
        }), NUPACK_BONUS);
    }
};

/******************************************************************************************/

/// Total partition function of [i, j]
struct BasicQ : Q_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {
        return c.sum(
            c.linear(c, i, j),
            IFF(i+3 < j, c(S, i, j)),
            IFF(i+4 < j, c.product(c.dot(c(Q, i, span(i, j-4)), c(S, span(i+1, j-3), j)), NUPACK_BONUS))
        );
    }
    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        return c.sum(
            c(S, i, j),
            c.product(c.sandwich(i, max(j-4, 0), [=, &c](auto const &L, auto const &R, span s, Ignore) {
                return c.dot(L(Q, i, s), R(S, s+1, j));
            }), NUPACK_BONUS)
        );
    }
};

/******************************************************************************************/

/// Partition function of [i, j] given that i, j are paired AND have no nested loops
template <class Constants>
auto single_pair(SingleStrand, Index i, Index j, Constants const &c) {
    return c.model.hairpin(c.left(i, j+1));
}

template <class Constants>
auto single_pair(MultiStrand, Index i, Index j, Constants const &c) {
    return IFF(c.model.can_close(c.left(i), c.right(j)), c.product(c.sum(
        IFF(c.adjacent(i, j), c.model.one()), // last base of left next to first base of right
        IFF(!c.on_left(i, +1) &&  c.on_right(j, -1), c.first_middle([=](auto const &R, Ignore, auto m) {NUPACK_QUICK_REQUIRE(m, ==, 0); return R(Q, m, j-1);})), // j with first strand after left
        IFF( c.on_left(i, +1) && !c.on_right(j, -1), c.last_middle([=](auto const &L, Ignore, auto m) {return L(Q, i+1, m);})), // i with last strand before right
        IFF( c.on_left(i, +1) &&  c.on_right(j, -1), c.nick_total([=, &c](auto const &L, auto const &R, Ignore, Ignore, auto last) { /// (0, 1)(2, 2), (0, 0)(1, 2)
            return c.product(L(Q, i+1, last), R(Q, 0, j-1)); // 43
        }))
    ), c.model.terminal(c.right(j), c.left(i))));
}

/******************************************************************************************/

/// B matrix - B(i, j) is partition function of [i, j] given that i, j are paired
struct BasicB : B_t {
    template <class W, class Constants>
    auto expression(W, Index i, Index j, Constants const &c) const {
        return c.product(c.pairing(i, j, [&] {return c.sum(
            single_pair(W(), i, j, c),
            inextensible_interior(W(), i, j, c),
            IFF(c.model.can_close(c.left(i), c.right(j)) && c.on_left(i, +1) && c.on_right(j, -1),
                c.product(c(M2, i+1, j-1), c.model.terminal(c.right(j), c.left(i)), c.model.multi1, c.model.multi2)),
            extensible_interior(W(), i, j, c)
        );}), NUPACK_BONUS);
    }
};

#ifndef NUPACK_ONLY_RECURSIONS
}
#undef IFF
#undef AND
#endif
