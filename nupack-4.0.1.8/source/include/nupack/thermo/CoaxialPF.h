#pragma once
#include "BasicPF.h"

#ifndef NUPACK_ONLY_RECURSIONS
#include "Common.h"
#define IFF(cond, expression) ((cond) ? iff_nullable(c, expression) : ZeroHelper<decltype(c)>{c})
#define AND &&

namespace nupack::thermo {

/******************************************************************************************/

/// Data has some cached information on the sequence, linear storage complexity
//  Generally, i pairs to s, s1 is 1 after s, s1 pairs to j
template <class T>
class CoaxialRows {
    /// Simple multidimensional array access helpers
    Index index_with_s(Base i, Base j) const noexcept {return n * (+i) + (+j);} // 16
    Index index_with_j(Base i, Base s, Base s1) const noexcept {return n * n * (+i + 1) + n * (+s) + (+s1);} // 64
    Index index_with_i(Base s, Base s1, Base j) const noexcept {return n * n * (+s + 1 + n) + n * (+s1) + (+j);} // 64
    // length is 2 * M^3 + M^2, i.e. 144

public:
    using value_type = T;
    Matrix<T> data;
    Index n = 0;
    NUPACK_REFLECT(CoaxialRows, data, n);

    CoaxialRows() = default;

    template <class I>
    auto operator()(I const &i, Base s, Base s1, Base j) const {
        static_assert(!std::is_same_v<I, Base>);
        return data(index_with_i(s, s1, j), i);
    }

    auto operator()(Base i, span s, Base j) const {
        return data(index_with_s(i, j), s); // s1 is implicit
    }

    template <class J>
    auto operator()(Base i, Base s, Base s1, J const &j) const {
        static_assert(!std::is_same_v<J, Base>);
        return data(index_with_j(i, s, s1), j);
    }

    template <class J>
    auto operator()(Base i, J const &j) const {
        static_assert(!std::is_same_v<J, Base>);
        return data(data.shape()[0] - n + +i, j);
    }

    template <class Model>
    void calculate(Sequence const &s, Model const &model) {
        auto const all = model.energy_model.alphabet().all();
        n = all.size();
        data.allocate(2 * cube(n) + sq(n) + n, len(s)); // 144 = 16 + 64 + 64
        if (s.size()) for (auto b : all) {
            for (auto c : all) {
                for (auto r : range(len(s) - 1))
                    data(index_with_s(b, c), r) = model.coaxial(b, s[r], s[r+1], c);
                for (auto d : all) for (auto r : indices(s)) {
                    data(index_with_j(b, c, d), r) = model.coaxial(b, c, d, s[r]);
                    data(index_with_i(b, c, d), r) = model.coaxial(s[r], b, c, d);
                }
            }
            izip(s, [&](auto i, Base c) {
                data(data.shape()[0] - n + +b, i) = (b == c) ? model.one() : model.zero();
            });
        }
    }
};

/******************************************************************************************/

template <class Constants, class F>
auto dangle_sum(Constants const &c, Index i, Index j, F &&f) {
    // scans through k, l in range(2), range(2) - compile time constants may or may not help
    return c.sum(
        f(d0, d0),
        IFF(c.on_right(j, -1), f(d0, d1)),
        IFF(c.on_left(i, +1), f(d1, d0)),
        IFF(c.on_left(i, +1) AND c.on_right(j, -1), f(d1, d1))
    );
}

/******************************************************************************************/

// dangle_state<true>(i, j, k, l) -> dangle(j-l, j, i, i+k)
// dangle_state<false>(i, j, k, l) -> dangle(i, i+k, j-l, j)
// dangle5 1 and 2 paired
// dangle3 2 and 3 paired
template <bool Outer, class Constants>
constexpr auto dangle_energy(Constants const &c, Index i, Index j, D0, D0) {
    return decltype(c.model.dangle3_energy(Base(), Base(), Base()))(0);
}

template <bool Outer, class Constants>
auto dangle_energy(Constants const &c, Index i, Index j, D0, D1) {
    if (Outer) return c.model.dangle3_energy(c.right(j-1), c.right(j), c.left(i));
    else return c.model.dangle5_energy(c.left(i), c.right(j-1), c.right(j));
}

template <bool Outer, class Constants>
auto dangle_energy(Constants const &c, Index i, Index j, D1, D0) {
    if (Outer) return c.model.dangle5_energy(c.right(j), c.left(i), c.left(i+1));
    else return c.model.dangle3_energy(c.left(i), c.left(i+1), c.right(j));
}

template <bool Outer, class Constants>
auto dangle_energy(Constants const &c, Index i, Index j, D1, D1) {
    if (Outer) return c.model.energy_model.terminal_mismatch(c.right(j-1), c.right(j), c.left(i), c.left(i+1));
    else return c.model.energy_model.terminal_mismatch(c.left(i), c.left(i+1), c.right(j-1), c.right(j));
}

template <bool Outer, class Constants, class ...Ts>
auto dangle_state(Constants const &c, Ts const &...ts) {return c.model.boltz(dangle_energy<Outer>(c, ts...));}

#endif

/******************************************************************************************/

struct CoaxialCD : CD_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {
        span R(i, j);
        return c.sum(
            dangle_sum(c, i, j, [=, &c](auto k, auto l) {
                return IFF(j >= Index(l) AND i + k < j - l AND c.model.can_close(c.left(i+k), c.right(j-l)),
                    c.product(dangle_state<0>(c, i, j, k, l), c(D, i+k, j-l)));
            }),
            IFF(i < j, c.dot(c.linear(c.left(i), R, c.right(j)), c(D, i, R), c(D, R+1, j)))
        );
    }
    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        return c.sum(
            c.sandwich(i, j, [=, &c](auto const &L, auto const &R, span s, auto const &coax) {
                return c.dot(coax(c.left(i), s, c.right(j)), L(D, i, s), R(D, s+1, j)); // ICS 6 (i & s, s+1 & j)
            }),
            dangle_sum(c, i, j, [=, &c](auto k, auto l) {
                return IFF(c.model.can_close(c.left(i+k), c.right(j-l)),
                    c.product(dangle_state<0>(c, i, j, k, l), c(D, i+k, j-l))); // D 5
            })
        );
    }
};

struct CoaxialMD : MD_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {
        return dangle_sum(c, i, j, [=, &c](auto k, auto l) {
            return IFF(j >= Index(l) AND i + k < j - l AND c.model.can_close(c.left(i+k), c.right(j-l)),
                c.product(c(D, i+k, j-l), c.model.multi3s(k+l), dangle_state<0>(c, i, j, k, l), c.model.multi2));
        });
    }
    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        return dangle_sum(c, i, j, [=, &c](auto k, auto l) {
            return IFF(c.model.can_close(c.left(i+k), c.right(j-l)),
                c.product(c(D, i+k, j-l), c.model.multi3s(k+l), dangle_state<0>(c, i, j, k, l), c.model.multi2)); // MD 5
        });
    }
};

struct CoaxialMC : MC_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {
        span R(i, j);
        return IFF(i < j, c.product(c.dot(c.linear(c.left(i), R, c.right(j)), c(D, i, R), c(D, R+1, j)), c.model.multi22));
    }
    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        return c.sandwich(i, j, [=, &c](auto const &L, auto const &R, span s, auto &coax) {
            return c.product(c.dot(coax(c.left(i), s, c.right(j)), L(D, i, s), R(D, s+1, j)), c.model.multi22); // MICS 6 (i & s, s+1 & j)
        });
    }
};

/******************************************************************************************/

struct CoaxialMCS : MCS_t { // Coaxial stack spanning i to <= j
    template <class S, class Constants>
    auto expression(S, Index i, Index j, Constants const &c) const {
        span R(S::value ? 0 : i, j+1);
        return c.dot(c(MC, i, R), c.model.multi3rs(j+1-R)); // not sure about j+1-R
    }
};

/******************************************************************************************/

struct CoaxialMS : MS_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {
        span R(i, j+1);
        return c.sum(c(MCS, i, j), c.dot(c.model.multi3rs(j+1-R), c(MD, i, R))); // 29
    }
    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        span R(0, j+1);
        return c.sum(c(MCS, i, j), c.dot(c.model.multi3rs(j+1-R), c(MD, i, R))); // 29
    }
};

/******************************************************************************************/

struct CoaxialM : M_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {
        span R(i, j+1), S(i, j);
        return c.sum(
            c.dot(c.model.multi3s(R-i), c(MS, R, j)),
            IFF(i < j, c.dot(c(M, i, S), c(MS, S+1, j)))
        );
    }
    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        span S(i, c.lsize());
        return c.sum(
            c.dot(c.model.multi3s(S-i), c(MS, S, j)), // 38
            c.sandwich(i, j, [=, &c](auto const &L, auto const &R, span s, Ignore) {
                return c.dot(L(M, i, s), R(MS, s+1, j));
            }) // 40
        );
    }
};

/// Partition function of [i, j] given that i, j close a multiloop
// We want another version with constraint that i+1 is not paired or j-1 is not paired
struct CoaxialM2 : M2_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {
        span s(i+5, j-5);
        return IFF(i+10 < j, c.sum(
            c.dot(c.model.multi3s(span(0, j-i-2)), c(MCS, span(i+1, j-1), j-1)),
            c.dot(c(M, i+1, s), c(MS, s+1, j-1))
        ));
    }

    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        return IFF(c.on_left(i, +1) AND c.on_right(j, -1), c.sum(
            c.dot(c.model.multi3s(span(0, c.lsize()-i-1)), c(MCS, span(i+1, c.lsize()), j-1)),
            c.sandwich(i+1, j-1, [=, &c](auto const &L, auto const &R, span s, Ignore) {
                return c.dot(L(M, i+1, s), R(MS, s+1, j-1));
            })
        ));
    }
};

/******************************************************************************************/

/// c.B contributions from closing pair dangles
template <class Constants>
auto closing_pair_dangles(Index i, Index j, Constants const &c) {
    NUPACK_QUICK_ASSERT(c.model.can_close(c.left(i), c.right(j)));
    return dangle_sum(c, i, j, [=, &c](auto k, auto l) {
        return c.product(c.sum(
            // Both sides can dangle
            IFF(c.on_left(i, k+1) AND c.on_right(j, -1-l), 
                c.product(c(M2, i+k, j-l), c.model.multi12, c.model.multi3s(k+l)) // MCPD 5 and 7
            ),
            c(N, i+k, j-l),
            // ECPD 8 Neither side dangles but they're consecutive
            IFF(c.adjacent(i+k, j-l), c.model.one()) // c(Q, i+k+1, j-l-1)=1 (empty recursion)
        ), dangle_state<1>(c, i, j, k, l));
    });
}

// c.B contributions from closing pair stacking
template <class Constants>
auto closing_pair_stacks(Index i, Index j, Constants const &c) {
    NUPACK_QUICK_ASSERT(c.model.can_close(c.left(i), c.right(j)));
    return c.sum(
        IFF(c.on_left(i, +1) AND c.on_right(j, -1), c.sum(
            c.sandwich(i+1, j-1, [=, &c](auto const &L, auto const &R, span s, auto const &coax) {
                return c.sum(
                    c.product(c.dot(coax(c.right(j), c.left(i), c.left(i+1), s), L(D, i+1, s), R(M, s+1, j-1)), c.model.multi122),  // MCPCS 6 (i & j, i+1 & s)
                    c.product(c.dot(coax(s+1, c.right(j-1), c.right(j), c.left(i)), L(M, i+1, s), R(D, s+1, j-1)), c.model.multi122),  // MCPCS 8 (i & j, s+1 & j-1)
                    // reserved values: i can not be at end of strand, j can not be at beginning of strand
                    c.dot(coax(c.right(j), c.left(i), c.left(i+1), s), L(D, i+1, s), R(N, s, j)), // ECPCS 8 (i & j, i+1 & s)
                    c.dot(coax(s+1, c.right(j-1), c.right(j), c.left(i)), L(N, i, s+1), R(D, s+1, j-1))  // ECPCS 14 (i & j, s+1 & j-1)
                );
            }),
            c.nick_total([=, &c](auto const &L, auto const &R, auto const &lcoax, auto const &rcoax, auto last) {
                return c.sum(
                    c.product(L(D, i+1, last), R(Q, 0, j-1), lcoax(c.right(j), c.left(i), c.left(i+1), last)), // ECPCS 10B (i & j, i+1 & last)
                    c.product(L(Q, i+1, last), R(D, 0, j-1), rcoax(0, c.right(j-1), c.right(j), c.left(i))) // ECPCS 16B (i & j, 0 & j-1)
                );
            })
        )),
        IFF(!c.on_left(i, +1) AND  c.on_right(j, -1), c.sum(
            c.first_middle([=, &c](auto const &R, auto const &coax, auto m) {
                return c.product(R(D, m, j-1), coax(m, c.right(j-1), c.right(j), c.left(i))); // ECPCS 16C (i & j, m & j-1)
            }),
            c.sandwich(i, j-1, [=, &c](auto const &L, auto const &R, span s, auto &coax) {
                // we want Q(i+1, s)
                return c.dot(coax(s+1, c.right(j-1), c.right(j), c.left(i)), L(N, i, s+1), R(D, s+1, j-1)); // ECPCS 16A (i & j, s+1, j-1)
            })
        )),
        IFF( c.on_left(i, +1) AND !c.on_right(j, -1), c.sum(
            c.last_middle([=, &c](auto const &L, auto const &coax, auto m) {
                return c.product(L(D, i+1, m), coax(c.right(j), c.left(i), c.left(i+1), m)); // ECPCS 10C (i & j, i+1 & m)
            }), // ECPCS 16C
            c.sandwich(i+1, j, [=, &c](auto const &L, auto const &R, span s, auto &coax) {
                // we want Q(s+1, j-1)
                return c.dot(coax(c.right(j), c.left(i), c.left(i+1), s), L(D, i+1, s), R(N, s, j)); // ECPCS 10A (i & j, i+1 & s)
            })
        ))
    );
}

struct CoaxialS : S_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {
        return c.dot(c(CD, i, span(i, j+1)));
    }
    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        return c.dot(c(CD, i, span(0, j+1)) /* * dangle(i, Span, Span + 1)  */);
    }
};

template <class Constants>
auto multi_closing_pairs(Index i, Index j, Constants const &c) {
    span const s(i+1, j-1);
    return IFF(i+10 < j, c.sum(
        c.product(c.sum( // TODO optimize these spans a bit, they can probably be smaller
            c.dot(c.linear(c.right(j), c.left(i), c.left(i+1), s), c(D, i+1, s), c(M, s+1, j-1)), // MULTILOOP_CLOSING_PAIR_COAXIAL_STACKING (i & j, i+1 s)
            c.dot(c.linear(s+1, c.right(j-1), c.right(j), c.left(i)), c(M, i+1, s), c(D, s+1, j-1)) // MULTILOOP_CLOSING_PAIR_COAXIAL_STACKING (i & j, j-1, R)
        ), c.model.multi122),
        dangle_sum(c, i, j, [=, &c](auto k, auto l) {
            return c.product(c.sum(c(M2, i+k, j-l)), c.model.multi3s(k+l), c.model.multi12, dangle_state<1>(c, i, j, k, l));
        })
    ));
}

struct CoaxialB : B_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {
        return c.pairing(i, j, [&] {
            return c.sum(
                single_pair(SingleStrand(), i, j, c),
                inextensible_interior(SingleStrand(), i, j, c),
                IFF(c.model.can_close(c.left(i), c.right(j)),
                    c.product(multi_closing_pairs(i, j, c), c.model.terminal(c.right(j), c.left(i)))),
                extensible_interior(SingleStrand(), i, j, c)
            );
        });
    }
    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        return c.pairing(i, j, [&] {
            return c.sum(
                inextensible_interior(MultiStrand(), i, j, c),
                // single_pair is taken care of by closing_pair_dangles
                IFF(c.model.can_close(c.left(i), c.right(j)), c.product(c.sum(
                    closing_pair_stacks(i, j, c),
                    closing_pair_dangles(i, j, c)
                ), c.model.terminal(c.right(j), c.left(i)))),
                extensible_interior(MultiStrand(), i, j, c)
            );
        });
    }
};


/// c recursion for single and multiple strands
struct CoaxialQ : Q_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {
        return c.sum(
            c.model.one(),
            IFF(i+3 < j, c(S, i, j)),
            IFF(i+4 < j, c.dot(c(Q, i, span(i, j-4)), c(S, span(i+1, j-3), j)))
        );
    }
    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        return c.sum(
            c(S, i, j),
            c.sandwich(i, max(j, 4) - 4, [=, &c](auto const &L, auto const &R, span s, Ignore) {
                return c.dot(L(Q, i, s), R(S, s+1, j));
            })
        );
    }
};

/******************************************************************************************/

/// c.N (exterior loop) recursions
// N is the pfunc of [i+1:j-1] given there is a single nick strictly between i and j 
struct CoaxialN : N_t {
    template <class Constants>
    auto expression(SingleStrand, Index i, Index j, Constants const &c) const {return c.model.zero();}

    template <class Constants>
    auto expression(MultiStrand, Index i, Index j, Constants const &c) const {
        return c.sum(
            // IFF(!c.on_left(i, +1) AND !c.on_right(j, -1), 1e300),
            IFF(!c.on_left(i, +1) AND  c.on_right(j, -1), c.first_middle([=](auto const &R, Ignore, auto m) {NUPACK_QUICK_REQUIRE(m, ==, 0); return R(Q, m, j-1);})), // j with first strand after left
            IFF( c.on_left(i, +1) AND !c.on_right(j, -1), c.last_middle([=](auto const &L, Ignore, auto m) {return L(Q, i+1, m);})), // i with last strand before right
            IFF( c.on_left(i, +1) AND  c.on_right(j, -1), c.nick_total([=, &c](auto const &L, auto const &R, Ignore, Ignore, auto last) { /// (0, 1)(2, 2), (0, 0)(1, 2)
                return c.product(L(Q, i+1, last), R(Q, 0, j-1)); // 43
            }))
        );
    }
};

#ifndef NUPACK_ONLY_RECURSIONS
}
#undef IFF
#undef AND
#endif
