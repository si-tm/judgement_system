/**
 * @brief Dynamic program algebras for forward and backward recursions
 *
 * @file Algebras.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "Rigs.h"
#include "Big.h"
#include "../standard/Optional.h"
#include "../types/Matrix.h"

namespace nupack::thermo {

/******************************************************************************************/


template <class Rig>
struct AlgebraBase {
    using rig_type = Rig;

    template <class T>
    static constexpr bool valid(T const &t) noexcept {return Rig::valid(t);}

    // struct zero_t : decltype(Rig::zero()) {};
    constexpr auto zero() const noexcept {return Rig::zero();}

    // struct one_t : decltype(Rig::one()) {};
    constexpr auto one() const noexcept {return Rig::one();}

    template <class M, class ...Is>
    static constexpr decltype(auto) element(M &&m, Is &&...is) {return static_cast<M &&>(m)(static_cast<Is &&>(is)...);}

    template <class T, class U>
    decltype(auto) ldexp(T const &t, U const &u) const noexcept {return Rig::ldexp()(t, u);}

    template <class ...Ts>
    decltype(auto) plus(Ts const &...ts) const noexcept {return fold(Rig::plus(), ts...);}

    template <class ...Ts>
    decltype(auto) times(Ts const &...ts) const noexcept {return fold(Rig::times(), ts...);}
};

/******************************************************************************************/

template <class Rig>
struct ScalarAlgebra : AlgebraBase<Rig> {
    using is_forward = True;

    /// (Ts +...) where all the Ts are scalars
    template <class ...Ts>
    auto sum(Ts const &...ts) const noexcept {return fold(Rig::plus(), ts...);}

    /**************************************************************************************/

    /// (Ts *...) where all the Ts are scalars
    template <class T, class ...Ts>
    auto product(T const &t, Ts const &...ts) const noexcept {
        return fold(Rig::times(), t, ts...);
    }

    /**************************************************************************************/

    // sum(range, f) sums the scalar outputs of f on the elements of the range
    template <class T, class F>
    auto total(T const &t, F &&f) const noexcept {
        std::decay_t<decltype(f(*begin_of(t)))> out = Rig::zero();
        for (auto &&i : t) Rig::plus_eq()(out, f(i));
        return out;
    }

    // sum(Ts[:] *...)
    template <class T, class ...Ts>
    auto dot(T const &t, Ts const &...ts) const noexcept {
        // This seems faster locally at least on an M1 mac mini
        if constexpr(std::is_same_v<Rig, PF<double>> && sizeof...(Ts) == 1) {
            if (auto const n = len(t); n > 64)
                return arma::dot(arma::Col<double>(const_cast<double *>(&t(0)), n, false), 
                                 arma::Col<double>(const_cast<double *>(&ts(0))..., n, false));
        }
        return simd::map_reduce(Rig::plus(), indices(t), [&](auto i) {return fold(Rig::times(), t(i), ts(i)...);}, Rig::sum());
    }

    /**************************************************************************************/

    template <class E>
    constexpr E nullable(E &&e) const noexcept {return static_cast<E &&>(e);}
};

/******************************************************************************************/

template <class Rig>
struct BigAlgebra : AlgebraBase<Rig> {
    using is_forward = True;

    template <class T, class E, NUPACK_IF(std::is_scalar_v<T>)>
    auto evaluate_mantissa(T const &t, E const &e) const noexcept {return Rig::ldexp()(t, e);}

    template <class T, class E, NUPACK_IF(std::is_invocable_v<T const &, E const &>)>
    auto evaluate_mantissa(T const &t, E const &e) const noexcept {return t(e);}

    template <class T, class E>
    auto evaluate_mantissa(Big<T> const &t, E const &e) const noexcept {return std::ldexp(t.mantissa, t.exponent + e);}////Rig::ldexp()(t, e);

    template <class T, class E>
    auto evaluate_mantissa(std::optional<T> const &t, E const &e) const noexcept {return t ? evaluate_mantissa(*t, e) : Rig::zero();}////Rig::ldexp()(t, e);

    /// (Ts +...) where all the Ts are scalars --> returns an expression with no natural exponent
    template <class ...Ts>
    auto sum(Ts const &...ts) const noexcept {
        return [=] (auto const &e) {
            return fold(Rig::plus(), evaluate_mantissa(ts, e)...);
        };
    }

    // /**************************************************************************************/

    /// (Ts *...) where all the Ts are scalars
    template <class T, class ...Ts>
    auto product(T const &t, Ts const &...ts) const noexcept {
        // one option is to multiply all the mantissas, then ldexp the result to the requested exponent
        // instead, we scale the first result down and then multiply
        return [=] (auto const &e) {
            return fold(Rig::times(), mantissa(ts)...,
                evaluate_mantissa(mantissa(t), fold(Rig::plus(), e, exponent(t), exponent(ts)...)));
        };
    }

    // /**************************************************************************************/

    // sum(range, f) sums the scalar outputs of f on the elements of the range --> returns expression with no natural exponent
    template <class T, class F> auto total(T const &t, F &&f) const noexcept {
        return [=] (auto const &e) {
            decay<decltype(evaluate_mantissa(f(*begin_of(t)), e))> out = Rig::zero();
            for (auto &&i : t) Rig::plus_eq()(out, evaluate_mantissa(f(static_cast<decltype(i) &&>(i)), e));
            return out;
        };
    }

    // sum(Ts[:] *...) --> returns expression with no natural exponent
    template <class T, class ...Ts> auto dot(T const &t, Ts const &...ts) const {
        return [=] (auto const &e) {
            auto map = [&](auto const &i) {
                return Rig::ldexp()(fold(Rig::times(), mantissa_at(t, i), mantissa_at(ts, i)...),
                                    fold(Rig::plus(), e, exponent_at(t, i), exponent_at(ts, i)...));
            };
            return simd::map_reduce(Rig::plus(), indices(t), std::move(map), Rig::sum());
        };
    }

    /**************************************************************************************/

    struct zero_t : decltype(Rig::zero()) {};
    // constexpr auto zero() const {return zero_t();}

    template <class T>
    struct Conditional : Optional<T> {
        Conditional(T const &t) : Optional<T>(t) {}
        Conditional(T &&t) : Optional<T>(std::move(t)) {}
        Conditional(zero_t) {}

        template <class E=simd::Zero>
        friend auto mantissa(Conditional const &c, E e={}) noexcept {return c ? mantissa(*c, e) : Rig::zero();}
        template <class E=simd::Zero>
        friend auto exponent(Conditional const &c, E e={}) noexcept {return c ? exponent(*c, e) : simd::Zero();}
    };

    struct Maybe {
        template <class T>
        friend constexpr Conditional<decay<T>> operator&(Maybe, T &&t) noexcept {return fw<T>(t);}
    };

    constexpr auto maybe() const noexcept {return Maybe();}
    auto zero() const noexcept {return std::nullopt;}

    template <class E>
    constexpr auto nullable(E &&e) const noexcept {
        static_assert(!std::is_scalar_v<std::decay_t<E>>);
        return std::make_optional(static_cast<E &&>(e));
    }
};

/**************************************************************************************/

struct Recursible_Base {};

template <class F>
struct Recursible : Recursible_Base {
    F expression;
    Recursible(F f) : expression(std::move(f)) {}
};

template <class F>
auto recursible(F &&f) noexcept {
    static_assert(is_class<decay<F>>, "should pass in functor");
    return Recursible<decay<F>>{fw<F>(f)};
}

NUPACK_DETECT(is_recursible, void_if<(derives_from<T, Recursible_Base>)>);

/**************************************************************************************/

NUPACK_DETECT(has_recursion_strand, decltype(std::declval<T>().recursion.strand));

template <class R>
struct StrandRecursion : R {
    int strand;
    constexpr StrandRecursion(int s) : strand(s) {}
};

template <class R>
constexpr StrandRecursion<R> strand_recursion(R, int s) noexcept {static_assert(std::is_empty_v<R>); return s;}

/**************************************************************************************/

template <class Rig, bool ShortCircuit>
struct BackwardAlgebra : AlgebraBase<Rig> {
    using is_forward = False;
    using Bool = if_t<ShortCircuit, bool, False>;
    using type = typename Rig::value_type;

    /// This returns if we should stop or not. If not short circuiting, return compile time False
    static constexpr Bool short_circuit(bool b) noexcept {return if_c<ShortCircuit>(b, False());}

    // In the following, fun is assumed to be (result, expressions...) -> bool
    // "recurse" is called with fun and a variadic number of expressions
    // The first expression will not be directly evaluated if it is recursible
    // All other expressions will be directly evaluated

    template <class T, class I, class J>
    static constexpr auto element(T &&t, I const &i, J const &j) noexcept {return t.addressed(t.recursion(), i, j);}

    template <class T, class I, class J>
    static constexpr auto element(T &&t, int s, I const &i, J const &j) noexcept {return t.addressed(strand_recursion(t.recursion(), s), i, j);}

    // template <class R, class T, class L, class I, class J>
    // static constexpr auto element(nupack::thermo::Table<R, T, L> &&t, I const &i, J const &j, uint material) noexcept {return t.table[material].addressed(t.recursion(), i, j);}

    // template <class R, class T, class L, class I, class J>
    // static constexpr auto element(nupack::thermo::Table<R, T, L> &&t, int s, I const &i, J const &j, uint material) noexcept {return t.table[material].addressed(strand_recursion(t.recursion(), s), i, j);}


    template <class R, class Q, class ...Ts>
    Bool recurse(Address<R> const &t, type s, Q const &q, Ts const &...ts) const noexcept {return q(s, t, ts...);}

    template <class F, class Q, class T, class ...Ts, NUPACK_IF(std::is_base_of_v<Addressed, F>)>
    Bool recurse(F const &fun, type s, Q const &q, T const &t, Ts const &...ts) const noexcept {
        return recurse(t, this->times(s, fun.value), q, ts..., fun.address());
    }

    template <class F, class Q, NUPACK_IF(std::is_base_of_v<Addressed, F>)>
    Bool recurse(F const &fun, type s, Q const &q) const noexcept {
        return q(this->times(fun.value, s), fun.address());
    }

    template <class F, class Q, class ...Ts, NUPACK_IF(std::is_invocable_v<F, type, Q, Ts...>)>
    Bool recurse(F const &fun, type s, Q const &q, Ts const &...ts) const noexcept {return fun(s, q, ts...);}

    template <class F, class Q, class ...Ts>
    Bool recurse(std::optional<F> const &fun, type s, Q const &q, Ts const &...ts) const noexcept {
        return fun ? recurse(*fun, s, q, ts...) : False();
    }

    template <class F, class Q, NUPACK_IF(is_numeric<F>::value)>
    Bool recurse(F const &f, type s, Q const &q) const noexcept {return q(this->times(s, f));}

    template <class F, class Q, class T, class ...Ts, NUPACK_IF(is_numeric<F>::value)>
    Bool recurse(F const &f, type s, Q const &q, T const &t, Ts const &...ts) const noexcept {
        return recurse(t, this->times(s, f), q, ts...);
    }

    template <class T, class ...Ts>
    auto product(T const &t, Ts const &...ts) const noexcept {
        return [=](type s, auto const &q, auto const &...factors) noexcept -> Bool {
            return recurse(t, s, q, ts..., factors...);
        };
    }

    // From variadic possibilities, recurse at most a single contribution
    template <class ...Ts>
    auto sum(Ts const &...ts) const noexcept {
        return [=](type s, auto const &q, auto const &...factors) noexcept -> Bool {
            return short_circuit((recurse(ts, s, q, factors...) || ...));
        };
    }

    // From an Iterable recurse at most a single contribution
    template <class T, class F>
    auto total(T const &t, F &&map) const noexcept {
        return [=](type s, auto const &q, auto const &...factors) noexcept -> Bool {
            for (auto const &i : t) if (recurse(map(i), s, q, factors...)) return short_circuit(true);
            return short_circuit(false);
        };
    }

    // From a range of products recurse a single contribution
    template <class T, class ...Ts>
    auto dot(T const &t, Ts const &...ts) const noexcept {
        return [=](type s, auto const &q, auto const &...factors) noexcept -> Bool {
            for (auto const i : indices(t))
                if (product(t(i), ts(i)...)(s, q, factors...)) return short_circuit(true);
            return short_circuit(false);
        };
    }

    auto zero() const noexcept {return std::nullopt;}

    template <class E>
    constexpr auto nullable(E &&e) const noexcept {
        static_assert(!std::is_scalar_v<std::decay_t<E>>, "just use normal boolean logic");
        return std::make_optional(static_cast<E &&>(e));
    }
};

/******************************************************************************************/

template <class T>
struct ForwardLSE {
    using is_forward = True;
    using rig_type = LSE<T>;

    T one() const noexcept {return rig_type::one();}
    T zero() const noexcept {return rig_type::zero();}

    bool valid(T x) const noexcept {return rig_type::valid(x);}

    template <class ...Ts>
    auto plus(Ts &&...ts) const noexcept {return rig_type::plus()(fw<Ts>(ts)...);}

    template <class ...Ts>
    auto sum(Ts &&...ts) const noexcept {return rig_type::plus()(fw<Ts>(ts)...);}

    template <class ...Ts>
    auto product(Ts &&...ts) const noexcept {return fold(rig_type::times(), fw<Ts>(ts)...);}
    
    template <class ...Ts>
    auto times(Ts &&...ts) const noexcept {return fold(rig_type::times(), fw<Ts>(ts)...);}

    template <class T1, class ...Ts>
    T dot(T1 &&t, Ts const &...ts) const {
        T const m = simd::map_reduce(simd::max, indices(t), [&](auto i) {return fold(simd::plus, t(i), ts(i)...);}, simd::reduce_max);
        if (m <= zero()) return zero();
        T const e = simd::map_reduce(simd::plus, indices(t), [&, mi=-m](auto i) {return simd::exp2(simd::tag<LSEVersion>(), fold(simd::plus, mi, t(i), ts(i)...));}, simd::reduce_sum);
        return simd::log2(simd::tag<LSEVersion>(), e) + m;
    }

    template <class V, class F>
    T total(V const &v, F &&f) const {
        T out = zero();
        for (auto const &x : v) out = plus(out, f(x));
        return out;
    }

    template <class M, class ...Is>
    static constexpr decltype(auto) element(M &&m, Is &&...is) {return static_cast<M &&>(m)(static_cast<Is &&>(is)...);}

    template <class E>
    constexpr E nullable(E &&e) const noexcept {return static_cast<E &&>(e);}
};

/******************************************************************************************/

}
