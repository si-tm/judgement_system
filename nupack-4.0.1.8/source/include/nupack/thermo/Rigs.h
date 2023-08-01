/**
 * @brief Algebraic rigs used in NUPACK dynamic programs
 *
 * @file Rigs.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "../common/Error.h"
#include "../reflect/Print.h"
#include "../math/SIMD.h"

namespace nupack::thermo {

/******************************************************************************************/

struct first_arg {
    template <class T, class U>
    constexpr auto operator()(T &&t, U) const {
        static_assert(!is_same<decay<T>, simd::Zero>, "");
        static_assert(is_same<decay<U>, simd::Zero>, "");
        return fw<T>(t);}
};

/******************************************************************************************/

/// In the PF ring, + is plus() and * is times()
template <class T>
struct PF {
    // static constexpr uint value = 0;
    using value_type = T;

    static constexpr T log_zero() {return -std::numeric_limits<T>::infinity();}
    static constexpr T zero() {return 0;}
    static constexpr T one() {return 1;}


    static constexpr auto plus() {return simd::plus;}
    static constexpr auto times() {return simd::multiplies;}

    static constexpr auto plus_eq() {return ::nupack::plus_eq;}
    static constexpr auto sum() {return simd::reduce_sum;}
    static constexpr auto ldexp() {return simd::ldexp;}
    static constexpr auto ifrexp() {return simd::ifrexp;}

    static constexpr auto invert() {return simd::reciprocal;}

    static constexpr bool valid(T const t) {
        // NUPACK_QUICK_REQUIRE(t, >=, 0); SIMD operations seem to have this fail in event of overflow.
        return likely(!std::signbit(t) && std::isfinite(t));
    }

    template <class U>
    static constexpr real as_logarithm(U const &u) {
        real o = 0;
        simd::set_logarithm(o, u);
        return o;
    }
    static constexpr T energy_scale(T const beta) {return -beta / std::log(2.0);}
    static constexpr T boltz(T const energy, T const energy_scale) {return std::exp2(energy * energy_scale);}
    static constexpr T free_energy(T const factor, T const energy_scale) {return std::log2(factor) / energy_scale;}
    
    template <class U>
    static real as_exponential(U const &factor) {return static_cast<real>(factor);}

    template <class U>
    static void adjust_exponential(U &value, T scale) {value = times()(value, scale);}
};

static constexpr int LSEVersion = 4;

template <class T>
struct LSE {
    using value_type = T;

    struct impl {
        template <class ...Ts>
        auto operator()(Ts const &...ts) {
            auto const m = simd::eval_if_simd(fold(simd::max, ts...));
            return simd::eval_if_simd(m + simd::log2(simd::tag<LSEVersion>(), fold(simd::plus, simd::exp2(simd::tag<LSEVersion>(), ts - m)...)));
        }
    };

    // Because of assumptions in log2 implementation, zero must be kept finite
    static constexpr T zero() {return 0.00390625 * std::numeric_limits<T>::lowest();}      // log(0)
    static constexpr T one() {return 0;}       // log(1)
    static constexpr bool valid(T const t) {return std::isfinite(t);}

    static constexpr auto plus() {return impl();}           // log(exp(a) + exp(b))
    static constexpr auto times() {return simd::plus;}          // log(exp(a) * exp(b)) = a + b

    static constexpr auto invert() {return simd::negate;}  // log(1/exp(a)) = -a
    static constexpr auto ldexp() {return first_arg();}

    static constexpr T energy_scale(T const beta) {return -beta / std::log(2.0);}
    static constexpr T boltz(T const energy, T const energy_scale) {return max(energy * energy_scale, zero());}
    
    static constexpr T truncate(T const factor) {return factor <= zero() ? -inf<T>() : factor;}

    static constexpr T free_energy(T const factor, T const energy_scale) {return truncate(factor) * std::log(2.0);}
    static constexpr T as_logarithm(T const factor) {return truncate(factor) * std::log(2.0);}
    static constexpr T as_exponential(T const factor) {return std::exp2(factor);}
    // template <class M, class E, class F>
    // static auto set_element(M &m, E e, F &&rule) {return m = mantissa(fw<F>(rule), -e), False();}

    static void adjust_exponential(T &value, real scale) {value += std::log2(scale);}
};

/// In the MFE ring, + is min() and * is plus()
template <class T>
struct MFE {
    // static constexpr uint value = 1;
    using value_type = T;

    static constexpr T log_zero() {return std::numeric_limits<T>::infinity();}
    static constexpr T zero() {return std::numeric_limits<T>::infinity();}
    static constexpr T one() {return 0;}

    static constexpr auto plus() {return simd::min;}
    static constexpr auto times() {return simd::plus;}

    static constexpr auto plus_eq() {return min_eq;}
    static constexpr auto invert() {return simd::negate;}

    static constexpr auto sum() {return simd::reduce_min;}
    static constexpr auto ldexp() {return first_arg();}

    static constexpr bool valid(T const t) {
        NUPACK_QUICK_REQUIRE(t, >, -inf<T>(), "Invalid MFE element result");
        return true;
    }

    template <class E, class F>
    static auto element_value(bool const &err, F &&rule, E e) {return mantissa(fw<F>(rule), -e);}

    static constexpr T boltz(T const energy, Ignore) {return energy;}
    static constexpr T free_energy(T const factor, Ignore) {return factor;}
    static constexpr T energy_scale(Ignore) {return one();}
    static constexpr T as_logarithm(T const t) {return t;}
};

/******************************************************************************************/

}
