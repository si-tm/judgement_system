/**
 * @brief SIMD functions and Boost.SIMD wrapper for dynamic programs
 *
 * @file SIMD.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "../iteration/Range.h"
#include "../common/Error.h"
#include "../reflect/Print.h"

#include <boost/align/aligned_allocator.hpp>
#include <type_traits>

/******************************************************************************************/

#ifndef NUPACK_NO_SIMD
#    include <simdpp/simd.h>
#endif

// If the eval is not used, it might be faster on some FMA architectures.
// But it gets unsafe due to returning internal references so the file with LIBSIMDPP_SIMD_EXPR_H 
// defined has to be changed so all the members are not references.
// There isn't any real difference I can see on arm64 mac so I opted against modifying it.
// #define NUPACK_SIMD_EVAL(x) x

#define NUPACK_SIMD_EVAL(x) (x).eval()

/******************************************************************************************/

namespace nupack::math {

    template <class T, class From, class SFINAE=void>
    struct Splat {using valid = False;};
    
    template <class T, class From>
    struct Splat<T, From, void_if<std::is_constructible_v<T, From>>> {
        using valid = True;
        static constexpr T from(From t) {return static_cast<T>(t);}
    };

    template <class To, class From, NUPACK_IF(sizeof(To) == sizeof(From))>
    To bit_cast(From const &src) noexcept {
        To dst;
        std::memcpy(&dst, &src, sizeof(To));
        return dst;
    }

    static_assert(std::numeric_limits<float>::max_exponent == 128);
    static_assert(std::numeric_limits<double>::max_exponent == 1024);

    // ldexp: build float of 2^e, then multiply by t. Denormals are not returned.
    // the min max is necessary so that overflow is captured as inf and underflow as 0
    // now adjusted so that overflow will result in a negative number
    // Incurred operations: 1 each of *, max, min, +, <<

    template <class T, NUPACK_IF(is_same<T, float>)>
    float ldexp(T t, std::int32_t e) noexcept {
        using I = std::int32_t;
        return 2.0f * bit_cast<float>(std::max<I>(I(0), std::min<I>(255, e + 126)) << 23) * t;
    }

    template <class T, NUPACK_IF(is_same<T, double>)>
    inline double ldexp(T t, std::int64_t e) noexcept {
        using I = std::int64_t;
        return 2.0 * bit_cast<double>(std::max<I>(I(0), std::min<I>(2047, e + 1022)) << 52) * t;
    }    

#ifndef NUPACK_NO_SIMD
    template <class T>
    using scalar_type = std::decay_t<decltype(simdpp::reduce_add(std::declval<T>()))>;

    /******************************************************************************************/

    template <class T, std::size_t ...Is>
    std::array<scalar_type<T>, sizeof...(Is)> unpack_impl(T const &t, std::index_sequence<Is...>) {return {simdpp::extract<Is>(t)...};}

    template <class T>
    auto unpack_impl(T const &t) {return unpack_impl(t, std::make_index_sequence<T::length>());}
    
    template <class T>
    auto unpack(T const &t) {return unpack_impl(t.eval());}

    /******************************************************************************************/
    
    template <class T, class E, NUPACK_IF(is_same<scalar_type<T>, float>)>
    inline auto ldexp(T const &t, E const &e) noexcept {
        using F = decltype(simdpp::to_float32(t).eval());
        using I = std::int32_t;
        // By excluding one factor of 2, we can make sure it overflows when oversaturated
        return NUPACK_SIMD_EVAL(2.0f * simdpp::bit_cast<F>(simdpp::shift_l<23>(simdpp::max(I(0), simdpp::min(I(255), e + I(126))))) * t);
    }

    template <class T, class E, NUPACK_IF(is_same<scalar_type<T>, double>)>
    inline auto ldexp(T const &t, E const &e) noexcept {
        using F = decltype(simdpp::to_float64(t).eval());
        using I = std::int64_t;
        return NUPACK_SIMD_EVAL(2.0 * simdpp::bit_cast<F>(simdpp::shift_l<52>(simdpp::max(I(0), simdpp::min(I(2047), e + I(1022))))) * t);
    }
#endif

template <int N, class T, NUPACK_IF(is_same<T, float> && N == 0)> 
inline T exp2(T const p) noexcept {return std::exp2(p);}

template <int N, class T, NUPACK_IF(is_same<T, float> && N == 0)>
inline T log2(T const x) noexcept {return std::log2(x);}

#   define NUPACK_TMP(X, EXPR) \
        template <int N, class T, NUPACK_IF(is_same<T, float> && N == X)> \
        inline T exp2(T const p) noexcept { \
            T const clipp = max(-126.0f, p); \
            T const z = clipp - std::floor(clipp); \
            return bit_cast<T>(std::uint32_t(T(1 << 23) * (clipp + EXPR))); \
        }
        NUPACK_TMP(1, 121.2740575f + -27.7280233f / (-4.84252568f + z) + -1.49012907f * z)
        NUPACK_TMP(2, 127.00000714443904f + z*(-0.3070674948344716f + z*(0.24170758463001507f + (0.05166834240891416f + 0.013676819273390137f*z)*z)))
        NUPACK_TMP(3, 126.99999977622763f + z*(-0.3068433142078966f + z*(0.2401320454848741f + z*(0.0558760569512532f + (0.008940786605023437f + 0.0018944130673471293f*z)*z))))
        NUPACK_TMP(4, 127.0f + z*(-0.307004344899854f + z*(0.24156559758227208f + (0.05175227608915405f + 0.013686471228427873f*z)*z)))
        NUPACK_TMP(5, 127.0f + z*(-0.30684641060675133f + z*(0.24014418894498427f + z*(0.055858567707279615f + (0.008948442438529115f + 0.0018952115159583394f*z)*z))))
        NUPACK_TMP(6, 127.0f + z*(-0.30685304015390386f + z*(0.24023037305770092f + z*(0.055480683521299735f + z*(0.009684322168588857f + (0.0012388863515365087f + 0.00021877505477783075f*z)*z)))))
#   undef NUPACK_TMP

#   define NUPACK_TMP(X, EXPR) \
        template <int N, class T, NUPACK_IF(is_same<T, float> && N == X)> \
        inline T log2(T const x) noexcept { \
            std::uint32_t const vx = bit_cast<std::uint32_t>(x); \
            T const z = bit_cast<T>(std::uint32_t((vx & 0x007FFFFF) | 0x3f000000)); \
            return T(vx) * 1.1920928955078125e-7f + EXPR; \
        }
        NUPACK_TMP(1, -1.498030302f * z + -1.72587999f / (0.3520887068f + z) + -124.22551499f)
        NUPACK_TMP(2, z*(6.058254889018033f + z*(-8.326893316539039f + (5.032519894385721f - 1.2669426630966842f*z)*z)) + -128.49708524582465f)
        NUPACK_TMP(3, z*(8.095732190409734f + z*(-13.974814880422908f + z*(12.756970696290203f + z*(-6.4812061032749755f + 1.3904702507142956f*z)))) + -128.7871295909422f)
        NUPACK_TMP(4, z*(6.11967589831446f + z*(-8.430450595367814f + (5.105839491450637f - 1.2849168614940558f*z)*z)) + -128.5101479329032f)
        NUPACK_TMP(5, z*(8.15122574121134f + z*(-14.107055571458751f + z*(12.909337804754145f + z*(-6.5654796298682445f + 1.4081500715664472f*z)))) + -128.79617841620495f)
        NUPACK_TMP(6, z*(10.184113054825591f + z*(-21.218354263516257f + z*(26.01118233219783f + z*(-19.977967611464095f + (8.644315249525148f - 1.6080367062131866f*z)*z)))) + -129.03525205535502f)
#   undef NUPACK_TMP


#ifndef NUPACK_NO_SIMD

#   define NUPACK_TMP(X, EXPR) \
        template <int N, class T, NUPACK_IF(is_same<scalar_type<T>, float> && N == X)> \
        auto log2(T const &x) noexcept { \
            using M = std::decay_t<decltype(x.eval())>; \
            using I = std::decay_t<decltype(simdpp::to_uint32(x).eval())>; \
            auto const vx = simdpp::bit_cast<I>(x).eval(); \
            auto const z = simdpp::bit_cast<M>((vx & std::uint32_t(0x007FFFFF)) | std::uint32_t(0x3f000000)).eval(); \
            return (simdpp::to_float32(vx) * 1.1920928955078125e-7f + EXPR).eval(); \
        }
        NUPACK_TMP(1, -1.498030302f * z + -1.72587999f / (0.3520887068f + z) + -124.22551499f)
        NUPACK_TMP(2, z*(6.058254889018033f + z*(-8.326893316539039f + (5.032519894385721f - 1.2669426630966842f*z)*z)) + -128.49708524582465f)
        NUPACK_TMP(3, z*(8.095732190409734f + z*(-13.974814880422908f + z*(12.756970696290203f + z*(-6.4812061032749755f + 1.3904702507142956f*z)))) + -128.7871295909422f)
        NUPACK_TMP(4, z*(6.11967589831446f + z*(-8.430450595367814f + (5.105839491450637f - 1.2849168614940558f*z)*z)) + -128.5101479329032f)
        NUPACK_TMP(5, z*(8.15122574121134f + z*(-14.107055571458751f + z*(12.909337804754145f + z*(-6.5654796298682445f + 1.4081500715664472f*z)))) + -128.79617841620495f)
        NUPACK_TMP(6, z*(10.184113054825591f + z*(-21.218354263516257f + z*(26.01118233219783f + z*(-19.977967611464095f + (8.644315249525148f - 1.6080367062131866f*z)*z)))) + -129.03525205535502f)
#   undef NUPACK_TMP

#   define NUPACK_TMP(X, EXPR) \
        template <int N, class T, NUPACK_IF(is_same<scalar_type<T>, float> && N == X)> \
        auto exp2(T const &p) noexcept { \
            auto clipp = simdpp::max(-126.0f, p).eval(); \
            auto const z = (clipp - simdpp::floor(clipp)).eval(); \
            return simdpp::bit_cast<decltype(clipp)>(simdpp::to_uint32(float(1 << 23) * (clipp + EXPR))).eval(); \
        }
        NUPACK_TMP(1, 121.2740575f + -27.7280233f / (-4.84252568f + z) + -1.49012907f * z)
        NUPACK_TMP(2, 1.0000071444388348f + z*(-0.3070674948344716f + z*(0.24170758463001507f + (0.05166834240891416f + 0.013676819273390137f*z)*z)))
        NUPACK_TMP(3, 0.9999997762274585f + z*(-0.3068433142078966f + z*(0.2401320454848741f + z*(0.0558760569512532f + (0.008940786605023437f + 0.0018944130673471293f*z)*z))))
        NUPACK_TMP(4, 127.0f + z*(-0.307004344899854f + z*(0.24156559758227208f + (0.05175227608915405f + 0.013686471228427873f*z)*z)))
        NUPACK_TMP(5, 127.0f + z*(-0.30684641060675133f + z*(0.24014418894498427f + z*(0.055858567707279615f + (0.008948442438529115f + 0.0018952115159583394f*z)*z))))
        NUPACK_TMP(6, 127.0f + z*(-0.30685304015390386f + z*(0.24023037305770092f + z*(0.055480683521299735f + z*(0.009684322168588857f + (0.0012388863515365087f + 0.00021877505477783075f*z)*z)))))
#   undef NUPACK_TMP

    template <class T, class F, std::size_t ...Is>
    T apply_impl(T const &t, F const &f, std::index_sequence<Is...>) {return simdpp::make_float(f(simdpp::extract<Is>(t))...);}

    template <class T, class F>
    T apply_impl(T const &t, F const &f) {return apply_impl(t, f, std::make_index_sequence<T::length>());}
    
    template <class T, class F>
    auto apply(T const &t, F const &f) {return apply_impl(t.eval(), f);}

    template <int N, class T, NUPACK_IF(is_same<scalar_type<T>, float> && N == 0)> 
    inline auto exp2(T const p) noexcept {return apply(p, [](auto t) {return std::exp2(t);});}

    template <int N, class T, NUPACK_IF(is_same<scalar_type<T>, float> && N == 0)>
    inline auto log2(T const p) noexcept {return apply(p, [](auto t) {return std::log2(t);});}

#endif

}

/******************************************************************************************/

namespace nupack::simd {

template <class X, NUPACK_IF(std::is_scalar_v<X>)>
decltype(auto) eval_if_simd(X const &x) {return x;}

template <class X, NUPACK_IF(!std::is_scalar_v<X>)>
decltype(auto) eval_if_simd(X const &x) {return x.eval();}

/******************************************************************************************/

/* Chunk is used as an element accessor, i.e. array[Chunk()]. It tells the container
   to return N elements starting at a given position  */
template <int N>
struct Chunk {
    static constexpr auto length = N;
    int value;
    explicit constexpr Chunk(int v) : value(v) {}

    friend auto operator-(Chunk c, int i) noexcept {return Chunk(c.value-i);}
    friend auto operator-(int i, Chunk c) noexcept {return Chunk<-N>(i-c.value);}
    friend auto operator+(int i, Chunk c) noexcept {return Chunk(c.value+i);}
    friend auto operator+(Chunk c, int i) noexcept {return Chunk(c.value+i);}

    NUPACK_REFLECT(Chunk, value); 
};

NUPACK_DEFINE_TEMPLATE(is_chunk, Chunk, int);

/******************************************************************************************/

template <class T, class SFINAE=void>
struct SingleDispatch;

#define NUPACK_TMP(name) \
    struct name##_t { \
        template <class T> \
        auto operator()(T &&t) const noexcept { \
            return SingleDispatch<no_qual<T>>::name(std::forward<T>(t)); \
        } \
    }; \
    static constexpr name##_t name{};

NUPACK_TMP(ifrexp);
NUPACK_TMP(reciprocal);
NUPACK_TMP(negate);
NUPACK_TMP(reduce_sum);
NUPACK_TMP(reduce_min);
NUPACK_TMP(reduce_max);

#undef NUPACK_TMP

template <class T>
struct SingleDispatch<T, std::enable_if_t<std::is_arithmetic_v<T>>> {
    static T reciprocal(T const &t) noexcept {return T(1) / t;}
    static T negate(T const &t) noexcept {return -t;}
    static auto ifrexp(T const &t) noexcept {
        int e;
        T m = std::frexp(t, &e);
        return std::pair<T, int_of_size<sizeof(T)>>(m, e);
    }
    static constexpr T const &reduce_sum(T const &t) noexcept {return t;}
    static constexpr T const &reduce_min(T const &t) noexcept {return t;}
    static constexpr T const &reduce_max(T const &t) noexcept {return t;}
};

/******************************************************************************************/

template <class T, class U, class SFINAE=void>
struct DoubleDispatch;

#define NUPACK_TMP(name) \
    struct name##_t { \
        template <class T, class U> \
        auto operator()(T &&t, U &&u) const noexcept { \
            return DoubleDispatch<no_qual<T>, no_qual<U>>::name(std::forward<T>(t), std::forward<U>(u)); \
        } \
    }; \
    static constexpr name##_t name{};

NUPACK_TMP(load);
NUPACK_TMP(store);
NUPACK_TMP(plus);
NUPACK_TMP(min);
NUPACK_TMP(max);
NUPACK_TMP(multiplies);
NUPACK_TMP(ldexp);
NUPACK_TMP(set_logarithm);
NUPACK_TMP(log2);
NUPACK_TMP(exp2);

#undef NUPACK_TMP

/******************************************************************************************/

template <class T>
struct DoubleDispatch<T, T, std::enable_if_t<std::is_arithmetic_v<T>>> {
    static T plus(T const &t, T const &u) noexcept {return t + u;}
    static T multiplies(T const &t, T const &u) noexcept {return t * u;}
    static T min(T const &t, T const &u) noexcept {return std::min(t, u);}
    static T max(T const &t, T const &u) noexcept {return std::max(t, u);}

    static void set_logarithm(T &t, T const &u) noexcept {t = std::log(u);}
};

template <class T, class U>
struct DoubleDispatch<T, U, std::enable_if_t<std::is_floating_point_v<T> && std::is_floating_point_v<U> && !std::is_same_v<T, U>>> {
    static void set_logarithm(T &t, U const &u) noexcept {
        t = std::log(static_cast<std::common_type_t<T, U>>(u));
    }
};

template <class T>
struct DoubleDispatch<T, int_of_size<sizeof(T)>, std::enable_if_t<std::is_floating_point_v<T>>> {
    static T ldexp(T const &t, int_of_size<sizeof(T)> e) noexcept {return math::ldexp(t, e);}
};

template <int N> struct tag {};

template <int N, class T>
struct DoubleDispatch<tag<N>, T> {
    static decltype(auto) exp2(tag<N>, T const &t) noexcept {return math::exp2<N>(t);}
    static decltype(auto) log2(tag<N>, T const &t) noexcept {return math::log2<N>(t);}
};

/******************************************************************************************/

#ifndef NUPACK_NO_SIMD
template <class T, int N>
using simd_type = if_t<is_same<T, float>, simdpp::float32<N>,
                    if_t<is_same<T, double>, simdpp::float64<N>,
                      if_t<is_same<T, std::int64_t>, simdpp::int64<N>,
                        if_t<is_same<T, std::int32_t>, simdpp::int32<N>, void>>>>;

NUPACK_DETECT(is_simd, decltype(std::declval<T>().eval()));

/******************************************************************************************/

template <class T>
struct SingleDispatch<T, std::enable_if_t<traits::is_simd<T>>> {
    static auto reduce_sum(T const &t) noexcept {return simdpp::reduce_add(t);}
    static auto reduce_min(T const &t) noexcept {return simdpp::reduce_min(t);}
    static auto reduce_max(T const &t) noexcept {return simdpp::reduce_max(t);}
    static auto negate(T const &t) noexcept {return NUPACK_SIMD_EVAL(simdpp::neg(t));}
};

/******************************************************************************************/

template <class T, int N>
struct DoubleDispatch<Chunk<N>, T *, std::enable_if_t<(N >= 1)>> {
    static auto load(Chunk<N> i, T *t) noexcept {return simd_type<std::remove_cv_t<T>, N>(simdpp::load_u(t+i.value));}
};

// Loading reversed chunk includes the pointer t but goes backwards from there
// Should experiment to see if the following permute options are actually faster.

// template <class T>
// struct DoubleDispatch<Chunk<-4>, T *> {
//     using O = simd_type<std::remove_cv_t<T>, 4>;
//     static O load(Ignore, T *t) noexcept {
//         return simdpp::permute4<3,2,1,0>(O(simdpp::load_u(t+(1-4))));
//     }
// };

// template <class T>
// struct DoubleDispatch<Chunk<-2>, T *> {
//     using O = simd_type<std::remove_cv_t<T>, 2>;

//     static O load(Ignore, T *t) noexcept {return simdpp::permute2<1,0>(O(simdpp::load_u(t+(1-2))));}
// };

template <class T>
struct DoubleDispatch<Chunk<-1>, T *> {
    using O = simd_type<std::remove_cv_t<T>, 1>;

    static O load(Chunk<-1> i, T *t) noexcept {return O(simdpp::load_u(t+i));}
};


template <int N, class T>
struct DoubleDispatch<Chunk<N>, T *, std::enable_if_t<N <= -1>> {
    using O = simd_type<std::remove_cv_t<T>, -N>;

    template <std::size_t ...Is>
    static O impl(Chunk<N> i, T *t, std::index_sequence<Is...>) noexcept {
        if constexpr(std::is_integral_v<T> && std::is_signed_v<T>) {
            return simdpp::make_int(t[i.value - int(Is)]...);
        } else if constexpr(std::is_integral_v<T> && !std::is_signed_v<T>) {
            return simdpp::make_uint(t[i.value - int(Is)]...);
        } else if constexpr(std::is_floating_point_v<T>) {
            return simdpp::make_float(t[i.value - int(Is)]...);
        } else {
            static_assert(T::no_load_possible);
        }
        
    }

    static O load(Chunk<N> i, T *t) noexcept {return impl(i, t, std::make_index_sequence<-N>());}
};

template <class X, class T>
struct DoubleDispatch<X, T *, std::enable_if_t<traits::is_simd<X>>> {
    static void store(X const &x, T *t) noexcept {
        simdpp::store(t, x);
    }
};

/******************************************************************************************/

template <class T, class U>
struct DoubleDispatch<T, U, std::enable_if_t<traits::is_simd<T> && traits::is_simd<U>>> {
    static auto plus(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::add(t, u));}
    static auto multiplies(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::mul(t, u));}
    static auto min(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::min(t, u));}
    static auto max(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::max(t, u));}
    static auto ldexp(T const &t, U const &e) noexcept {return math::ldexp(t, e);}
};

/******************************************************************************************/

template <class T, class U>
struct DoubleDispatch<T, U, std::enable_if_t<traits::is_simd<T> && std::is_arithmetic_v<U>>> {
    static auto plus(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::add(t, u));}
    static auto multiplies(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::mul(t, u));}
    static auto min(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::min(t, u));}
    static auto max(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::max(t, u));}
};

template <class T, class U>
struct DoubleDispatch<T, U, std::enable_if_t<traits::is_simd<U> && std::is_arithmetic_v<T>>> {
    static auto plus(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::add(t, u));}
    static auto multiplies(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::mul(t, u));}
    static auto min(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::min(t, u));}
    static auto max(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::max(t, u));}
};

template <class T> struct OptimalSize {};
template <> struct OptimalSize<float> {static constexpr auto value = SIMDPP_FAST_FLOAT32_SIZE;};
template <> struct OptimalSize<double> {static constexpr auto value = SIMDPP_FAST_FLOAT64_SIZE;};
template <> struct OptimalSize<std::int32_t> {static constexpr auto value = SIMDPP_FAST_FLOAT32_SIZE;};
template <> struct OptimalSize<std::int64_t> {static constexpr auto value = SIMDPP_FAST_FLOAT64_SIZE;};

template <class T>
static constexpr auto optimal_size = OptimalSize<T>::value;

template <class T>
static constexpr std::uintptr_t alignment = sizeof(T) * optimal_size<T>;

template <class T>
using allocator = boost::alignment::aligned_allocator<T, alignment<T>>;


#else

template <class T>
static constexpr std::uintptr_t alignment = alignof(T);

template <class T>
using allocator = std::allocator<T>;

#endif

/******************************************************************************************/

// Compile-time zero object. Avoids computation when an operand is known to be 0.
struct Zero {
    template <class T>
    friend T const &operator+(T const &t, Zero) {return t;}

    template <class T>
    friend T const &operator+(Zero, T const &t) {return t;}

    friend Zero operator+(Zero, Zero) {return {};}

    friend False operator<(Zero, Zero) {return {};}

    template <class T>
    friend Zero operator*(T const &t, Zero) {return {};}

    template <class T>
    friend Zero operator*(Zero, T const &t) {return {};}

    friend Zero operator*(Zero, Zero) {return {};}

    template <class T>
    friend T const &operator-(T const &t, Zero) {return t;}

    template <class T>
    friend T operator-(Zero, T const &t) {return -t;}

    friend Zero operator-(Zero, Zero) {return {};}

    Zero operator-() const {return {};}

    constexpr operator std::int32_t() const {return 0;}
    constexpr operator std::int64_t() const {return 0;}
};

NUPACK_UNARY_FUNCTOR(always_zero, Zero());

/******************************************************************************************/

template <class T>
struct DoubleDispatch<T, Zero> {
    static T const &ldexp(T const &t, Zero) noexcept {return t;}
    static T const &plus(T const &t, Zero) noexcept {return t;}
};

template <>
struct DoubleDispatch<Zero, Zero> {
    static constexpr Zero plus(Zero, Zero) noexcept {return {};}
    static constexpr Zero min(Zero, Zero) noexcept {return {};}
    static constexpr Zero max(Zero, Zero) noexcept {return {};}
};

template <>
struct SingleDispatch<Zero> {
    static constexpr Zero negate(Zero) noexcept {return {};}
};

/******************************************************************************************/

/// Perform a map-reduce operation with SIMD, where the output is modified in place
template <class R, class D, class M, class Op>
auto map_reduce(R reduce, D domain, M map, Op op) noexcept {
    auto it = begin_of(domain);
    auto const end = end_of(domain);
    auto ret = map(*it++);

#ifndef NUPACK_NO_SIMD
    constexpr auto Z = optimal_size<decltype(ret)>;
    if (it + Z <= end) {
        auto sum = map(Chunk<Z>(*it)).eval();
        for (it += Z; it + Z <= end; it += Z) sum = reduce(sum, map(Chunk<Z>(*it)));
        ret = reduce(ret, op(sum));
    }
#endif

    for (; it != end; ++it) ret = reduce(ret, map(*it));
    return ret;
}

/******************************************************************************************/

/// Variadic max function
template <class T, class U, class ...Ts>
constexpr auto max_of(T const &t, U const &u, Ts const &...ts) noexcept {
    if constexpr(sizeof...(Ts) == 0) {
        return max(t, u);
    } else {
        return max_of(max(t, u), ts...);
    }
}

/******************************************************************************************/

}

#ifndef NUPACK_NO_SIMD

namespace nupack::math {

template <class T, class From>
struct Splat<T, From, void_if<(simd::traits::is_simd<T> && std::is_arithmetic_v<From>)>> {
    using valid = True;
    static constexpr T from(From t) {return simdpp::splat<T>(t);}
};

}

#endif

