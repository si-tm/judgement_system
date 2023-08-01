/**
 * @brief Compile-time constants; pretty rudimentary and should be swapped out (maybe for boost hana)
 *
 * @file Constants.h
 * @author Mark Fornace
 * @date 2018-05-18
 */
#pragma once
#include "Traits.h"
#include <limits>

namespace nupack {

/******************************************************************************************/

template <int N>
using size_constant = std::integral_constant<std::size_t, N>;

template <int N>
using int_constant = std::integral_constant<int, N>;

/******************************************************************************************/

// template <class T, auto N, class=void>
// struct ConvertConstant {
//     template <class U=T, NUPACK_IF(std::is_same_v<T, U> && can_convert<decltype(N), U>)>
//     constexpr U operator()() const {return static_cast<U>(N);}
// };

// // template <class T, auto N> using ConvertConstant = if_t<derives_from<ConvertConstant<T, N>, std::false_type>,
// //     detail::DefaultConstantConvert<T, N>, ConvertConstant<T, N>>;

// /// compile-time constant that can be explicitly converted to numeric types, generally
// template <auto N>
// struct Constant {
//     static constexpr auto value = N;

//     template <class T, NUPACK_IF(can_construct<decltype(ConvertConstant<T, N>()()), T>)>
//     constexpr operator T() const {return static_cast<T>(ConvertConstant<T, N>()());}

//     /// provide another type that is implicitly convertible in the same way
//     // struct implicit {
//     //     template <class T, NUPACK_IF(can_convert<decltype(ConvertConstant<T, N>()()), T>)>
//     //     constexpr operator T() const {return static_cast<T>(ConvertConstant<T, N>()());}
//     //     constexpr operator Constant() const {return Constant();}

//     constexpr Constant<+N> operator+() const {return {};}
//     constexpr Constant<-N> operator-() const {return {};}

//     // constexpr auto operator*() const {return implicit();}
// };

// template <auto T, auto U>
// constexpr Constant<T + U> operator+(Constant<T>, Constant<U>) {return {};}

// #define NUPACK_TEMP(op) \
// template <class T, auto C> \
// constexpr auto operator op(T const t, Constant<C> i) -> decltype(t op C) {return t op C;} \
// template <class T, auto C> \
// constexpr auto operator op(Constant<C> i, T const t) -> decltype(C op t) {return C op t;}

//    NUPACK_TEMP(+); NUPACK_TEMP(-); NUPACK_TEMP(*); NUPACK_TEMP(/); NUPACK_TEMP(|);
//    NUPACK_TEMP(<<); NUPACK_TEMP(>>); NUPACK_TEMP(<); NUPACK_TEMP(>);
//    NUPACK_TEMP(<=); NUPACK_TEMP(>=); NUPACK_TEMP(==); NUPACK_TEMP(!=);
// #undef NUPACK_TEMP

// /******************************************************************************************/

// struct True : Constant<true> {};
// static constexpr auto const true_c = True{};

// struct False : Constant<false> {};
// static constexpr auto const false_c = False{};

// using Zero = Constant<0>;
// static constexpr auto const zero = Zero{};

template <class T>
struct Always {
    T t;
    template <class ...Ts>
    constexpr T operator()(Ts const &...ts) const {return t;}
};

// struct AlwaysZero {template <class ...Ts> constexpr Zero operator()(Ts const &...) const {return {};}};

// using One = Constant<1>;
// static constexpr auto const one = One{};

// using Two = Constant<2>;
// static constexpr auto const two = Two{};

// using Three = Constant<3>;
// static constexpr auto const three = Three{};


// using Inf = Constant<99>;
// static constexpr auto const inf = Inf{};

// using MinusInf = Constant<-99>;
// static constexpr auto const minus_inf = MinusInf{};

// using Blank = Constant<100>;
// static constexpr auto const blank = Blank{};

// using NaN = Constant<-100>;
// static constexpr auto const nan = NaN{};

// //#define NUPACK_COP(op) template <class M, class N> \
// //constexpr auto operator op(Constant<M> c, Constant<N> d) {return bool_t<decltype(c)::value op decltype(d)::value>{};}
// //NUPACK_COP(<); NUPACK_COP(>); NUPACK_COP(<=); NUPACK_COP(>=); NUPACK_COP(==); NUPACK_COP(!=);
// //#undef NUPACK_COP

// /******************************************************************************************/

// template <class T>
// struct ConvertConstant<T, Inf::value, std::enable_if_t<std::is_integral_v<T>>> {
//     constexpr T operator()() const {return std::numeric_limits<T>::max();}
// };

// template <class T>
// struct ConvertConstant<T, Inf::value, std::enable_if_t<std::is_floating_point_v<T>>> {
//     constexpr T operator()() const {return std::numeric_limits<T>::infinity();}
// };

// template <class T>
// struct ConvertConstant<T, MinusInf::value, std::enable_if_t<std::is_integral_v<T>>> {
//     constexpr T operator()() const {return std::numeric_limits<T>::lowest();}
// };

// template <class T>
// struct ConvertConstant<T, MinusInf::value, std::enable_if_t<std::is_floating_point_v<T>>> {
//     constexpr T operator()() const {return -std::numeric_limits<T>::infinity();}
// };

// template <class T>
// struct ConvertConstant<T, NaN::value, std::enable_if_t<std::is_floating_point_v<T>>> {
//     constexpr T operator()() const {
//         using L = std::numeric_limits<T>;
//         return L::has_quiet_NaN ? L::quiet_NaN() : L::signaling_NaN();
//     }
// };

// static_assert(std::is_convertible_v<Zero, uint>);

using False = std::false_type;
using True = std::true_type;


template <bool B>
using bool_t  = std::conditional_t<B, True, False>;


template <class T, std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
constexpr T inf() {return std::numeric_limits<T>::infinity();}

template <class T, std::enable_if_t<std::is_integral_v<T>, int> = 0>
constexpr T inf() {return std::numeric_limits<T>::max();}

/******************************************************************************************/

}
