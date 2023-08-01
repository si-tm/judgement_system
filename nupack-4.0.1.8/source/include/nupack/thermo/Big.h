#pragma once
#include "Tensor.h"

namespace nupack {

template <class T>
struct native_type_t {using type = std::decay_t<T>;};

template <class T>
using native_type = typename native_type_t<T>::type;

/******************************************************************************************/

template <class T>
using Exponent = std::conditional_t<std::is_const_v<T>, int_of_size<sizeof(T)> const, int_of_size<sizeof(T)>>;

template <class T>
struct Big {
    static_assert(!std::is_reference_v<T>);
    static_assert(!std::is_const_v<T>);

    using mantissa_type = T;
    using exponent_type = Exponent<T>;

    constexpr Big(T t) : mantissa(t), exponent(0) {}
    constexpr Big(T t, exponent_type e) : mantissa(t), exponent(e) {}
    constexpr Big() = default;

    template <class U, NUPACK_IF(std::is_floating_point_v<U>)>
    explicit operator U() const noexcept {return std::ldexp(static_cast<U>(mantissa), exponent);}

    mantissa_type mantissa;
    exponent_type exponent;

    // template <class U, NUPACK_IF(std::is_floating_point_v<U>)>
    // explicit operato

    friend std::ostream &operator<<(std::ostream &os, Big<T> const &b) {return os << b.mantissa << 'o' << b.exponent;}
};

template <class T>
struct native_type_t<Big<T>> {using type = std::decay_t<T>;};

template <class T>
struct is_numeric : std::is_arithmetic<T> {};

template <class T>
struct is_numeric<Big<T>> : True {};

/******************************************************************************************/

template <class T>
struct MaxExponent {
    constexpr int operator()() const noexcept {return std::numeric_limits<T>::max_exponent;}
};

template <class T>
struct MaxExponent<Big<T>> {
    constexpr int operator()() const noexcept {return std::numeric_limits<int>::max();}
};

/******************************************************************************************/

NUPACK_TEMPLATE_FALLBACK(Mantissa, typename T::mantissa_type, T);

/******************************************************************************************/

template <class T>
struct BigIterator : TotallyOrdered {
    using value_type = Big<std::decay_t<T>>;

    T *mantissa;
    Exponent<T> *exponent;

    BigIterator() = default;
    BigIterator(T *m, Exponent<T> *e) : mantissa(m), exponent(e) {}

    operator BigIterator<T const>() const {return {mantissa, exponent};}

    template <class I>
    value_type operator[](I const &i) const {return {mantissa[i], exponent[i]};}

    BigIterator & operator++() {++mantissa; ++exponent; return *this;}

    template <class I>
    friend BigIterator operator+(BigIterator const &t, I const &i) {return {t.mantissa + i, t.exponent + i};}

    bool operator<(BigIterator const &o) const {return mantissa < o.mantissa;}
    bool operator==(BigIterator const &o) const {return mantissa == o.mantissa;}
};

/******************************************************************************************/

template <class T>
struct thermo::Slice<BigIterator<T>> {
    using value_type = Big<std::decay_t<T>>;
    Slice<T *> mantissa;
    Slice<Exponent<T> *> exponent;

    Slice(Slice<T *> const &m, Slice<Exponent<T> *> const &e) : mantissa(m), exponent(e) {}

    Slice(BigIterator<T> const &b, BigIterator<T> const &e) : Slice({b.mantissa, e.mantissa}, {b.exponent, e.exponent}) {}

    value_type operator()(Index i) const {return {mantissa(i), exponent(i)};}
    
    Slice operator()(span i) const {return {mantissa(i), exponent(i)};}

    void set(Index i, Big<T> const &t) {mantissa(i) = t.mantissa; exponent(i) = t.exponent;}

    void fill(Big<T> const &t) {mantissa.fill(t.mantissa); exponent.fill(t.exponent);}

    auto size() const {return mantissa.size();}
};

/******************************************************************************************/

template <class T>
struct thermo::SharedArray<Big<T>> {
    using value_type = Big<T>;
    using iterator = BigIterator<T>;
    using const_iterator = BigIterator<T const>;

    SharedArray<typename Big<T>::mantissa_type> mantissa;
    SharedArray<typename Big<T>::exponent_type> exponent;

    bool allocate(std::size_t n) {
        if (size() == n) return false;
        mantissa.allocate(n);
        exponent.allocate(n);
        return true;
    }

    void clear() noexcept {mantissa.clear(); exponent.clear();}

    void check() const {NUPACK_QUICK_REQUIRE(mantissa.size(), ==, exponent.size());}

    iterator begin() {check(); return {mantissa.begin(), exponent.begin()};}
    iterator end() {check(); return {mantissa.end(), exponent.end()};}

    const_iterator begin() const {check(); return {mantissa.begin(), exponent.begin()};}
    const_iterator end() const {check(); return {mantissa.end(), exponent.end()};}

    auto data() {check(); return begin();}
    auto data() const {check(); return begin();}

    auto size() const {check(); return mantissa.size();}

    explicit operator bool() const {check(); return bool(mantissa);}

    void fill(Big<T> const &t) {
        mantissa.fill(t.mantissa);
        exponent.fill(t.exponent);
    }
};

template <class From, class To>
struct thermo::ConvertArray<From, Big<To>, std::enable_if_t<std::is_scalar_v<From>>> {
    template <class Converter>
    static void copy(SharedArray<From> const &t, SharedArray<Big<To>> &u, Converter const &c) {
        zip(t, u.mantissa, u.exponent, c);
    }
};


template <class From, class To>
struct thermo::ConvertArray<Big<From>, To, std::enable_if_t<std::is_scalar_v<To>>> {
    template <class Converter>
    static void copy(SharedArray<Big<From>> const &t, SharedArray<To> &u, Converter const &c) noexcept {
        zip(t.mantissa, t.exponent, u, c);
    }
};

template <class From, class To>
struct thermo::ConvertArray<Big<From>, Big<To>> {
    template <class Converter>
    static void copy(SharedArray<From> const &t, SharedArray<To> &u, Converter const &c) {
        zip(t.mantissa, t.exponent, u.mantissa, u.exponent, c);
    }
};

/******************************************************************************************/

NUPACK_DETECT(has_exponent, decltype(std::declval<T>().exponent));

template <class T, class I, NUPACK_IF(traits::has_exponent<T>)>
decltype(auto) mantissa_at(T const &t, I const &i) {return t.mantissa(i);}

template <class T, class I, NUPACK_IF(!traits::has_exponent<T>)>
decltype(auto) mantissa_at(T const &t, I const &i) {return t(i);}

template <class T, class I, NUPACK_IF(traits::has_exponent<T>)>
decltype(auto) exponent_at(T const &t, I const &i) {return t.exponent(i);}

template <class T, NUPACK_IF(!traits::has_exponent<T>)>
constexpr simd::Zero exponent_at(T const &t, Ignore) {return {};}

/******************************************************************************************/

template <class T, NUPACK_IF(traits::has_exponent<T>)>
decltype(auto) mantissa(T const &t) {return t.mantissa;}

template <class T, NUPACK_IF(!traits::has_exponent<T>)>
decltype(auto) mantissa(T const &t) {return t;}

template <class T, NUPACK_IF(traits::has_exponent<T>)>
decltype(auto) exponent(T const &t) {return t.exponent;}

template <class T, NUPACK_IF(!traits::has_exponent<T>)>
constexpr simd::Zero exponent(T const &t) {return {};}

/******************************************************************************************/

}

namespace nupack::simd {

template <class T>
struct SingleDispatch<Big<T>> {
    static Big<T> exp(Big<T> const &t) noexcept {return {std::exp(T(t)), 0};}
    static T log(Big<T> const &t) noexcept {return std::log(t.mantissa) + t.exponent * LogOf2;}
    static Big<T> reciprocal(Big<T> const &t) noexcept {return {1/t.mantissa, -t.exponent};}
};


/******************************************************************************************/

template <class T, class U>
struct DoubleDispatch<Big<T>, U> {
    static_assert(std::is_floating_point_v<U>);
    static Big<T> multiplies(Big<T> const &t, U const &u) noexcept {return Big<T>(t.mantissa * u, t.exponent);}
};

template <class T, class U>
struct DoubleDispatch<U, Big<T>> {
    static_assert(std::is_floating_point_v<U>);
    static Big<T> multiplies(U const &u, Big<T> const &t) noexcept {return Big<T>(u * t.mantissa, t.exponent);}
    static void set_logarithm(U &u, Big<T> const &t) noexcept {
        u = std::log(static_cast<std::common_type_t<T, U>>(t.mantissa)) + t.exponent * LogOf2; // use higher precision type
    }
};

template <class T, class U>
struct DoubleDispatch<Big<T>, Big<U>> {
    static Big<T> multiplies(Big<T> const &t, Big<U> const &u) noexcept {return Big<T>(t.mantissa * u.mantissa, t.exponent + u.exponent);}

    static Big<std::remove_cv_t<T>> plus(Big<T> const &t, Big<U> const &u) noexcept {
        if (t.exponent < u.exponent)
             return Big<std::remove_cv_t<T>>(u.mantissa + std::ldexp(t.mantissa, t.exponent-u.exponent), u.exponent);
        else return Big<std::remove_cv_t<T>>(t.mantissa + std::ldexp(u.mantissa, u.exponent-t.exponent), t.exponent);
    }
};

template <class M, class E>
struct BigChunk {
    M mantissa;
    E exponent;
};

template <class M, class E> 
BigChunk(M, E) -> BigChunk<M, E>;

template <class T, int N>
struct DoubleDispatch<Chunk<N>, BigIterator<T>> {
    static auto load(Chunk<N> i, BigIterator<T> const &t) noexcept {
        return BigChunk{simd::load(i, t.mantissa), simd::load(i, t.exponent)};
    }
};

/******************************************************************************************/

// How many elements to go until the next aligned position. 
template <class T>
int next_alignment(T const *t) {
    constexpr std::uintptr_t a = alignment<T>;
    static_assert(0 == (a & (a-1)), "should be power of 2");
    std::uintptr_t const p = reinterpret_cast<std::uintptr_t>(t);
    std::uintptr_t const aligned = p & ~(a - 1); // erase all the lesser bits for the previous aligned position
    return (aligned == p) ? 0 : (aligned + a - p) / sizeof(T);
}

/******************************************************************************************/

/// Perform a map operation with SIMD, where the output is modified in place
// With SIMD in simdpp the store has to be aligned, so some of the first and last elements are handled with non-SIMD
template <bool SIMD=true, class O, class F>
void map(O &&out, int i, int stop, F &&f) noexcept {

    if constexpr(std::is_scalar_v<value_type_of<O>>) {
        auto const m = begin_of(out);
#       ifndef NUPACK_NO_SIMD
            if constexpr(SIMD) {
                auto const first = std::min(i + next_alignment(std::addressof(*(m + i))), stop);
                for (; i < first; ++i) 
                    *(m + i) = f(i);

                constexpr auto Z = simd::optimal_size<value_type_of<O>>;
                for (; i + Z <= stop; i += Z)
                    simd::store(f(Chunk<Z>(i)), std::addressof(*(m + i)));
            }

#       endif
        for (; i < stop; ++i) {
            *(m + i) = f(i);
        }
    } else {
        auto const m = begin_of(out.mantissa);
        auto const e = begin_of(out.exponent);
#       ifndef NUPACK_NO_SIMD
            if constexpr(SIMD) {
                constexpr auto Z = simd::optimal_size<value_type_of<decltype(m)>>;
                static_assert(!is_same<value_type_of<decltype(m)>, double>, "Cannot use Big<real64> with simdpp right now");
                auto const first = std::min(i + next_alignment(std::addressof(*(m + i))), stop);
                for (; i < first; ++i) 
                    std::tie(*(m + i), *(e + i)) = f(i);
                
                for (; i + Z <= stop; i += Z) {
                    auto p = f(Chunk<Z>(i));
                    simd::store(std::move(p.first), std::addressof(*(m + i)));
                    simd::store(std::move(p.second), std::addressof(*(e + i)));
                }
            }
#       endif
        for (; i < stop; ++i) {
            std::tie(*(m + i), *(e + i)) = f(i);
        }
    }
}

template <bool SIMD=true, class O, class F>
void map(O &&out, span s, F &&f) noexcept {
    return map<SIMD>(fw<O>(out), s.start(), s.stop(), fw<F>(f));
}

/******************************************************************************************/

}
