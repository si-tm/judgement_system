#pragma once

namespace nupack {

/******************************************************************************************/

// Array descriptor of compile time dimensions
template <std::size_t I, std::size_t ...Is>
struct ParameterDescriptor {
    static constexpr std::size_t begin = I;
    static constexpr std::size_t ndim = sizeof...(Is);
    constexpr std::size_t size() const noexcept {return (1 * ... * Is);}
    constexpr std::size_t back() const noexcept {return begin + size() - 1;}

    constexpr auto start(Ignore={}) const {return begin;}
    constexpr auto stop(Ignore={}) const {return begin + size();}

    template <class ...Js>
    std::size_t index(Ignore, Js ...js) const noexcept {
        static_assert(sizeof...(Js) == sizeof...(Is));
        std::size_t out = begin, stride = 1;
        ((out += js * stride, stride *= Is), ...); // column major
        return out;
    }

    std::size_t safe_index(Ignore, if_t<Is == 0, uint, uint>... is) const {
        if (((is >= Is) || ...)) NUPACK_ERROR("invalid parameter index", std::make_tuple(is...), std::make_tuple(Is...));
        return index(is...);
    }
};

/******************************************************************************************/

template <std::size_t N>
struct ParameterGrid {
    static constexpr std::size_t ndim = N;

    std::uint32_t begin;

    constexpr auto start(Ignore) const {return begin;}
    auto stop(std::uint32_t n) const {return begin + calculate_size(n);}

    template <class ...Js>
    std::size_t index(std::uint32_t n, Js ...js) const noexcept {
        static_assert(sizeof...(Js) == N);
        std::uint32_t out = begin, stride = 1;
        ((out += js * stride, stride *= n), ...); // column major
        return out;
    }

    std::uint32_t calculate_size(std::uint32_t n) const noexcept {
        std::uint32_t out = 1;
        for (std::size_t i = 0; i != N; ++i) out *= n;
        return out;
    }
    void set_length(std::uint32_t &b, std::uint32_t n) noexcept {begin = b; b += calculate_size(n);}

    ParameterGrid const &operator()() const {return *this;}
};

/******************************************************************************************/

template <class V, class I, class T, std::size_t ...Js>
std::size_t array_index(V const &v, I const &i, T const &t, std::index_sequence<Js...>) {return v.index(i, t[Js]...);}

template <class V, class I, class T>
std::size_t array_index(V const &v, I const &i, T const &t) {return array_index(v, i, t, std::make_index_sequence<V::ndim>());}

/******************************************************************************************/

}
