#pragma once
#include "../standard/Vec.h"
#include "../reflect/Reflection.h"

namespace nupack {

/******************************************************************************************/

template <bool Inclusive, class T>
constexpr auto upper_index(T i, T j) noexcept {return j * (Inclusive ? j+1 : j-1) / 2 + i;}

template <class T, bool Inclusive=true>
struct Triangle {
    vec<T> values;
    uint m_length = 0;

    NUPACK_REFLECT(Triangle, values, m_length);

    Triangle() = default;
    explicit Triangle(uint n) : values(n * (Inclusive ? n+1 : n-1) / 2), m_length(n) {}

    auto length() const noexcept {return m_length;}

    T & operator()(uint i, uint j) {
        NUPACK_QUICK_REQUIRE(i, <, length());
        NUPACK_QUICK_REQUIRE(i, <, j + Inclusive);
        return values[upper_index<Inclusive>(i, j)];
    }

    T const &operator()(uint i, uint j) const {return const_cast<Triangle &>(*this)(i, j);}
};

/******************************************************************************************/

template <class T>
struct Square {
    vec<T> values;
    uint m_length = 0;

    Square() = default;
    explicit Square(uint n) : values(n * n), m_length(n) {}

    NUPACK_REFLECT(Square, values, m_length);

    auto length() const noexcept {return m_length;}

    T & operator()(uint i, uint j) {
        NUPACK_QUICK_REQUIRE(i, <, length());
        NUPACK_QUICK_REQUIRE(j, <, length());
        return values[i * m_length + j];
    }

    T const & operator()(uint i, uint j) const {return const_cast<Square const &>(*this)(i, j);}
};

/******************************************************************************************/

}