/** \file Sequence.h
 * @brief Contains Base class for nucleotides
 */
#pragma once
#include <iostream>

#include "../common/Random.h"
#include "../standard/Array.h"
#include "Matrix.h"
#include "../iteration/Search.h"

namespace nupack {

/******************************************************************************************/

using BaseIndex = std::uint8_t;

struct Base {
    static constexpr unsigned capacity = 64; // compatible with WildcardIndex below
    BaseIndex value;

    Base() = default;

    template <class T, NUPACK_IF(std::is_same_v<T, BaseIndex>)>
    constexpr explicit Base(T t) noexcept : value(t) {}

    struct from_index_t {
        constexpr Base operator()(BaseIndex value) const noexcept {return Base(value);}
    };

    static constexpr auto const from_index = from_index_t{};

    constexpr BaseIndex operator+() const noexcept {return value;}

    constexpr bool operator==(Base b) const noexcept {return value == b.value;}
    constexpr bool operator!=(Base b) const noexcept {return value != b.value;}
    constexpr bool operator< (Base b) const noexcept {return value <  b.value;}
    constexpr bool operator> (Base b) const noexcept {return value >  b.value;}
    constexpr bool operator<=(Base b) const noexcept {return value <= b.value;}
    constexpr bool operator>=(Base b) const noexcept {return value >= b.value;}

    static constexpr Base null() noexcept {return from_index(capacity);}
    static constexpr Base invalid() noexcept {return from_index(capacity+1);}

    constexpr bool is_null() const noexcept {return *this == null();}

    char raw_char() const noexcept;
    static Base from_raw_char(char c) noexcept;

    friend std::ostream & operator<<(std::ostream &os, Base b) {return os << b.raw_char();}

    BaseIndex save_repr() const noexcept {return value;}
    void load_repr(BaseIndex c) noexcept {value = c;}
};

static_assert(sizeof(Base) == sizeof(BaseIndex), "Should be same size");
static_assert(std::is_pod<Base>::value, "Should be POD type");

/******************************************************************************************/

using WildcardIndex = std::uint64_t;

struct Wildcard {
    WildcardIndex value;

    template <class T, NUPACK_IF(std::is_same_v<T, WildcardIndex>)>
    constexpr explicit Wildcard(T t) noexcept : value(t) {}
    
    Wildcard() = default;

    struct from_index_t {
        constexpr Wildcard operator()(WildcardIndex i) const noexcept {return Wildcard(i);}
    };

    static constexpr auto const from_index = from_index_t{};
    
    struct from_base_t {
        constexpr Wildcard operator()(Base b) const noexcept {return Wildcard(WildcardIndex(1) << (+b));}
    };

    static constexpr auto const from_base = from_base_t{};

    struct from_mask_t {
        template <class V>
        constexpr Wildcard operator()(V const &v) const noexcept {
            WildcardIndex w = 0;
            izip(v, [&w](auto i, auto const &b) {if (b) w |= (WildcardIndex(1) << i);});
            return from_index(w);
        }
    };

    static constexpr auto const from_mask = from_mask_t{};

    struct from_indices_t {
        template <class V>
        constexpr Wildcard operator()(V const &v) const noexcept {
            WildcardIndex w = 0;
            for (WildcardIndex const &i : v) w |= (WildcardIndex(1) << i);
            return from_index(w);
        }
    };

    static constexpr auto const from_indices = from_indices_t{};


    static constexpr Wildcard null() noexcept {return from_index(0);}

    // use std::bitset if you want the boolean mask.

    // indices of allowed bases
    template <class Out=vec<BaseIndex>>
    Out indices() const {
        Out out;
        auto v = value;
        for (BaseIndex i = 0; v; ++i, v >>= 1) if (v & 1) out.emplace_back(i); 
        return out;
    }

    // Return first possible base, or else null()
    constexpr Base first() const noexcept {
        auto v = value;
        for (BaseIndex i = 0; v; ++i, v >>= 1) if (v & 1) return Base::from_index(i);
        return Base::null();
    }

    constexpr explicit operator bool() const noexcept {return value;}

    // Return if the base satisfies the wildcard
    constexpr bool includes(Base b) const noexcept {return value & (1 << b.value);}

    // Return if all bases in b satisfy this
    constexpr bool includes(Wildcard b) const noexcept {return (value | b.value) == value;}

    // Return if there is exactly one possible value
    constexpr bool is_determined() const noexcept {return value && !(value & (value - 1));}

    constexpr WildcardIndex operator+() const noexcept {return value;}
    constexpr Wildcard operator~() const noexcept {return from_index(~value);}
    constexpr bool operator==(Wildcard b) const noexcept {return value == b.value;}
    constexpr bool operator!=(Wildcard b) const noexcept {return value != b.value;}
    constexpr bool operator< (Wildcard b) const noexcept {return value <  b.value;}
    constexpr bool operator> (Wildcard b) const noexcept {return value >  b.value;}
    constexpr bool operator<=(Wildcard b) const noexcept {return value <= b.value;}
    constexpr bool operator>=(Wildcard b) const noexcept {return value >= b.value;}
    
    friend constexpr Wildcard operator^(Wildcard a, Wildcard b) noexcept {return Wildcard(a.value ^ b.value);}
    friend constexpr Wildcard operator|(Wildcard a, Wildcard b) noexcept {return Wildcard(a.value | b.value);}
    friend constexpr Wildcard operator&(Wildcard a, Wildcard b) noexcept {return Wildcard(a.value & b.value);}

    char raw_char() const noexcept;
    static Wildcard from_raw_char(char c) noexcept;

    friend std::ostream & operator<<(std::ostream &os, Wildcard b);

    WildcardIndex save_repr() const noexcept {return value;}
    void load_repr(WildcardIndex c) noexcept {value = c;}
};


static_assert(sizeof(Wildcard) == sizeof(WildcardIndex), "Should be same size");
static_assert(std::is_pod<Wildcard>::value, "Should be POD type");

NUPACK_UNARY_FUNCTOR(is_determined, t.is_determined());

/******************************************************************************************/

struct AreCompatible {
    constexpr bool operator()(Base a, Base b) const noexcept {return a == b;}
    constexpr bool operator()(Wildcard a, Base b) const noexcept {return a.includes(b);}
    constexpr bool operator()(Base a, Wildcard b) const noexcept {return b.includes(a);}
    constexpr bool operator()(Wildcard a, Wildcard b) const noexcept {return a.value & b.value;}

    template <class A, class B>
    constexpr bool operator()(A const &a, B const &b) const noexcept {
        return std::equal(a.begin(), a.end(), b.begin(), b.end(), *this);
    }
};

static constexpr AreCompatible are_compatible{};

/******************************************************************************************/

struct BasePairing {
    using is_member_ordered = True;
    
    small_vec<bool, 16> pairing, closing;
    small_vec<small_vec<Base, 4>, 4> possible_pairs;

    std::uint32_t length;
    BasePairing() = default;
    BasePairing(std::uint32_t n) : pairing(n * n), closing(n * n), possible_pairs(n), length(n) {}

    NUPACK_REFLECT(BasePairing, pairing, closing, possible_pairs);

    constexpr uint turn() const {return 3;}
    constexpr uint multi_gap() const {return turn() + turn() + 2 + 2;} // 10. min multiloop end-end distance in single strand

    auto pairs(Base b) const noexcept {return view(possible_pairs[+b]);}
    auto all() const noexcept {return indirect_view(range<BaseIndex>(length), Base::from_index);}

    bool can_close(Base b, Base c) const noexcept {return closing[length * b.value + c.value];} 
    bool can_pair(Base b, Base c) const noexcept {return pairing[length * b.value + c.value];} 

    bool operator()(Base b, Base c) const {return can_pair(b, c);}

    template <class Iter>
    bool operator()(bool diff_strand, Iter b, Iter c) const {
        if (!diff_strand) NUPACK_QUICK_REQUIRE(b, <=, c, "same-strand bases should be ordered for this function");
        return (diff_strand || b + turn() < c) && can_pair(*b, *c);
    }

    template <class V>
    bool check_loop(V const &v) const {
        bool ok = true;
        for_circularly_adjacent(v, [&](auto const &s1, auto const &s2) {
            Base const b = s1.back(), c = s2.front();
            ok = ok && ((b == Base::null() && c == Base::null()) || (*this)(b, c));
        });
        return ok;
    }
};

/******************************************************************************************/

template <> struct is_character_t<Base> : std::true_type {};
template <> struct is_character_t<Wildcard> : std::true_type {};

}

namespace std {
    template <> struct hash<nupack::Base> {
        size_t operator()(nupack::Base b) const noexcept {return hash<typename nupack::BaseIndex>()(b.value);}
    };
    template <> struct hash<nupack::Wildcard> {
        size_t operator()(nupack::Wildcard b) const noexcept {return hash<typename nupack::WildcardIndex>()(b.value);}
    };
}
