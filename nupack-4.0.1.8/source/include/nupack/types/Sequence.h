/** \file Sequence.h
 * @brief Contains base definitions, enums, Sequence class, SubsequenceList class to hold Sequences
 */
#pragma once
#include <string>
#include <string_view>

#include "Base.h"
#include "../iteration/Search.h"
#include "../iteration/View.h"
#include "../iteration/Patterns.h"
#include "../types/IO.h"

namespace nupack {

/******************************************************************************************/

struct SingleStrand : std::false_type {};
struct MultiStrand : std::true_type {};

/******************************************************************************************/

using Subsequence = View<Base const *>;
using Subdomain = View<Wildcard const *>;

string raw_sequence_string(Subsequence, std::size_t n=0); 
string raw_domain_string(Subdomain, std::size_t n=0);

inline std::ostream& operator<<(std::ostream &os, Subsequence const &s) {return os << raw_sequence_string(s);}
inline std::ostream& operator<<(std::ostream &os, Subdomain const &s) {return os << raw_domain_string(s);}

/******************************************************************************************/

/// BaseIter was previously Sequence::const_iterator, but this is now elided since they're essentially equivalent
using BaseIter = Base const *;
static_assert(is_dumb_ptr<BaseIter>, "");

/******************************************************************************************/

template <class T>
struct SequenceData : CompareByKey<SequenceData<T>> {
    using iterator = T const*;
    using value_type = T;
    using const_iterator = T const*;

    using difference_type = std::int32_t;
    using size_type = std::uint32_t;
    using ID = std::int32_t;

    std::shared_ptr<T const> ptr;
    std::uint32_t length = 0;
    ID id = 0;

    static std::shared_ptr<T> allocate(std::uint32_t n) {return std::shared_ptr<T>(new T[n], [](T *t) {delete[] t;});}

    SequenceData() = default;
    void clear() {ptr.reset(); length = 0; id = 0;}

    template <class F>
    SequenceData(size_type n, F &&f, ID i=0) : length(n), id(i) {
        auto p = allocate(length);
        f(p.get(), n);
        ptr = std::move(p);
    }

    template <class V, NUPACK_IF(!std::is_base_of_v<SequenceData, V> && std::is_convertible_v<value_type_of<V>, T>)>
    SequenceData(V const &v, ID i=0) : length(std::size(v)), id(i) {
        auto p = allocate(length);
        std::copy(std::begin(v), std::end(v), p.get());
        ptr = std::move(p); 
    }

    SequenceData(SequenceData const &) = default;
    SequenceData& operator=(SequenceData const &) = default;

    SequenceData(SequenceData &&s) noexcept 
        : ptr(std::move(s.ptr)), length(std::exchange(s.length, 0)), id(std::exchange(s.id, 0)) {}

    SequenceData& operator=(SequenceData &&s) noexcept {
        ptr = std::move(s.ptr);
        length = std::exchange(s.length, 0);
        id = std::exchange(s.id, 0);
        return *this;
    }
    
    T const* data() const {return ptr.get();}
    std::uint32_t size() const {return length;}

    bool empty() const {return size() == 0;}
    auto begin() const noexcept {return data();}
    auto end() const noexcept {return data() + size();}
    
    T const &operator[](size_type i) const noexcept {return begin()[i];}

    T const &front() const noexcept {return *begin();}
    T const &back() const noexcept {return *std::prev(end());}

    void slice(std::size_t start, std::size_t stop) noexcept {
        T const *data = ptr.get() + start;
        ptr = std::shared_ptr<T const>(std::move(ptr), data);
        length = stop - start;
    }

    void pop_back() noexcept {slice(0, length-1);}
    void pop_front() noexcept {slice(1, length);}

    template <class It>
    void insert(iterator i, It b, It e) {
        if (b == e) return;
        auto const n = std::distance(b, e);
        auto p = allocate(length + n);
        auto it = std::copy(this->begin(), i, p.get());
        it = std::copy(b, e, it);
        std::copy(i, this->end(), it);
        ptr = std::move(p);
        length += n;
    }

    auto compare_key() const noexcept {return std::make_pair(id, view(*this));}
};

/******************************************************************************************/

// Immutable sequence
struct Sequence : SequenceData<Base> {
    using base_type = SequenceData<Base>;
    using base_type::base_type;
    
    Sequence & operator+=(Subsequence const &);
    friend Sequence operator+(Sequence a, Subsequence const &b) {a += b; return a;}

    std::string save_repr() const {return raw_sequence_string(*this, 3);}
    void load_repr(std::string_view);
};

static_assert(std::is_convertible_v<Subsequence, Sequence>);

Sequence operator "" _4(char const *, std::size_t);

/******************************************************************************************/

struct Domain : SequenceData<Wildcard> {
    using base_type = SequenceData<Wildcard>;
    using base_type::base_type;

    Domain & operator+=(Subdomain const &);
    friend Domain operator+(Domain a, Subdomain const &b) {a += b; return a;}

    template <class V>
    static Domain from_sequence(V const &v) {return Domain(indirect_view(v, Wildcard::from_base));}

    std::string save_repr() const {return raw_domain_string(*this, 3);}
    void load_repr(std::string_view);
};

static_assert(std::is_convertible_v<Subdomain, Domain>);

/******************************************************************************************/

struct SeqDistance {
    template <class A, class B>
    constexpr std::size_t operator()(A const &a, B const &b) const noexcept {
        auto const min = std::min(std::size(a), std::size(b));
        auto const max = std::max(std::size(a), std::size(b));
        return max - std::inner_product(std::begin(a), std::begin(a) + min, 
            std::begin(b), std::size_t(0), std::plus<std::size_t>(), are_compatible);
    }
};

static constexpr SeqDistance seq_distance{};

/******************************************************************************************/

using DomainList = vec<Domain, 4>;
using SubsequenceList = vec<Subsequence, 8>;
using SequenceList = vec<Sequence, 4>;

/******************************************************************************************/

template <class V, NUPACK_IF(is_same<value_type_of<V>, Base, Wildcard>)>
std::size_t nt(V const &v) {return len(v);}

template <class V, NUPACK_IF(is_same<value_type_of<value_type_of<V>>, Base, Wildcard>)>
std::size_t nt(V const &v) {return sum(v, len);}

template <class V, NUPACK_IF(is_same<value_type_of<value_type_of<V>>, Base, Wildcard>)>
std::size_t n_strands(V const &v) {return len(v);}

// template <class Out=vec<bool>, class S=Sequence>
// auto one_hot_sequence(S const &sequence);
 // {
 //    Out out(Base::all().size() * len(sequence));
 //    auto it = begin_of(out);
 //    for (auto &&c : sequence) for (auto b : Base::all()) it++ == b;
 //    return out;
// }

/******************************************************************************************/

// template <class RNG=decltype(StaticRNG) &>
// Sequence sample(Sequence s, RNG &&rng=StaticRNG);
// {
//     for (auto &b : s) b = b.sample(rng);
//     return s;
// }

// template <class V>
// inline bool is_palindromic(V const &seq) {return seq == seq.reverse_complement();}

// template <class RNG=decltype(StaticRNG) &>
// Sequence random_sequence(iseq n, real gc=0.5, RNG &&gen=StaticRNG) {
//     auto dist = Base::distribution(gc);
//     Sequence out; out.reserve(n);
//     for (iseq i=0; i != n; ++i) out.push_back(Base::from_index(dist(gen)));
//     return out;
// }

/******************************************************************************************/

template <class T> 
struct memory::impl<T, std::enable_if_t<is_same<T, Sequence, Domain>>> {
    auto operator()(T const &s) const {return sizeof(T) + sizeof(typename T::value_type) * s.size() / s.ptr.use_count();}
    void erase(T &s) {T s_; swap(s, s_);}
};

/******************************************************************************************/

using Nick = int;
static constexpr Nick const NoNick = -1;

/// Find sequence index of nick: the index is to the sequence after the nick
template <class V>
Nick find_nick(V const &v) {
    auto b = std::begin(v);
    for (Nick i = 0; b != std::end(v); ++b, ++i)
        if (front(*b) == Base::null()) return i;
    return NoNick;
}

/******************************************************************************************/

/// Convert a single string e.g. "ACTGTA+ACTGAT" into a vector of strings
vec<string> split_sequence_string(string_view s);

/******************************************************************************************/

/**
 * @brief Given a vector of views, split it into 3 at bb of sequence b and ee of sequence e
 * @return std::pair<V, V> representing {1 + 3, 2}
 */
template <class V>
std::pair<V, V> split_midway(V const &v, const_iterator_of<V> b, const_iterator_of<V> e, BaseIter bb, BaseIter ee) {
    NUPACK_ASSERT(e >= b, e - b);
    // v1 should get len(v) - e + b + 1
    V v1; v1.assign(begin_of(v), b + 1);
    v1.back().set_end(bb + 1);
    extend(v1, e, end_of(v));
    v1[b + 1 - begin_of(v)].set_begin(ee);
    // v2 should get e + 1 - b sequences
    V v2; v2.assign(b, e + 1);
    v2.front().set_begin(bb);
    v2.back().set_end(ee + 1);
    return std::make_pair(v1, v2);
}

/******************************************************************************************/

/// Return new loop sequences after a base pair deletion
template <class V> V merged_seqs(V const &p, V const &k, size_type_of<V> pk, size_type_of<V> kp) {
    auto kpm = (!kp ? len(p) : kp) - 1;
    auto pkm = (!pk ? len(k) : pk) - 1;

    V ret; ret.reserve(len(p) + len(k) - 2);
    circular_cat(ret, p, begin_of(p) + kp, begin_of(p) + kpm);
    ret.front().set_begin(begin_of(k[pkm]));
    if (len(k) > 1) {ret.emplace_back(k[pk]); back(ret).set_begin(begin_of(p[kpm]));}
    else ret.front().set_begin(begin_of(p[kpm]));
    circular_cat(ret, k, begin_of(k) + pk + 1, begin_of(k) + pkm);
    return ret;
}

/******************************************************************************************/

/// Return new loop sequences after a dissociation event
template <class V> std::pair<V, V> get_split_seqs(V const &pseqs, V const &kseqs,
    size_type_of<V> pnick, size_type_of<V> knick, size_type_of<V> pk, size_type_of<V> kp) {
    V new_pseqs, new_kseqs;

    circular_cat(new_pseqs, pseqs, begin_of(pseqs) + pnick, begin_of(pseqs) + kp);
    new_pseqs.back().set_end(end_of(kseqs[pk]));
    circular_cat(new_pseqs, kseqs, begin_of(kseqs) + pk + 1, begin_of(kseqs) + knick);

    circular_cat(new_kseqs, kseqs, begin_of(kseqs) + knick, begin_of(kseqs) + pk);
    new_kseqs.back().set_end(end_of(pseqs[kp]));
    circular_cat(new_kseqs, pseqs, begin_of(pseqs) + kp + 1, begin_of(pseqs) + pnick);

    return std::make_pair(std::move(new_pseqs), std::move(new_kseqs));
}

/******************************************************************************************/

using Edge = int;
static constexpr Edge Ether = -1;
using EdgeList = vec<Edge>;

/******************************************************************************************/

}

namespace std {
    template <> struct hash<nupack::Sequence> {
        size_t operator()(nupack::Sequence const &s) const noexcept;
    };
    template <> struct hash<nupack::Domain> {
        size_t operator()(nupack::Domain const &s) const noexcept;
    };
}
