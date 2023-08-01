#pragma once
#include "Sequence.h"
#include "../standard/Map.h"
#include "nupack/algorithms/Operators.h"

namespace nupack {

/******************************************************************************************/

struct MaterialData : MemberOrdered {
    char prefix;
    std::array<Base, CharCapacity> bases; // Base for each input character
    std::array<Wildcard, CharCapacity> wildcards; // Wildcard for each input character
    string name;
    NUPACK_REFLECT(MaterialData, prefix, bases, wildcards, name);
};

/******************************************************************************************/

struct AlphabetData : CompareByKey<AlphabetData> {
    std::array<Base, CharCapacity> bases; // canonical base for each letter
    std::array<char, Base::capacity+1> letters; // canonical name of each base
    std::array<Base, Base::capacity+1> complements; // complement for each base
    vec<MaterialData> materials;
    vec<std::string> base_names; // names of each base
    std::unordered_map<Wildcard, char> wildcard_letters; // assigned names for specific wildcards

    json save_repr() const;
    void load_repr(json const &);
    std::size_t length() const {return base_names.size();}

    auto const &compare_key() const {return base_names;}

    NUPACK_REFLECT(AlphabetData, bases, letters, complements, wildcard_letters, base_names, materials);

};

/******************************************************************************************/

struct Alphabet : TotallyOrdered {
    std::shared_ptr<AlphabetData> data;
    BaseIndex m_length = 0; // number of possible bases

    Alphabet() = default;
    Alphabet(std::shared_ptr<AlphabetData> d) : data(std::move(d)), m_length(data->length()) {}

    NUPACK_REFLECT(Alphabet, data);

    AlphabetData const &get() const {NUPACK_QUICK_ASSERT(data); return *data;}

    char operator()(Base b) const noexcept {return get().letters[+b];}
    char operator()(Wildcard b) const noexcept;

    string to_string(View<Base const *>, std::size_t n=0) const;
    string to_string(View<Wildcard const *>, std::size_t n=0) const;

    auto to_string(std::size_t n=0) const {return [a=*this, n](auto const &x) {return a.to_string(x, n);};}

    template <class V>
    string to_delimited_string(V const &v, char delim='+') const {
        string o;
        o.reserve(sum(v, len));
        for (auto const &s : v) {if (!o.empty()) o.push_back(delim); o += to_string(s);}
        return o;
    }

    auto const &material(uint i) const {
        NUPACK_QUICK_ASSERT(data && i < len(data->materials), "material index out of range");
        return data->materials[i];
    }

    Base first_base(char) const;
    Wildcard first_wildcard(char) const;
    Sequence sequence(std::string_view) const;
    Domain domain(std::string_view) const;

    /// Convert multiple strings e.g. ["ACTGTA", "ACTGAT"] into a vector of strands
    template <class V=SequenceList, class S=vec<string>, NUPACK_IF(!can_construct<string_view, S>)>
    V sequences(S const &strs) const {return vmap<V>(strs, [&](auto const &s) {return sequence(s);});}

    /// Convert a single string e.g. "ACTGTA+ACTGAT" into a vector of strands
    template <class V=SequenceList>
    V sequences(string_view s) const {return sequences<V>(split_sequence_string(s));}

    /// Convert multiple strings e.g. ["ACTGTA", "ACTGAT"] into a vector of domains
    template <class V=DomainList, class S=vec<string>, NUPACK_IF(!can_construct<string_view, S>)>
    V domains(S const &strs) const {return vmap<V>(strs, [&](auto const &s) {return domain(s);});}

    /// Convert a single string e.g. "ACTGTA+ACTGAT" into a vector of domains
    template <class V=DomainList>
    V domains(string_view s) const {return domains<V>(split_sequence_string(s));}

    Base complement(Base b) const noexcept {return get().complements[b.value];}

    // union of the complement of each possible base in b
    Wildcard complement(Wildcard b) const noexcept;

    template <class V>
    V reverse_complement(V const &v) const {return indirect_view(reversed(v), [&](auto b) {return complement(b);});}

    BaseIndex length() const noexcept {return m_length;}

    auto all() const noexcept {return indirect_view(range(length()), Base::from_index);}

    bool operator==(Alphabet const &b) const noexcept {
        return (data == b.data) || (data && b.data && *data == *b.data);
    }

    bool operator<(Alphabet const &b) const noexcept {
        return ((!data || !b.data) && data < b.data) || (data && b.data && *data < *b.data);
    }

    bool is_palindromic(Subdomain) const noexcept;
    bool is_palindromic(Subsequence) const noexcept;

    json save_repr() const;
    void load_repr(json const &);
};

/******************************************************************************************/

extern Alphabet const DNA;
extern Alphabet const RNA;
extern Alphabet const RNADNA;
extern Alphabet const RNA2OMeRNA;

void disable_noncomplement_closing(BasePairing &, Alphabet const &);

/******************************************************************************************/

struct AlphabetConverter {
    std::array<Base, Base::capacity> map;

    NUPACK_REFLECT(AlphabetConverter, map);

    Base operator()(Base b) const noexcept {return map[+b];}

    AlphabetConverter() = default;
    AlphabetConverter(Alphabet const &from, Alphabet const &to);
};

/******************************************************************************************/

Base standardize_base(Base b, Alphabet const &a_orig, Alphabet const &a_new);

template <class S>
auto convert_sequence(S const &seq, Alphabet const &a_orig, Alphabet const &a_new){
    return indirect_view(seq, [&a_orig, &a_new](Base b) {return standardize_base(b, a_orig, a_new);});
}

/******************************************************************************************/

}

namespace std {

template <>
struct hash<nupack::AlphabetData> {
    size_t operator()(nupack::AlphabetData const &) const noexcept;
};


template <>
struct hash<nupack::Alphabet> {
    size_t operator()(nupack::Alphabet const &a) const noexcept {
        return a.data ? hash<nupack::AlphabetData>()(*a.data) : 0;
    }
};

}