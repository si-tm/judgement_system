#pragma once
#include "Sequence.h"
#include "Structure.h"
#include "Complex.h"
#include "Alphabet.h"
#include "../reflect/Hash.h"

namespace nupack {

/******************************************************************************************/

inline string complement_name(string s) {
    if (s.empty() || s.back() != '*') s.push_back('*');
    else s.pop_back();
    return s;
}

/******************************************************************************************/

struct NamedSequence : Sequence {
    using base_type = Sequence;
    using is_member_ordered = True;
    using load_repr = void;
    using save_repr = void;

    string name;
    Alphabet alphabet;

    NamedSequence() = default;
    NamedSequence(Sequence s, string_view n, Alphabet a) 
        : Sequence(std::move(s)), name(n), alphabet(std::move(a)) {}

    NUPACK_REFLECT_BASE(NamedSequence, Sequence, name, alphabet);

    NamedSequence operator~() const;

    auto compare_key() const {return std::tie(name, static_cast<Sequence const &>(*this), alphabet);}
};

/******************************************************************************************/

struct DomainMap {
    std::unordered_map<std::string, Domain> map;
    
    template <class V>
    void apply(V &v) const {for (auto &t : v) t.replace(*this);}
};

struct NamedDomain : Domain {
    using base_type = Domain;
    using is_member_ordered = True;
    using load_repr = void;
    using save_repr = void;

    string name;
    Alphabet alphabet;

    NamedDomain() = default;
    NamedDomain(Domain s, string_view n, Alphabet a) 
        : Domain(std::move(s)), name(n), alphabet(std::move(a)) {}

    NUPACK_REFLECT_BASE(NamedDomain, Domain, name, alphabet);

    NamedDomain operator~() const;
    void replace(DomainMap const &map) {
        if (auto it = map.map.find(name); it != map.map.end()) static_cast<Domain &>(*this) = it->second;
    }

    auto compare_key() const {return std::tie(name, static_cast<Domain const &>(*this), alphabet);}
};

/******************************************************************************************/

struct NamedComplex {
    using is_member_ordered = True;    

    vec<NamedSequence> strands;
    string name;
    real bonus = 0;

    NamedComplex() = default;
    NamedComplex(vec<NamedSequence>, string_view n, real b=0);

    NamedComplex lowest_rotation() const;
    auto nt() const noexcept {return sum(strands, len);}

    operator Complex() const {return vmap<Complex>(strands);} 

    NUPACK_REFLECT(NamedComplex, strands, name, bonus);
};

/******************************************************************************************/

struct ComplexSet {
    using is_member_ordered = True;

    vec<NamedSequence> strands;
    vec<NamedComplex> complexes;

    ComplexSet() = default;
    ComplexSet(vec<NamedSequence>, vec<NamedComplex>);

    static ComplexSet join(View<ComplexSet const *>);

    // SetSpec minimal_spec() const;

    NUPACK_REFLECT(ComplexSet, strands, complexes);
};

/******************************************************************************************/

struct Tube : ComplexSet {
    using is_member_ordered = True;
    using base_type = ComplexSet;

    vec<real> concentrations;
    string name;

    Tube() = default;
    Tube(ComplexSet c, vec<real> x, string_view n);

    static Tube join(View<Tube const *>);

    NUPACK_EXTEND_REFLECT(Tube, ComplexSet, concentrations, name);
};

/******************************************************************************************/

using NamedDomainList = vec<NamedDomain>;

NamedDomainList operator+(NamedDomainList, NamedDomainList);
NamedDomainList operator+(NamedDomainList, NamedDomain);
NamedDomainList operator+(NamedDomain, NamedDomainList);
NamedDomainList operator+(NamedDomain, NamedDomain);
NamedDomainList operator~(NamedDomainList const &v);
inline std::size_t nt(NamedDomainList const &v) {return sum(v, len);}

/******************************************************************************************/

// A Sequence with a list of component domains
struct TargetStrand {
    using is_member_ordered = True;

    vec<NamedDomain> domains;
    string name;

    NUPACK_REFLECT(TargetStrand, domains, name);

    void replace(DomainMap const &map) {map.apply(domains);}
    auto nt() const noexcept {return sum(domains, len);} 
};

/******************************************************************************************/

struct TargetComplex {
    using is_member_ordered = True;

    vec<TargetStrand> strands;
    Structure structure;
    string name;
    real bonus = 0;

    TargetComplex lowest_rotation() const;

    void resolve(TargetComplex const &);

    TargetComplex() = default;
    TargetComplex(vec<TargetStrand>, PairList, string_view, real);

    TargetComplex(vec<TargetStrand> v, Structure s, string_view n, real b)
        : TargetComplex(std::move(v), std::move(static_cast<PairList &>(s)), n, b) {
            if (!structure.empty()) NUPACK_REQUIRE(s.nicks, ==, structure.nicks, "Structure strand lengths do not match");
        }

    NUPACK_REFLECT(TargetComplex, strands, structure, name, bonus);
    
    void replace(DomainMap const &map) {map.apply(strands);}
    auto nt() const noexcept {return sum(strands, [](auto const &s) {return s.nt();});} 
};

/******************************************************************************************/

struct TargetTube {
    using is_member_ordered = True;

    vec<TargetComplex> complexes;
    vec<real> concentrations;
    string name;
    uint n_on_targets = 0;

    TargetTube() = default;
    TargetTube(vec<TargetComplex>, vec<real>, string_view);    

    auto on_targets() const {return view(complexes, 0, n_on_targets);}
    auto off_targets() const {return view(complexes, n_on_targets, complexes.size());}
    void replace(DomainMap const &map) {map.apply(complexes);}

    NUPACK_REFLECT(TargetTube, complexes, concentrations, name, n_on_targets);
};

/******************************************************************************************/

}

namespace std {
    template <> struct hash<nupack::NamedSequence>  {
        size_t operator()(nupack::NamedSequence const &) const noexcept;
    };

    template <> struct hash<nupack::NamedDomain>  {
        size_t operator()(nupack::NamedDomain const &) const noexcept;
    };

    template <> struct hash<nupack::NamedComplex> : nupack::MemberHash {};
    template <> struct hash<nupack::TargetComplex> : nupack::MemberHash {};
    template <> struct hash<nupack::TargetStrand> : nupack::MemberHash {};
    template <> struct hash<nupack::TargetTube> : nupack::MemberHash {};
    template <> struct hash<nupack::Tube> : nupack::MemberHash {};
    template <> struct hash<nupack::ComplexSet> : nupack::MemberHash {};
}

/******************************************************************************************/
