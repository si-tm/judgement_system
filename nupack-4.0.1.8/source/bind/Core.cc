#include "Bind.h"
#include <nupack/types/Sequence.h>
#include <nupack/types/PairList.h>
#include <nupack/types/Named.h>

namespace nupack {

/******************************************************************************************/

std::optional<PairList> request(Type<PairList>, rebind::Variable const &r, rebind::Dispatch &msg) {
    if (auto z = r.request<std::string_view>()) return PairList(*z);
    return msg.error("Cannot convert to PairList");
}

std::optional<Sequence> request(Type<Sequence>, rebind::Variable const &r, rebind::Dispatch &msg) {
    if (auto s = r.request<std::string_view>()) return DNA.sequence(*s);
    return msg.error("Cannot convert to Sequence");
}

std::optional<Domain> request(Type<Domain>, rebind::Variable const &r, rebind::Dispatch &msg) {
    if (auto s = r.request<std::string_view>()) return DNA.domain(*s);
    return msg.error("Cannot convert to Domain");
}


void render(Document &doc, Type<Sparsity> t) {
    doc.type(t, "core.Sparsity");
    doc.method(t, "new", rebind::construct<double, std::size_t>(t));
    doc.method(t, "new", rebind::construct<Sparsity>(t));
    render_public(doc, t);
}

void render(Document &doc, Type<SharedExecutor> t) {
    doc.type(t, "core.SharedExecutor");
    doc.method(t, "new", [](uint n) {return SharedExecutor(n);});
    doc.method(t, "threads", [](SharedExecutor const &e) {return e ? e.impl->num_workers() : 0;});
}

/******************************************************************************************/

void render(Document &doc, Type<PairList> t) {
    doc.type(t, "core.PairList");
    render_public(doc, t);
    static_assert(has_lt<PairList>);
    render_comparisons(doc, t);
    render_json(doc, t);
    render_hash(doc, t);
    doc.method(t, "new", rebind::construct<std::string_view>(t));
    doc.method(t, "new", rebind::construct<typename PairList::data_type>(t));
    doc.method(t, "^", [](PairList const &v, PairList const &w) {
        NUPACK_REQUIRE(len(v), ==, len(w));
        return v ^ w;
    });
    doc.method(t, "dp", &PairList::dp<pair_data_type const &>);
    doc.method(t, "pseudoknots", &PairList::pseudoknots);
}

/******************************************************************************************/

void render(Document &doc, Type<BasePairing> t) {
    doc.type(t, "core.BasePairing");
    render_comparisons(doc, t);
    render_json(doc, t);
    doc.method(t, "can_pair", [](BasePairing const &p, Base b, Base c) {
        NUPACK_REQUIRE(+b, <, p.length);
        NUPACK_REQUIRE(+c, <, p.length);
        return p.can_pair(b, c);
    });
    doc.method(t, "can_close", [](BasePairing const &p, Base b, Base c) {
        NUPACK_REQUIRE(+b, <, p.length);
        NUPACK_REQUIRE(+c, <, p.length);
        return p.can_close(b, c);
    });
    doc.method(t, "length", &BasePairing::length);
}

/******************************************************************************************/

void render(Document &doc, Type<Base> t) {
    doc.type(t, "core.Base");
    doc.method(t, "new", [](uint n) {NUPACK_REQUIRE(n, <=, Base::capacity); return Base::from_index(n);});
    doc.method(t, "new", rebind::construct<Base>(t));
    doc.method(t, "__pos__", [](Base b) {return uint(+b);});

    render_comparisons(doc, t);
    render_json(doc, t);
    render_hash(doc, t);
}

void render(Document &doc, Type<Wildcard> t) {
    doc.type(t, "core.Wilcard");
    doc.method(t, "new", [](WildcardIndex n) {return Wildcard::from_index(n);});
    doc.method(t, "new", [](Base b) {return Wildcard::from_index(WildcardIndex(1) << b.value);});
    doc.method(t, "new", rebind::construct<Wildcard>(t));
    doc.method(t, "__pos__", [](Wildcard b) {return uint(+b);});
    doc.method(t, "indices", &Wildcard::indices<vec<BaseIndex>>);
    doc.method(t, "first", &Wildcard::first);
    doc.method(t, "__bool__", &Wildcard::operator bool);
    doc.method(t, "is_determined", &Wildcard::is_determined);
    doc.method(t, "includes", [](Wildcard w, Base b) {return w.includes(b);});
    doc.method(t, "includes", [](Wildcard w, Wildcard b) {return w.includes(b);});

    doc.method(t, "~", &Wildcard::operator~);
    doc.method(t, "+", &Wildcard::operator+);
    doc.method(t, "^", [](Wildcard a, Wildcard b) {return a ^ b;});
    doc.method(t, "&", [](Wildcard a, Wildcard b) {return a & b;});
    doc.method(t, "|", [](Wildcard a, Wildcard b) {return a | b;});

    doc.function("core.are_compatible", [](Base a, Base b) {return are_compatible(a, b);});
    doc.function("core.are_compatible", [](Wildcard a, Base b) {return are_compatible(a, b);});
    doc.function("core.are_compatible", [](Base a, Wildcard b) {return are_compatible(a, b);});
    doc.function("core.are_compatible", [](Wildcard a, Wildcard b) {return are_compatible(a, b);});
    doc.function("core.are_compatible", [](Domain const &a, Domain const &b) {return are_compatible(a, b);});
    doc.function("core.are_compatible", [](Domain const &a, Sequence const &b) {return are_compatible(a, b);});
    doc.function("core.are_compatible", [](Sequence const &a, Domain const &b) {return are_compatible(a, b);});
    doc.function("core.are_compatible", [](Sequence const &a, Sequence const &b) {return are_compatible(a, b);});

    render_comparisons(doc, t);
    render_json(doc, t);
    render_hash(doc, t);
}

/******************************************************************************************/

void render(Document &doc, Type<Alphabet> t) {
    doc.type(t, "core.Alphabet");
    render_json(doc, t);
    render_hash(doc, t);
    render_comparisons(doc, t);

    //doc.method(t, "wildcard", &Alphabet::wildcard);
    //doc.method(t, "base", &Alphabet::base);
    doc.method(t, "complement", [](Alphabet const &a, Base b) {return a.complement(b);});
    doc.method(t, "complement", [](Alphabet const &a, Wildcard b) {return a.complement(b);});
    doc.method(t, "sequence", &Alphabet::sequence);
    doc.method(t, "domain", &Alphabet::domain);
    doc.method(t, "length", &Alphabet::length);
    doc.method(t, "letter", [](Alphabet const &a, Base b) {return a(b);});
    doc.method(t, "letter", [](Alphabet const &a, Wildcard b) {return a(b);});
    doc.method(t, "reverse_complement", &Alphabet::reverse_complement<Sequence>);
    doc.method(t, "reverse_complement", &Alphabet::reverse_complement<Domain>);

    doc.method(t, "new", [](std::string_view s, std::optional<std::string_view> mat) -> Alphabet {
        if (mat && !mat->empty()) {
            if (io::iequals(*mat, "rna")) return RNA;
            if (io::iequals(*mat, "dna")) return DNA;
            if (io::iequals(*mat, "rnadna")) return RNADNA;
            NUPACK_ERROR("Material is not supported", *mat);
        } else {
            bool t = s.find('T') != std::string_view::npos || s.find('t') != std::string_view::npos;
            bool u = s.find('U') != std::string_view::npos || s.find('u') != std::string_view::npos;
            bool r = s.find('r') != std::string_view::npos;
            bool d = s.find('d') != std::string_view::npos;
            if (!r && !d){
                if (t && !u) return DNA;
                if (u && !t) return RNA;
                if (u && t) return RNADNA;
                return DNA;
            } else{
                if (r && !d) return RNA;
                if (d && !r) return DNA;
                return RNADNA;
            }
        }
    });

    doc.function("core.standardize_alphabet", [](vec<NamedComplex> const &complexes, Alphabet const & new_alpha) {
        return vmap(complexes, [&](auto const &nc) {
            return vmap(nc.strands, [&](auto const &strand) {
                return NamedSequence(indirect_view(strand, [&](Base b) {
                    return standardize_base(b, strand.alphabet, new_alpha);
                }), strand.name, new_alpha);
            });
        });
    });
}

/******************************************************************************************/

template <class S>
void render_sequence_methods(Document &doc, Type<S> t) {
    render_comparisons(doc, t);
    render_json(doc, t);
    render_hash(doc, t);

    doc.method(t, "nt", [](S const &s) {return s.size();});
    doc.method(t, "[]", [](S const &s, std::ptrdiff_t i) {
        std::ptrdiff_t n = len(s);
        if (i < n && i >= -n) return s[i < 0 ? i + n : i];
        NUPACK_ERROR("invalid sequence index", i, len(s));
    });
    doc.method(t, "+", [](S s, S const &t) {s += t; return s;});
    doc.method(t, "__contains__", [] (S const &s, value_type_of<S> b) {return contains(s, b);});
    doc.method(t, "__len__", [](S const &s) {return s.size();});
    doc.method(t, "^", [](S const &x, S const &y) {
        NUPACK_REQUIRE(len(x), ==, len(y));
        return hamming_distance(x, y);
    });
}

/******************************************************************************************/

void render(Document &doc, Type<Sequence> t) {
    doc.type(t, "core.Sequence");
    render_sequence_methods(doc, t);
    doc.method(t, "new", [](string_view s, Alphabet const &a) {return a.sequence(s);});
    doc.method(t, "new", rebind::construct<Sequence>(t));
    doc.method(t, "{}", [](Sequence const &t) {return raw_sequence_string(t, 0);});

    pack<SequenceList, DomainList>::for_each([&](auto t) {
        pack<SequenceList, DomainList>::for_each([&](auto u) {
            doc.function("core.seq_distance", [](decltype(*t) const &a, decltype(*u) const &b) -> std::size_t {
                NUPACK_REQUIRE(len(a), ==, len(b), "Complexes are differently sized");
                return std::inner_product(a.begin(), a.end(), b.begin(), std::size_t(0), std::plus<>(), seq_distance);
            });
        });
    });
}

void render(Document &doc, Type<Domain> t) {
    doc.type(t, "core.RawDomain");
    render_sequence_methods(doc, t);
    doc.method(t, "new", [](string_view s, Alphabet const &a) {return a.domain(s);});
    doc.method(t, "new", rebind::construct<Domain>(t));
    doc.method(t, "new", [](Sequence const &s) {return Domain::from_sequence(s);});
    doc.method(t, "{}", [](Domain const &t) {return raw_domain_string(t, 0);});
}

/******************************************************************************************/

void render(Document &doc, Type<NamedSequence> t) {
    doc.type(t, "core.Strand");
    NUPACK_PUBLIC(t, name, alphabet);
    render_comparisons(doc, t);
    render_json(doc, t);
    render_hash(doc, t);

    doc.method(t, "new", rebind::construct<Sequence, string_view, Alphabet>(t));
    doc.method(t, "nt", [](NamedSequence const &s) {return len(s);});
    doc.method(t, "join", [](vec<NamedSequence> const &v, string_view name) {
        NUPACK_ASSERT(!v.empty());
        NamedSequence o({}, name, v[0].alphabet);
        for (auto const &s : v) o += s;
        return o;
    });

    doc.method(t, "~", std::bit_not<NamedSequence>());
    doc.method(t, "to_string", [](NamedSequence const &s, std::size_t n) {return s.alphabet.to_string(s, n);});
}

void render(Document &doc, Type<NamedDomain> t) {
    doc.type(t, "core.Domain");
    NUPACK_PUBLIC(t, name, alphabet);
    render_comparisons(doc, t);
    render_json(doc, t);
    render_hash(doc, t);

    doc.method(t, "new", rebind::construct<Domain, string_view, Alphabet>(t));
    doc.method(t, "nt", [](NamedDomain const &d) {return len(d);});
    doc.method(t, "join", [](vec<NamedDomain> const &v, string_view name) {
        NUPACK_ASSERT(!v.empty());
        NamedDomain o({}, name, v[0].alphabet);
        for (auto const &s : v) o += s;
        return o;
    });
    doc.method(t, "+", [](NamedDomain a, NamedDomain b) {return std::move(a) + std::move(b);});
    doc.method(t, "+", [](NamedDomain a, NamedDomainList b) {return std::move(a) + std::move(b);});

    doc.method(t, "~", std::bit_not<NamedDomain>());
    doc.method(t, "to_string", [](NamedDomain const &s, std::size_t n) {return s.alphabet.to_string(s, n);});
}

/******************************************************************************************/

void render(Document &doc, Type<NamedDomainList> t) {
    doc.type(t, "core.DomainList");
    render_comparisons(doc, t);
    render_json(doc, t);
    render_hash(doc, t);

    doc.method(t, "new", rebind::construct<NamedDomainList>(t));
    doc.method(t, "nt", [](NamedDomainList const &v) {return nt(v);});
    doc.method(t, "~", std::bit_not<NamedDomainList>());
    rebind::VectorRenderer<NamedDomainList>()(doc, false);

    doc.method(t, "+", [](NamedDomainList a, NamedDomain b) {return std::move(a) + std::move(b);});
    doc.method(t, "+", [](NamedDomainList a, NamedDomainList b) {return std::move(a) + std::move(b);});
}

/******************************************************************************************/

template <class F, class T>
void render_compare_key(Document &doc, Type<T> t) {
    doc.method(t, "<",  [](T const &a, T const &b) {return F::get(a) <  F::get(b);});
    doc.method(t, ">",  [](T const &a, T const &b) {return F::get(a) >  F::get(b);});
    doc.method(t, "<=", [](T const &a, T const &b) {return F::get(a) <= F::get(b);});
    doc.method(t, ">=", [](T const &a, T const &b) {return F::get(a) >= F::get(b);});
    doc.method(t, "==", [](T const &a, T const &b) {return F::get(a) == F::get(b);});
    doc.method(t, "!=", [](T const &a, T const &b) {return F::get(a) != F::get(b);});
    doc.method(t, "__hash__", [](T const &a) {return hash_of(F::get(a));});
}

struct NamedComplexCompare {
    static auto get(NamedComplex const &x) {return lowest_rotation(x.strands);}
};

void render(Document &doc, Type<NamedComplex> t) {
    doc.type(t, "core.Complex");
    doc.method(t, "new", rebind::construct<vec<NamedSequence>, string_view, real>(t));
    doc.method(t, "lowest_rotation", &NamedComplex::lowest_rotation);
    doc.method(t, "nt", &NamedComplex::nt);
    doc.method(t, "nstrands", [](NamedComplex const &c) {return len(c.strands);});
    doc.method(t, "nicks", [](NamedComplex const &c) {return complex_nicks(c.strands);});
    doc.method(t, "symmetry", [](NamedComplex const &c) {return rotational_symmetry(c.strands);});
    render_public(doc, t);
    render_json(doc, t);
    render_compare_key<NamedComplexCompare>(doc, t);
}

void render(Document &doc, Type<TargetStrand> t) {
    doc.render<DomainList>();

    doc.type(t, "core.TargetStrand");
    doc.method(t, "new", rebind::construct<vec<NamedDomain>, string>(t));
    doc.method(t, "nt", &TargetStrand::nt);
    render_public(doc, t);
    render_comparisons(doc, t);
    render_json(doc, t);
    render_hash(doc, t);
}

struct TargetComplexCompare {
    static auto get(TargetComplex const &x) {return lowest_rotation(x.strands);}
};

void render(Document &doc, Type<TargetComplex> t) {
    doc.type(t, "core.TargetComplex");
    doc.method(t, "new", rebind::construct<vec<TargetStrand>, Structure, string, real>(t));
    doc.method(t, "lowest_rotation", &TargetComplex::lowest_rotation);
    doc.method(t, "nt", &TargetComplex::nt);
    render_public(doc, t);
    render_json(doc, t);
    render_compare_key<TargetComplexCompare>(doc, t);
}

void render(Document &doc, Type<ComplexSet> t) {
    doc.type(t, "core.ComplexSet");
    doc.method(t, "new", rebind::construct<vec<NamedSequence>, vec<NamedComplex>>(t));
    doc.function("core.ComplexSet.union", [](vec<ComplexSet> const &v) {return ComplexSet::join(data_view(v));});
    doc.method(t, "+", [](ComplexSet const &t, ComplexSet const &u) {
        vec<ComplexSet> v{t, u};
        return ComplexSet::join(data_view(v));
    });
    // doc.method(t, "minimal_spec", &ComplexSet::minimal_spec);
    render_public(doc, t);
    render_comparisons(doc, t);
    render_json(doc, t);
    render_hash(doc, t);
};

void render(Document &doc, Type<Tube> t) {
    doc.render<ComplexSet>();
    doc.type(t, "core.Tube");
    doc.method(t, "new", rebind::construct<ComplexSet, vec<real>, string_view>(t));
    doc.function("core.Tube.union", [](vec<Tube> const &v) {return Tube::join(data_view(v));});
    doc.method(t, "+", [](Tube const &t, Tube const &u) {
        vec<Tube> v{t, u};
        return Tube::join(data_view(v));
    });
    render_public(doc, t);
    render_comparisons(doc, t);
    render_json(doc, t);
    render_hash(doc, t);
}

void render(Document &doc, Type<TargetTube> t) {
    doc.type(t, "core.TargetTube");
    doc.method(t, "new", rebind::construct<vec<TargetComplex>, vec<real>, string>(t));
    render_public(doc, t);
    render_comparisons(doc, t);
    render_json(doc, t);
    render_hash(doc, t);
}

template <class V>
void render_sequence_list(Document &doc, Type<V> t) {
    doc.method(t, "new", rebind::construct<V>(t));

    doc.method(t, "==", [](V const &a, V const &b) {return lowest_rotation(a) == lowest_rotation(b);});
    doc.method(t, "!=", [](V const &a, V const &b) {return lowest_rotation(a) != lowest_rotation(b);});
    rebind::VectorRenderer<V>()(doc, false);

    doc.method(t, "lowest_rotation", [](V const &v) {return lowest_rotation(v);});
    doc.method(t, "nt", [](V const &v) {return sum(v, len);});
    doc.method(t, "nicks", complex_nicks<V>);
    doc.method(t, "symmetry", rotational_symmetry<V>);
    render_json(doc, t);
    render_hash(doc, t);
}

void render(Document &doc, Type<SequenceList> t) {
    doc.type(t, "core.SequenceList");

    doc.method(t, "new", [](vec<std::string_view> const &v, Alphabet const &a) {return a.sequences(v);});
    doc.method(t, "new", [](string_view v, Alphabet const &a) {return a.sequences(v);});
    doc.method(t, "new", [](Sequence s) {return SequenceList{std::move(s)};});
    render_sequence_list(doc, t);
}

/******************************************************************************************/

void render(Document &doc, Type<DomainList> t) {
    doc.type(t, "core.RawDomainList");

    doc.method(t, "new", [](vec<std::string_view> const &v, Alphabet const &a) {return a.domains(v);});
    doc.method(t, "new", [](string_view v, Alphabet const &a) {return a.domains(v);});
    doc.method(t, "new", [](Domain s) {return DomainList{std::move(s)};});
    render_sequence_list(doc, t);
}

/******************************************************************************************/

std::uint32_t short_string_hash(std::string_view str) {
   std::uint32_t h = 0;
   for (auto p : str) h = 37 * h + p;
   return h;
}

/******************************************************************************************/

// void render(Document &doc, Type<Complex> t) {
//     doc.type(t, "core.RawComplex");

//     doc.method(t, "new", rebind::construct<Complex>(t));
//     doc.method(t, "[]", [](Complex &v, std::size_t i) {return v.at(i);});
//     doc.method(t, "__len__", [](Complex const &v) {return v.size();});
// //     doc.method(t, "==", compare_complexes<std::equal_to<>>);
// //     doc.method(t, "!=", compare_complexes<std::not_equal_to<>>);
// //     doc.method(t, "<", compare_complexes<std::less<>>);
// //     doc.method(t, ">", compare_complexes<std::greater<>>);
// //     doc.method(t, "<=", compare_complexes<std::less_equal<>>);
// //     doc.method(t, ">=", compare_complexes<std::greater_equal<>>);

// //     render_json(doc, t);
// //     render_hash(doc, t);

//     doc.method(t, "__contains__", [](Complex const &x, Strand const &s) {return contains(x, s);});
// //     doc.method(t, "__xor__", [](Complex const &a, Complex const &b) {
// //         return std::inner_product(a.begin(), a.end(), b.begin(), std::size_t(0),
// //             std::plus<std::size_t>(), [](Strand const &a, Strand const &b) {
// //                 NUPACK_REQUIRE(len(a), ==, len(b), "hamming_distance: unequal strand lengths");
// //                 return hamming_distance(a, b);
// //         });
// //     });
//     doc.method(t, "nt", [](Complex const &x) {return nt(x);});
//     doc.method(t, "symmetry", [](Complex const &x) {return rotational_symmetry(x);});
//     doc.method(t, "lowest_rotation", [](Complex x) {return lowest_rotation(x);});
//     doc.method(t, "nicks", [](Complex const &x) {return prefixes<Nicks>(false, indirect_view(x, len));});
// }

/******************************************************************************************/

// bool complex_eq(NamedComplex &t, NamedComplex &u) {
//     bool eq = t == u;
//     if (eq) {
//         if (t.name.empty() && !u.name.empty()) t.name = u.name;
//         if (u.name.empty() && !t.name.empty()) u.name = t.name;

//         if (t.bonus == 0 && u.bonus != 0) t.bonus = 0;
//         if (u.bonus == 0 && t.bonus != 0) u.bonus = 0;
//     }
//     return eq;
// }

// Complex const & get_strands(NamedComplex const &c) {return c;}
// auto get_strands(TargetComplex const &c) {return c.strands;}

// template <class C>
// std::string complex_name(C const &c) {
//     if (!c.name.empty()) return c.name;
//     string o;
//     o.push_back('(');
//     bool first = true;
//     for (auto const &s : get_strands(c)) {
//         if (!std::exchange(first, false)) o.push_back('+');
//         o += strand_name(s);
//     }
//     o.push_back(')');
//     return o;
// }

// /******************************************************************************************/

// void render(Document &doc, Type<NamedComplex> t) {
//     doc.type(t, "core.Complex");
//     doc.method(t, "new", rebind::construct<Complex, string, real>(t));
//     // doc.method(t, "new", rebind::construct<NamedComplex>(t));

//     doc.method(t, "name", complex_name<NamedComplex>);
// //     doc.method(t, "strands", &NamedComplex::strands);
//     doc.method(t, ".bonus", &NamedComplex::bonus);

//     // This is a nasty business below but is inescapable if dealing with colliding optional names...
//     doc.method(t, "==", [](NamedComplex &t, NamedComplex &u) {return complex_eq(t, u);});
//     doc.method(t, "!=", [](NamedComplex &t, NamedComplex &u) {return !complex_eq(t, u);});

//     doc.method(t, "<", std::less<NamedComplex>());
//     doc.method(t, ">", std::greater<NamedComplex>());

//     doc.method(t, "<=", std::less_equal<NamedComplex>());
//     doc.method(t, ">=", std::greater_equal<NamedComplex>());

//     render_json(doc, t);
//     render_hash(doc, t);
// }

// /******************************************************************************************/

// void render(Document &doc, Type<NamedStrand> t) {
//     doc.type(t, "named.Strand");
//     doc.method(t, "new", rebind::construct<Strand, Strand, string>(t));
//     doc.method(t, "new", rebind::construct<NamedStrand>(t));
//     doc.method(t, "~", &NamedStrand::operator~);
//     doc.method(t, "reverse_complement", &NamedStrand::reverse_complement);

//     doc.method(t, ".name", &NamedStrand::name);

//     render_comparisons(doc, t);
//     render_json(doc, t);
//     render_hash(doc, t);
// }

/******************************************************************************************/

// void render(Document &doc, Type<TargetStrand> t) {
//     doc.type(t, "core.TargetStrand");
//     doc.method(t, "new", rebind::construct<Strand, DomainList>(t));
//     doc.method(t, "new", rebind::construct<TargetStrand>(t));

//     doc.method(t, "reverse_complement", [](TargetStrand d, bool wobble) {
//         std::reverse(d.domains.begin(), d.domains.end());
//         for (auto &d : d.domains) d = Domain(strand_complement<Strand>(d, wobble));
// 	    static_cast<Strand&>(d) = strand_complement<Strand>(d, wobble);
//         return d;
//     });

//     doc.method(t, ".domains", &TargetStrand::domains);

//     render_comparisons(doc, t);
//     render_json(doc, t);
//     render_hash(doc, t);
// }

/******************************************************************************************/

// This is a nasty business below but is inescapable if dealing with colliding optional names and complexes D:
// bool target_complex_eq(TargetComplex &t, TargetComplex &u) {
//     bool eq = t == u;
//     if (eq) {
//         if (t.name.empty() && !u.name.empty()) t.name = u.name;
//         if (u.name.empty() && !t.name.empty()) u.name = t.name;

//         if (t.structure.empty() && !u.structure.empty()) t.structure = u.structure;
//         if (u.structure.empty() && !t.structure.empty()) u.structure = t.structure;

//         if (t.bonus == 0 && u.bonus != 0) t.bonus = 0;
//         if (u.bonus == 0 && t.bonus != 0) u.bonus = 0;
//     }
//     return eq;
// }

// void render(Document &doc, Type<TargetComplex> t) {
//     doc.type(t, "core.TargetComplex");
//     doc.method(t, "new", rebind::construct<vec<TargetStrand>, string, Structure, real>(t));
//     doc.method(t, "new", rebind::construct<TargetComplex>(t));
//     doc.method(t, "nt", &TargetComplex::nt);
//     doc.method(t, "__len__", [](TargetComplex const &d) {return len(d.strands);});
//     doc.method(t, "__getitem__", [](TargetComplex const &d, std::size_t i) {return d.strands.at(i);});

//     doc.method(t, "==", target_complex_eq);

//     doc.method(t, "!=", [](TargetComplex &t, TargetComplex &u) {return !target_complex_eq(t, u);});

//     doc.method(t, "<", std::less<TargetComplex>());
//     doc.method(t, ">", std::greater<TargetComplex>());

//     doc.method(t, "<=", std::less_equal<TargetComplex>());
//     doc.method(t, ">=", std::greater_equal<TargetComplex>());

//     doc.method(t, "name", complex_name<TargetComplex>);
//     doc.method(t, "is_named", [](TargetComplex const &x) {return !x.name.empty();});

//     NUPACK_PUBLIC(t, strands, structure, bonus);
//     render_json(doc, t);
//     render_hash(doc, t);
// }

/******************************************************************************************/

void render(Document &doc, Type<Structure> t) {
    doc.type(t, "core.Structure");
    doc.method(t, "new", rebind::construct(t));
    doc.method(t, "new", rebind::construct<string>(t));
    doc.method(t, "new", rebind::construct<Structure>(t));
    doc.method(t, "new", rebind::construct<PairList, Nicks>(t));
    doc.method(t, ".values", &PairList::values);

    doc.method(t, "dp",     &Structure::dp);
    doc.method(t, "nicks",             &Structure::nicks);

    render_comparisons(doc, t);
    // render_public(doc, t);
    render_json(doc, t);
    render_hash(doc, t);
}

/******************************************************************************************/

}
