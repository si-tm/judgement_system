#pragma once
#include "Bind.h"
#include <nupack/types/Sequence.h>
#include <nupack/types/Named.h>
#include <nupack/math/Sparse.h>
#include <nupack/execution/Executor.h>

namespace nupack {

/******************************************************************************/

void render(Document &doc, Type<Sparsity>);

template <class T>
void render(Document &doc, Type<PairMatrix<T>> t) {
    doc.type(t, "core.PairMatrix");
    doc.method(t, "new", rebind::construct<Mat<T> const &, Sparsity>(t));
    doc.method(t, "to_array", &PairMatrix<T>::full);
    doc.method(t, "defect", [](PairMatrix<T> const &m, PairList const &p) {return m.defect(p);});
    render_public(doc, t);
    render_json(doc, t);
}

/******************************************************************************/

std::optional<PairList> request(Type<PairList>, rebind::Variable const &r, rebind::Dispatch &msg);
std::optional<Sequence> request(Type<Sequence>, rebind::Variable const &r, rebind::Dispatch &msg);
std::optional<Domain> request(Type<Domain>, rebind::Variable const &r, rebind::Dispatch &msg);

void render(Document &doc, Type<SharedExecutor> t);

void render(Document &doc, Type<AlwaysFalse>);
void render(Document &doc, Type<AlwaysTrue>);

void render(Document &doc, Type<Alphabet> t);

void render(Document &doc, Type<Base> t);
void render(Document &doc, Type<Wildcard> t);

void render(Document &doc, Type<Sequence> t);
void render(Document &doc, Type<Domain> t);

void render(Document &doc, Type<NamedSequence> t);
void render(Document &doc, Type<NamedDomain> t);
void render(Document &doc, Type<NamedDomainList> t);

void render(Document &doc, Type<Complex> t);
void render(Document &doc, Type<NamedComplex> t);
void render(Document &doc, Type<ComplexSet> t);
void render(Document &doc, Type<Tube> t);

void render(Document &doc, Type<TargetStrand> t);
void render(Document &doc, Type<TargetComplex> t);
void render(Document &doc, Type<TargetTube> t);

void render(Document &doc, Type<Local> t);
void render(Document &doc, Type<PairList>);
void render(Document &doc, Type<BasePairing> t);
void render(Document &doc, Type<Structure> t);
void render(Document &doc, Type<DomainList> t);

/******************************************************************************/

}

namespace rebind {

// Have to put these as explicit specialization because the aliased types
// aren't actually nupack:: classes (they're std::vector)

template <>
struct Renderer<nupack::NamedDomainList> {
    void operator()(Document &doc) const {nupack::render(doc, Type<nupack::NamedDomainList>());}
};

template <>
struct ImplicitConversions<nupack::NamedComplex> {
    using Base = nupack::Complex;
    using types = decltype(concat(Pack<Base>(), typename ImplicitConversions<Base>::types()));
};

template <>
struct Renderer<nupack::DomainList> {
    void operator()(Document &doc) const {nupack::render(doc, Type<nupack::DomainList>());}
};

template <>
struct Renderer<nupack::Complex> {
    void operator()(Document &doc) const {nupack::render(doc, Type<nupack::Complex>());}
};

// template <>
// struct Request<nupack::Complex> {
//     std::optional<nupack::Complex> operator()(Variable const &r, Dispatch &msg) const;
// };

}
