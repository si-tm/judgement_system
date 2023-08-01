#pragma once
#include "Designer.h"
#include "Weights.h"
#include "nupack/design/Objectives.h"

namespace nupack::design {

/** @brief specification convertible to CompConstraint or IdentConstraint */
struct CompareConstraint {
    vec<NamedDomain> left; // name of left domains
    vec<NamedDomain> right; // name of right domains

    std::pair<vec<int>, vec<int>> get_variables(DesignSequence const &) const;
    void replace(DomainMap const &map) {map.apply(left); map.apply(right);}

    template <class F>
    void register_definitions(F const &f, Ignore) const {
        for (auto const &d : left) f(d);
        for (auto const &d : right) f(d);
    }

    CompareConstraint() = default;

    CompareConstraint(vec<NamedDomain> l, vec<NamedDomain> r)
        : left(std::move(l)), right(std::move(r)) {}

    void check_constraint() const;
    NUPACK_REFLECT(CompareConstraint, left, right);
};


struct ComplementarityConstraint : CompareConstraint {
    using base_type = CompareConstraint;

    bool wobble_mutations = false;

    NUPACK_EXTEND_REFLECT(ComplementarityConstraint, base_type, wobble_mutations);

    ComplementarityConstraint() = default;

    ComplementarityConstraint(vec<NamedDomain> l, vec<NamedDomain> r, bool wobble=false)
        : CompareConstraint(std::move(l), std::move(r)), wobble_mutations(wobble) {}

    void implement_hard(DesignSequence &, Model<> const &) const;
};


struct MatchConstraint : CompareConstraint {
    using base_type = CompareConstraint;
    using base_type::base_type;

    void implement_hard(DesignSequence &, Ignore) const;
};

struct PairingConstraint : CompareConstraint {
    using base_type = CompareConstraint;
    using base_type::base_type;

    void implement_hard(DesignSequence &, Model<> const &) const;
};



/** @brief specification convertible to PatternConstraint */
struct PatternConstraint {
    vec<NamedDomain> domains; // names of domains
    DomainList patterns; // pattern to match
    real weight = 1;

    void implement_hard(DesignSequence &, Ignore) const;
    void check_constraint() const; 
    void replace(DomainMap const &map) {map.apply(domains);}

    template <class F>
    void register_definitions(F const &f, Ignore) const {for (auto const &d : domains) f(d);}

    NUPACK_REFLECT(PatternConstraint, domains, patterns, weight);
};


struct DiversityConstraint {
    vec<NamedDomain> domains;
    std::uint32_t word_length;
    std::uint32_t min_nucleotide_types;

    void implement_hard(DesignSequence &, Ignore) const;
    void check_constraint() const; 
    void replace(DomainMap const &map) {map.apply(domains);}

    template <class F>
    void register_definitions(F const &f, Ignore) const {for (auto const &d : domains) f(d);}

    NUPACK_REFLECT(DiversityConstraint, domains, word_length, min_nucleotide_types);
};


/** @brief specification convertible to LibraryConstraint */
struct LibraryConstraint {
    vec<NamedDomain> domains; // domain names
    vec<DomainList> comparisons;

    void replace(DomainMap const &map) {map.apply(domains);}
    void implement_hard(DesignSequence &, Ignore) const;
    void check_constraint() const; 

    template <class F>
    void register_definitions(F const &f, Ignore) const {for (auto const &d : domains) f(d);}

    NUPACK_REFLECT(LibraryConstraint, domains, comparisons);
};

// This is identical to Library but only length 1...
struct WindowConstraint : LibraryConstraint {
    WindowConstraint() = default;

    void check_constraint() const; 

    WindowConstraint(vec<NamedDomain> d, DomainList c)
        : LibraryConstraint{std::move(d), {std::move(c)}} {}
};


/** @brief specification convertible to MatchConstraint */
struct SimilarityConstraint {
    vec<NamedDomain> domains;
    Domain reference;
    std::pair<real, real> range;
    real weight = 1;

    void implement_hard(DesignSequence &, Ignore) const;
    void check_constraint() const; 

    template <class F>
    void register_definitions(F const &f, Ignore) const {for (auto const &d : domains) f(d);}
    void replace(DomainMap const &map) {map.apply(domains);}

    NUPACK_REFLECT(SimilarityConstraint, domains, reference, range, weight);
};

using HardConstraint = std::variant<
    ComplementarityConstraint,
    MatchConstraint,
    PairingConstraint,
    PatternConstraint,
    DiversityConstraint,
    WindowConstraint,
    SimilarityConstraint,
    LibraryConstraint
>;

using SoftConstraint = std::variant<
    SSMObjective,
    EnergyEqualizationObjective,
    PatternConstraint,
    SimilarityConstraint
>;

struct Specification {
    vec<TargetTube> tubes;
    Model<> model;
    DesignParameters options;
    vec<HardConstraint> hard_constraints;
    vec<SoftConstraint> soft_constraints;
    vec<Weight> defect_weights;
    real objective_weight = 1;

    Designer create_designer() const;

    vec<uint> ensure_compatibility(SingleResult const &res) const;

    void replace(DomainMap const &map);

    NUPACK_REFLECT(Specification, tubes, model, options, hard_constraints, soft_constraints, defect_weights, objective_weight);
};

vec<int> extract_variables(vec<NamedDomain> const &names, DesignSequence const &seqs);
vec<int> extract_element(NamedDomain const &name, DesignSequence const &seqs);

}
