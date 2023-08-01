#pragma once
#include "Design.h"
#include "../reflect/Serialize.h"

#include "../standard/Optional.h"
#include "../standard/Set.h"
// #include <boost/optional/optional_io.hpp>

#define NUPACK_BLANK(cls)
// #define NUPACK_BLANK(cls) auto members() {return ::nupack::make_members();} \
//     static constexpr auto names() {return ::nupack::make_names();} \
//     static constexpr auto accesses() {return ::nupack::make_members();}

namespace nupack::design {

struct Weights;
struct EngineObserver;

struct MultitubeObjective {
    MultitubeObjective() = default;
    void initialize(Design const &) {};

    Defect evaluate(Env const &env, Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const;
    Optional<Defect> reevaluate(Env const &env, Design const &design, uint depth,
            EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const;

    // NUPACK_BLANK(MultitubeObjective);
};


struct TubeObjective {
    string tube_name;
    uint tube_id;

    TubeObjective() = default;
    TubeObjective(string name) : tube_name(name) {}

    void initialize(Design const &);

    Defect evaluate(Env const &env, Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const;
    Optional<Defect> reevaluate(Env const &env, Design const &design, uint depth,
            EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const;

    NUPACK_REFLECT(TubeObjective, tube_name, tube_id);
};


struct ComplexObjective {
    string complex_name;
    uint complex_id;

    ComplexObjective() = default;
    ComplexObjective(string name) : complex_name(name) {}

    void initialize(Design const &);

    Defect evaluate(Env const &env, Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const;
    Optional<Defect> reevaluate(Env const &env, Design const &design, uint depth,
            EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const;

    NUPACK_REFLECT(ComplexObjective, complex_name, complex_id);
};

using NucleotideRelationMap = std::map<uint, std::set<uint>>;

struct SSMObjective {
    vec<TargetComplex> complexes;
    uint word_size;
    vec<uint> complex_ids;
    NucleotideRelationMap identicals;
    NucleotideRelationMap complements;
    std::set<vec<uint>> complement_restricted;
    /* possibly replace with lazy generator later */
    vec<vec<uint>> words;
    real weight = 1;
    uint normalization {0};

    SSMObjective() = default;

    SSMObjective(vec<TargetComplex> xs, uint word_size, real w=1.0) : complexes(std::move(xs)), word_size(word_size), weight(w) {}

    void replace(DomainMap const &map) {map.apply(complexes);}
    void check_constraint() const {}

    template <class F>
    void register_definitions(Ignore, F const &f) const {for (auto const &c : complexes) f(c);}

    void initialize(Design const &);
    void add_identicals(DesignSequence const &);
    void add_complements(DesignSequence const &);
    void process_words(Design const &);
    void process_structures(Design const &);
    vec<Range<uint>> ranges(Nicks const &) const;

    bool identical(vec<uint> const &, vec<uint> const &) const;
    bool complementary(vec<uint> const &, vec<uint> const &) const;

    Defect evaluate(Env const &env, Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const;
    Optional<Defect> reevaluate(Env const &env, Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const { return optional(); };

    NUPACK_REFLECT(SSMObjective, complexes, word_size, complex_ids,
            identicals, complements, complement_restricted, normalization, weight);
};


struct PatternObjective {
    vec<string> component_names;
    DomainList patterns;
    vec<Variant<DomainView, StrandView>> elements;

    std::map<uint, DomainList> grouped_patterns;
    uint normalization {0}; // number of potentially problematic windows

    void check_constraint() const {}

    PatternObjective() = default;
    PatternObjective(vec<string> names, DomainList prevents) : component_names(names), patterns(prevents) {}

    void initialize(Design const &);

    Defect evaluate(Env const &env, Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const;
    Optional<Defect> reevaluate(Env const &env, Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const { return optional(); };

    NUPACK_REFLECT(PatternObjective, component_names, patterns, elements, grouped_patterns, normalization);
};


struct SimilarityObjective {
    vec<NamedDomain> domains;
    Domain reference_sequence;
    std::pair<real, real> limits;
    vec<DomainView> elements;

    SimilarityObjective() = default;
    SimilarityObjective(vec<NamedDomain> names, Domain ref, std::pair<real, real> lims) :
            domains(std::move(names)), reference_sequence(std::move(ref)), limits(lims) {}

    void initialize(Design const &);

    Defect evaluate(Env const &env, Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const;
    Optional<Defect> reevaluate(Env const &env, Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const { return optional(); };

    NUPACK_REFLECT(SimilarityObjective, domains, reference_sequence, limits, elements);
};


struct EnergyEqualizationObjective {
    vec<NamedDomain> domains;
    Optional<real> ref_energy;
    vec<DomainView> domain_views;
    Model<real> model;
    real weight = 1;

    EnergyEqualizationObjective() = default;
    EnergyEqualizationObjective(vec<NamedDomain> domains, Optional<real> ref, real w=1) :
            domains(std::move(domains)), ref_energy(ref), weight(w) {}

    template <class F>
    void register_definitions(F const &f, Ignore) const {for (auto const &d : domains) f(d);}
    
    void replace(DomainMap const &map) {map.apply(domains);}
    void initialize(Design const &);

    void check_constraint() const {}

    Defect evaluate(Env const &env, Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const;
    Optional<Defect> reevaluate(Env const &env, Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const { return optional(); };

    NUPACK_REFLECT(EnergyEqualizationObjective, domains, ref_energy, domain_views, model, weight);
};


struct Objective {
    using Var = Variant<MultitubeObjective, TubeObjective, ComplexObjective, SSMObjective,
            SimilarityObjective, EnergyEqualizationObjective, PatternObjective>;
    Var variant;

    Objective() = default;
    Objective(Var const &var) : variant(var) {};

    void initialize(Design const &);
    Defect evaluate(Env const &env, Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const;
    Optional<Defect> reevaluate(Env const &env, Design const &design, uint depth,
            EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const;

    NUPACK_REFLECT(Objective, variant);

    ::nupack::json save_repr() const;
    // {
    //     std::map<uint, string> mapped {
    //             {0, "Multitube"},
    //             {1, "Tube"},
    //             {2, "Complex"},
    //             {3, "SSM"},
    //             {4, "Similarity"},
    //             {5, "EnergyEqualization"},
    //             {6, "Pattern"}
    //             };

    //     string type = mapped.at(variant.index());
    //     return fork(variant, [&](auto const &x) {
    //         ::nupack::json j;
    //         j["type"] = type;
    //         j["info"] = x;
    //         return j;
    //     });
    // }

    void load_repr(::nupack::json value);
    // {
    //     std::map<string, Var> mapped {
    //             {"Multitube", MultitubeObjective()},
    //             {"Tube", TubeObjective()},
    //             {"Complex", ComplexObjective()},
    //             {"SSM", SSMObjective()},
    //             {"Similarity", SimilarityObjective()},
    //             {"EnergyEqualization", EnergyEqualizationObjective()},
    //             {"Pattern", PatternObjective()}
    //             };

    //     auto type = value["type"].get<string>();
    //     fork(mapped.at(type), [&](auto const &t) {
    //         using obj_type = decay<decltype(t)>;
    //         variant = value["info"].get<obj_type>();
    //     });
    // }
};


template <class Container>
Optional<uint> find_multitube(Container const &objectives) {
    Optional<uint> ret;
    izip(objectives, [&](auto i, auto const &o) {
        if (maybe_get<MultitubeObjective>(o.variant) != nullptr) ret = i;
    });
    return ret;
}

/* how to handle resolving which sets of indices are the same */
struct Ind {
    vec<uint> indices;
    Optional<int> group {};

    Ind() = default;
    Ind(vec<uint> inds, Optional<int> grp={}) : indices(inds), group(grp) {}

    bool assigned() const {return bool(group);}

    NUPACK_REFLECT(Ind, indices, group);
};


struct Index_Map {
    std::vector<Ind> used;
    uint num_violations {0};

    Index_Map() = default;

    void add(vec<uint> indices) { used.emplace_back(indices, Optional<int>()); }
    /* return number of groups found */
    void resolve_groups(SSMObjective const &obj);
    real assign_blame(vec<real> &) const;

    NUPACK_REFLECT(Index_Map, used, num_violations);
};

}
