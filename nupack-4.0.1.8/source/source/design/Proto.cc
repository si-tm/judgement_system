#include "nupack/algorithms/Macro.h"
#include "nupack/design/DesignParameters.h"
#include "nupack/design/Objectives.h"
#include "nupack/proto/Spec.h"
#include <nupack/proto/Design.h>
#include <nupack/proto/Convert.h>

namespace nupack::jobs {

/******************************************************************************/

design::DesignParameters from_spec(spec::DesignOptions const &in) {
    return {
        .seed = static_cast<std::size_t>(in.seed),
        .f_stop = in.f_stop,
        .f_passive = in.f_passive,
        .H_split = in.h_split,
        .N_split = int(in.n_split),
        .f_split = in.f_split,
        .f_stringent = in.f_stringent,
        .dG_clamp = in.dg_clamp,
        .M_bad = int(in.m_bad),
        .M_reseed = int(in.m_reseed),
        .M_reopt = int(in.m_reopt),
        .f_redecomp = in.f_redecomp,
        .f_refocus = in.f_refocus,
        .f_sparse = in.f_sparse,
        .max_time = in.max_time/real(1000),
        .wobble_mutations = in.wobble_mutations
    };
}

spec::DesignOptions to_spec(design::DesignParameters const &in) {
    return {
        .seed = static_cast<std::size_t>(in.seed),
        .f_stop = in.f_stop,
        .f_passive = in.f_passive,
        .h_split = in.H_split,
        .n_split = static_cast<std::uint32_t>(in.N_split),
        .f_split = in.f_split,
        .f_stringent = in.f_stringent,
        .dg_clamp = in.dG_clamp,
        .m_bad = static_cast<std::uint32_t>(in.M_bad),
        .m_reseed = static_cast<std::uint32_t>(in.M_reseed),
        .m_reopt = static_cast<std::uint32_t>(in.M_reopt),
        .f_redecomp = in.f_redecomp,
        .f_refocus = in.f_refocus,
        .f_sparse = in.f_sparse,
        .max_time = static_cast<std::uint64_t>(in.max_time*1000),
        .wobble_mutations = in.wobble_mutations
    };
}

/******************************************************************************/

design::Weight from_spec(spec::DesignWeight const &w) {
    return {.tube=w.tube, .complex=w.complex, .strand=w.strand, .domain=w.domain, .weight=w.weight};
};

spec::DesignWeight to_spec(design::Weight const &w) {
    return {.tube=w.tube, .complex=w.complex, .strand=w.strand, .domain=w.domain, .weight=w.weight};
};

/******************************************************************************/

// # defaults to 2 for RNA, 3 for DNA
uint get_H_split(spec::Parameters p) {
    switch (p) {
        case spec::Parameters::RNA: {return 2;}
        case spec::Parameters::DNA: {return 3;}
        case spec::Parameters::RNA95: {return 2;}
        case spec::Parameters::RNA95_NUPACK3: {return 2;}
        case spec::Parameters::RNA99_NUPACK3: {return 2;}
        case spec::Parameters::DNA04_NUPACK3: {return 3;}
    }
    NUPACK_ERROR("Bad parameters value", static_cast<int>(p));
}

/******************************************************************************/

NUPACK_UNARY_FUNCTOR(get_name, t.name);

struct Lookup {
    std::unordered_map<string, NamedDomain> domains;
    std::unordered_map<string, TargetStrand> strands;
    std::unordered_map<string, TargetComplex> complexes;

    NUPACK_REFLECT(Lookup, domains, strands, complexes);

    auto domain() const {return [&](string const &name) {
        if (auto k = domains.find(name); k != domains.end()) return k->second;
        NUPACK_ERROR("Constraint refers to unknown domain", name);
    };};
};

/******************************************************************************/

design::MatchConstraint from_spec(spec::MatchConstraint const &c, Lookup const &lu) {
    return {vmap(c.left, lu.domain()), vmap(c.right, lu.domain())};
}

spec::MatchConstraint to_spec(design::MatchConstraint const &c) {
    return {vmap(c.left, get_name), vmap(c.right, get_name)};
}

/******************************************************************************/

design::ComplementarityConstraint from_spec(spec::ComplementarityConstraint const &c, Lookup const &lu) {
    return {vmap(c.left, lu.domain()), vmap(c.right, lu.domain()), c.wobble_mutations};
}

spec::ComplementarityConstraint to_spec(design::ComplementarityConstraint const &c) {
    return {vmap(c.left, get_name), vmap(c.right, get_name), c.wobble_mutations};
}

/******************************************************************************/

design::DiversityConstraint from_spec(spec::DiversityConstraint const &c, Lookup const &lu) {
    return {.domains=vmap(c.domains, lu.domain()),
        .word_length=c.word_length, .min_nucleotide_types=c.min_nucleotide_types};
}

spec::DiversityConstraint to_spec(design::DiversityConstraint const &c) {
    return {.domains=vmap(c.domains, get_name),
        .word_length=c.word_length, .min_nucleotide_types=c.min_nucleotide_types};
}

/******************************************************************************/

design::WindowConstraint from_spec(spec::WindowConstraint const &c, Lookup const &lu) {
    return {vmap(c.domains, lu.domain()),
            vmap(c.sequences, [](auto const &s) {return DNA.domain(s);})};
}

spec::WindowConstraint to_spec(design::WindowConstraint const &c) {
    return {.domains=vmap(c.domains, get_name),
            .sequences=vmap(c.comparisons.at(0), [](Domain const &d) {return DNA.to_string(d);})};
}

/******************************************************************************/

design::LibraryConstraint from_spec(spec::LibraryConstraint const &c, Lookup const &lu) {
    return {vmap(c.domains, lu.domain()), vmap(c.libraries, [](auto const &v) {
        return vmap(v.sequences, [](auto const &s) {return DNA.domain(s);});
    })};
}

spec::LibraryConstraint to_spec(design::LibraryConstraint const &c) {
    return {.domains=vmap(c.domains, get_name),
            .libraries=vmap(c.comparisons, [](auto const &v) -> spec::SequenceList {
                return {vmap(v, [](Domain const &d) {return DNA.to_string(d);})};
            })};
}

/******************************************************************************/

design::PatternConstraint from_spec(spec::PatternConstraint const &c, Lookup const &lu, real weight=1) {
    return {c.domains.empty() ? item_view(lu.domains) : vmap(c.domains, lu.domain()),
            vmap(c.patterns, [](string const &s) {return DNA.domain(s);}), weight};
}

spec::PatternConstraint to_spec(design::PatternConstraint const &c) {
    return {.domains=vmap(c.domains, get_name),
            .patterns=vmap(c.patterns, [](Domain const &d) {return DNA.to_string(d);})};
}

/******************************************************************************/

design::SimilarityConstraint from_spec(spec::SimilarityConstraint const &c, Lookup const &lu, real weight=1) {
    return {.domains = vmap(c.domains, lu.domain()), .reference=DNA.domain(c.reference), .range={c.min, c.max}, .weight=weight};
}

spec::SimilarityConstraint to_spec(design::SimilarityConstraint const &c) {
    return {.domains=vmap(c.domains, get_name),
        .reference=DNA.to_string(c.reference),
        .min=c.range.first, .max=c.range.second};
}

/******************************************************************************/

design::SSMObjective from_spec(spec::SsmConstraint const &c, Lookup const &lu, real weight) {
    vec<TargetComplex> complexes;
    if (c.complexes.empty()) {
        for (auto const &[k, v] : lu.complexes)
            if (!v.structure.empty()) complexes.emplace_back(v);
    } else {
        for (auto const &k : c.complexes) complexes.emplace_back(lu.complexes.at(k));
    }

    return {std::move(complexes), c.word_size, weight};
}

spec::SsmConstraint to_spec(design::SSMObjective const &c) {
    return {vmap(c.complexes, [](auto const &c) {return c.name;}), c.word_size};
}

/******************************************************************************/

design::EnergyEqualizationObjective from_spec(spec::EnergyMatchConstraint const &c, Lookup const &lu, real weight) {
    return {vmap(c.domains, [&](auto const &d) {return lu.domains.at(d);}),
        std::isfinite(c.energy_ref) ? std::make_optional(c.energy_ref) : std::nullopt, weight};
}

spec::EnergyMatchConstraint to_spec(design::EnergyEqualizationObjective const &c) {
    return {vmap(c.domains, [](auto const &d) {return d.name;}), c.ref_energy ? *c.ref_energy : std::numeric_limits<real>::infinity()};
}

/******************************************************************************/

spec::HardConstraint to_spec(design::HardConstraint const &c) {
    return std::visit([](auto const &c) -> spec::HardConstraint {return to_spec(c);}, c);
}

/******************************************************************************/

spec::SoftConstraint to_spec(design::SoftConstraint const &c) {
    return std::visit([](auto const &c) -> spec::SoftConstraint {return {to_spec(c), c.weight};}, c);
}

/******************************************************************************/

design::Specification from_spec(spec::DesignJob const &in) {
    design::Specification s;
    s.objective_weight = in.options.objective_weight;
    s.model = from_spec(in.model);
    s.options = from_spec(in.options);

    if (s.options.H_split < 0)
        s.options.H_split = get_H_split(in.model.parameters);

    Lookup lu;

    for (auto const &d : in.domains) lu.domains.try_emplace(d.name, DNA.domain(d.sequence), d.name, DNA);

    for (auto const &s : in.strands) lu.strands.try_emplace(s.name, TargetStrand{
        .domains=vmap(s.domains, [&](auto const &d) {return lu.domains.at(d);}), .name=s.name});

    for (auto const &x : in.complexes) lu.complexes.try_emplace(x.name,
        vmap(x.strands, [&](auto const &s) {return lu.strands.at(s);}), Structure(x.structure), x.name, 0.0);

    s.hard_constraints = vmap(in.hard_constraints, [&](auto const &c) {
        return std::visit([&](auto const &c) -> design::HardConstraint {return from_spec(c, lu);}, c);
    });

    s.soft_constraints = vmap(in.soft_constraints, [&](auto const &c) {
        return std::visit([&, w=c.weight](auto const &c) -> design::SoftConstraint {return from_spec(c, lu, w);}, c.constraint);
    });

    s.defect_weights = vmap(in.weights, [](auto const &w) {return from_spec(w);});

    s.tubes = vmap(in.tubes, [&](auto const &t) {
        return TargetTube(vmap(t.complexes, [&](auto const &x) {return lu.complexes.at(x);}), t.concentrations, t.name);
    });

    return s;
}

/******************************************************************************/

spec::DesignJob to_spec(design::Specification const &d) {
    spec::DesignJob job{
        .model=to_spec(d.model),
        .options=to_spec(d.options),
        .hard_constraints=vmap(d.hard_constraints, [](auto const &c) {return to_spec(c);}),
        .soft_constraints=vmap(d.soft_constraints, [](auto const &c) {return to_spec(c);}),
        .weights=vmap(d.defect_weights, [](auto const &c) {return to_spec(c);})
    };

    std::set<string> complexes, strands, domains;
    for (auto const &t : d.tubes) {
        job.tubes.emplace_back(spec::DesignTube{.name=t.name, .complexes=vmap(t.complexes, get_name),
                .concentrations={t.concentrations.begin(), find_if(t.concentrations, [](auto x) {return x == 0;})},
                .includes={},
                .excludes={},
                .max_size=0
            });
        for (auto const &x : t.complexes) {
            if (!complexes.emplace(x.name).second) continue;
            job.complexes.emplace_back(spec::DesignComplex{.name=x.name, .strands=vmap(x.strands, get_name), .structure=x.structure.dp()});

            for (auto const &s : x.strands) {
                if (!strands.emplace(s.name).second) continue;
                job.strands.emplace_back(spec::DesignStrand{.name=s.name, .domains=vmap(s.domains, get_name)});

                for (auto const &d : s.domains) {
                    if (!domains.emplace(d.name).second) continue;
                    job.domains.emplace_back(spec::DesignDomain{.name=d.name, .sequence=DNA.to_string(d)});
                }
            }
        }
    }
    return job;
}

/******************************************************************************/

spec::DesignStats to_spec(design::DesignStats const &d) {
    spec::DesignStats o;
    o.num_leaf_evaluations = d.num_leaf_evaluations;
    o.num_reseeds = d.num_reseeds;
    o.num_redecompositions = view(d.num_redecompositions);
    o.offtargets_added_per_refocus = view(d.offtargets_added_per_refocus);
    o.design_time = d.design_time;
    o.analysis_time = d.analysis_time;
    o.final_psi = {view(d.final_Psi.mask), d.final_Psi.deflate};
    o.seed = d.seed;
    return o;
}

design::DesignStats from_spec(spec::DesignStats const &d) {
    design::DesignStats o;
    o.num_leaf_evaluations = d.num_leaf_evaluations;
    o.num_reseeds = d.num_reseeds;
    o.num_redecompositions = view(d.num_redecompositions);
    o.offtargets_added_per_refocus = view(d.offtargets_added_per_refocus);
    o.design_time = d.design_time;
    o.analysis_time = d.analysis_time;
    o.final_Psi.mask = view(d.final_psi.mask);
    o.final_Psi.deflate = d.final_psi.deflate;
    o.seed = d.seed;
    return o;
}

/******************************************************************************/

spec::DesignResult to_spec(design::SingleResult const &r, design::DesignStats const &stats, spec::DefectHistory history, Model<> const &model, std::uint32_t trial) {
    NUPACK_ASSERT(!r.defects.empty(), "no design results have been generated");
    spec::DesignResult o;
    o.trial = trial;
    o.weighted_defects = r.weighted_defects;
    o.unweighted_defects = r.defects;

    o.unweighted_defects[0] = sum(r.tubes, [](auto const &t) {return t.normalized_defect;}) / len(r.tubes);

    o.domains = vmap(r.domains, [](auto const &d) {
        return spec::DesignDomain{.name=d.first, .sequence=DNA.to_string(d.second)};
    });

    o.complexes = vmap(r.complexes, [&](auto c) {
        return spec::DesignComplexResult{
            .name=c.name,
            .log_partition_function=c.log_partition_function,
            .defect=c.defect,
            .normalized_defect=c.normalized_defect,
            .pair_probabilities=to_spec(c.pair_probabilities),
            .energy=c.structure.empty() ? std::numeric_limits<real>::quiet_NaN() : model.structure_energy(c.sequence, c.structure)
        };
    });

    o.tubes = vmap(r.tubes, [](auto const &t) {
        return spec::DesignTubeResult{
            .name=t.name,
            .nucleotide_concentration=t.nucleotide_concentration,
            .defect=t.defect,
            .normalized_defect=t.normalized_defect,
            .complexes=vmap(t.complexes, [](auto const &c) -> spec::DesignTubeComplexResult {
                return {
                    .name=c.name,
                    .concentration=c.concentration,
                    .target_concentration=c.target_concentration,
                    .defect=c.defect,
                    .structural_defect=c.structural_defect,
                    .concentration_defect=c.concentration_defect,
                    .normalized_defect_contribution=c.normalized_defect_contribution
                };
            })
        };
    });

    o.stats = to_spec(stats);
    o.history = std::move(history);

    return o;
}

design::DesignResult from_spec(spec::DesignResult const &r) {
    design::DesignResult o;
    // o.model =
    // o.parameters =
    NUPACK_ERROR("not impl");
    o.stats = from_spec(r.stats);
    // o.objectives =
    // o.results =
    // o.weights =
    return o;
}

/******************************************************************************/

}