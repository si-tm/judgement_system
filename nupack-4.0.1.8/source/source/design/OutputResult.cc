#include <nupack/design/OutputResult.h>
#include <nupack/design/Designer.h>
#include <nupack/reflect/Serialize.h>
#include <limits>
#include <cmath>

namespace nupack::design {

// DesignResult::DesignResult(Specification const &spec, Designer const &designer, Result const &result) :
DesignResult::DesignResult(Env const &env, Designer const &designer) :
        model(at(designer.design.complexes, 0).target.model), parameters(designer.parameters),
        stats(designer.stats), objectives(designer.objectives),
        results({designer.best.full.full_evaluation(env, designer)}),
        weights(designer.weights),
        success(designer.success()) {}
        // results(indirect_view(designer.archive.full.results, [&](auto const &res) {
        //     return res.evaluated;
        // })) {}


SingleResult::SingleResult(Env const &env, Designer const &designer, Result const &res) {
    auto const &design = designer.design;
    auto const &seqs = design.sequences;
    auto const &sequence = res.sequence;
    /* domains */
    for (auto const &el : seqs.domains) { domains.emplace(el.first, el.second.to_sequence(sequence)); }
    /* strands */
    for (auto const &el : seqs.strands) { strands.emplace(el.first, el.second.to_sequence(sequence)); }

    auto const &models = design.models;

    auto & engobs = const_cast<EngineObserver &>(designer.obs);
    /* complexes */
    for (auto const &d : design.complexes) {
        ComplexResult comp;
        comp.name = d.name;
        comp.sequence = to_nick_sequence(d.strands, sequence);
        comp.structure = d.target.structure;
        auto x = d.log_pfunc(env, models, sequence, 0, {}, engobs);
        x = std::isfinite(x) ? x : std::numeric_limits<decltype(x)>::lowest();
        comp.log_partition_function = x;
        if (d.is_on_target()) {
            comp.pair_probabilities = d.pair_probabilities(env, models, sequence, 0, {}, engobs);
        }
        comp.defect = d.defect(env, models, sequence, 0, {}, engobs).total();
        comp.normalized_defect = comp.defect / nt(d);

        complexes.emplace_back(comp);
    }

    auto log_pfuncs = design.log_pfuncs(env, 0, {}, {}, engobs);
    auto complex_defects = design.complex_defects(env, 0, {}, {}, engobs);

    // real weighted = 0;
    /* tubes */
    izip(design.tubes, [&](auto i, auto const &t) {
        TubeResult tube;
        tube.name = t.name;
        tube.nucleotide_concentration = t.nucleotide_concentration;
        tube.defect = t.defect(log_pfuncs, complex_defects).total();
        tube.normalized_defect = tube.defect / tube.nucleotide_concentration;

        // weighted += t.normalized_defect(log_pfuncs, complex_defects, {}, designer.weights ? designer.weights.per_tube.at(i) : ComplexWeights()).total();
        auto concentrations = t.concentrations(log_pfuncs);

        zip(t.targets, concentrations, [&](auto const &c, auto conc) {
            auto const &ref = at(complexes, c.complex_index);
            TubeComplex comp;
            comp.name = ref.name;
            comp.concentration = conc;
            comp.target_concentration = c.target_conc;

            auto comp_defect = at(complex_defects, c.complex_index);
            comp.structural_defect = structural_defect(c, comp_defect, comp.concentration).total();
            comp.concentration_defect = concentration_defect(c, comp.concentration).total();
            comp.defect = comp.structural_defect + comp.concentration_defect;
            comp.normalized_defect_contribution = comp.defect / tube.nucleotide_concentration;

            tube.complexes.emplace_back(comp);
        });

        tubes.emplace_back(tube);
    });

    // Earlier version of code had diverging values of weighted vs unweighted objective when weight=1
    // Seemed like some mis-caching issue so I recalculated the defects here instead.
    NUPACK_REQUIRE(len(designer.objectives), ==, len(designer.weights.objective_weights));
    zip(designer.objectives, designer.weights.objective_weights, [&](auto &&o, auto const &w) {
        {
            auto opt = o.reevaluate(env, design, 0, {}, designer.weights, engobs);
            auto defect = opt ? std::move(*opt) : o.evaluate(env, design, 0, {}, designer.weights, engobs);
            weighted_defects.emplace_back(w * defect.total());
        }
        {
            auto opt = o.reevaluate(env, design, 0, {}, {}, engobs);
            auto defect = opt ? std::move(*opt) : o.evaluate(env, design, 0, {}, {}, engobs);
            defects.emplace_back(defect.total());
        }
    });
    // BEEP(defects, weighted_defects, weighted); Was testing that weighted_defect equaled weighted.
}





}
