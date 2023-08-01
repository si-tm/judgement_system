/** @brief contains classes for holding all the conceptual objects in a design
  (Design) and an object representing the design logic*/
#pragma once
#include "Design.h"
#include "DesignParameters.h"
#include "DesignComponents.h"
#include "Objectives.h"
#include "Result.h"
#include "Weights.h"

namespace nupack {
namespace design {

Logs get_logs(DesignParameters const &);
EngineObserver get_thermo_log(DesignParameters const &);

struct Designer {
    Design design;
    vec<Objective> objectives;

    DesignParameters parameters;
    Weights weights;
    EnsemblePartition Psi;

    int max_depth;

    DesignStats stats;
    Timer timer;
    Logs logs;
    EngineObserver obs;

    ResultState best {inf_result};
    // ArchiveState archive;

    std::set<Sequence> known_bads;

    std::function<bool(Designer &, bool)> checkpoint {AlwaysFalse()}; // checks if should checkpoint and does it. return whether design should finish.

    Designer() = default;
    Designer(Design d, vec<Objective> objs, Weights weights, DesignParameters params={}) :
            design(std::move(d)), objectives(std::move(objs)),
            parameters(std::move(params)),
            weights(std::move(weights)),
            Psi(design.complexes, parameters.f_passive * parameters.f_stop),
            logs(get_logs(parameters)),
            obs(get_thermo_log(parameters)) {}
            
    NUPACK_REFLECT(Designer, design, objectives, parameters, Psi, max_depth, weights, stats, best, known_bads)

    void initialize(Env const &, bool decompose=true);
    void subset_decompose(Env const &, vec<uint> subset, uint depth=0);

    void redecompose_active(Env const &env, uint depth);
    bool redecompose(Env const &env, uint depth, Sequence const &sequence);
    void refocus(Env const &env, Sequence const &sequence);


    /* alternate stuff */
    // Result alternate_optimize_tubes(Env const &env);
    // bool length_extrapolation_refocus(Env const &env);
    // bool sum_pf_refocus(Env const &env);
    /* end stuff */


    Result optimize_tubes(Env const &env);
    Result optimize_tubes_impl(Env const &env);
    Result optimize_forest(Env const &env, Sequence seq);
    Result optimize_leaves(Env const &env, Sequence seq);
    Result mutate_leaves(Env const &env, Sequence seq);

    /* multiobjective */
    Result evaluate_objectives(Env const &env, uint depth, EnsemblePartition const &part, Weights const &weights);
    Result reevaluate_objectives(Env const &env, Result const &res, uint depth, EnsemblePartition const &part, Weights const &weights);

    Sequence best_sequence(Env const &env);

    bool improvement_slowing(vec<uint> const &x, vec<real> const &y);

    auto time_elapsed() const { return stats.design_time + timer.elapsed(); }

    bool success() const {return best.full.weighted_total() <= parameters.f_stop;}
    void time_analysis(Env const &env);

    static constexpr auto repr_names() {return make_names("design", "parameters", "weights", "Psi", "stats", "timer", "best", "max_depth", "known_bads");}

    auto save_repr() const {return make_members(design, parameters, weights, Psi, stats,
            timer, best, max_depth, known_bads);}

    void load_repr(Design design, DesignParameters parameters, Weights weights, EnsemblePartition Psi, DesignStats stats,
            Timer timer, ResultState best, int max_depth, std::set<Sequence> known_bads) {
        *this = Designer();
        this->design = std::move(design);
        this->parameters = std::move(parameters);
        this->weights = std::move(weights);
        this->Psi = std::move(Psi);
        this->stats = std::move(stats);
        this->timer = std::move(timer);
        this->best = std::move(best);
        this->max_depth = std::move(max_depth);
        this->known_bads = std::move(known_bads);
    }
};




}}
