#include <nupack/design/ThermoWrapper.h>
#include <nupack/design/Models.h>
#include <nupack/design/DesignComponents.h>
#include <nupack/model/Model.h>
#include <nupack/thermo/Engine.h>

namespace nupack::design {

EngineObserver NullEngineObserver {};

template <class F>
auto maybe_repeat(EngineObserver &engobs, string_view c, uint n, bool b, F &&f) {
    double t = 0;
    if (!engobs.slowdown) return f();
    std::optional<decltype(f())> out;
    auto time = time_it(engobs.slowdown, [&] {out = f();});
    if (engobs.log) engobs.log(c, n, time, b);
    return std::move(*out);
}

/**
 * @brief An adapter for thermo::dynamic_program()
 * @details Used to avoid multiple recomputation of dynamic programming
 *     algorithm code.
 *
 * @param seqs The sequence for which to compute the partition function
 * @param models References to two models with the same conditions but with
 *     single- and double-precision datatypes
 *
 * @return The logarithm of the partition function
 */
// real partition_function(Env const &env, ::nupack::Complex const &seqs, Model<> const &model, EngineObserver &engobs) {
//     return maybe_repeat(engobs, "partition function", nt(seqs), false, [&]{
//         return thermo::compute(thermo::Job{seqs, thermo::Job::PF()}, model).pfunc->logq;
//         // return thermo::dynamic_program<3, 0, 0, 1, 1>(env, seqs, models);
//     });
// }


real log_partition_function(Env const &env, ::nupack::Complex const &seqs, ThermoEnviron &t_env, EngineObserver &engobs) {
    return maybe_repeat(engobs, "partition function", nt(seqs), true, [&]{
        std::unique_lock lock(t_env.mut.mut);
        return thermo::compute(thermo::Job{seqs, thermo::Job::PF()}, t_env.computers, t_env.compute_options(env)).pfunc.value().logq;
        // return thermo::dynamic_program<3, 0, 0, 1, 1>(env, seqs, t_env.doubled(), underlying_cache);
    });
}



// /**
//  * @brief An adapter for thermo::pair_probability()
//  * @details Used to avoid multiple recomputation of dynamic programming
//  *     algorithm code.
//  *
//  * @param seqs The sequence for which to compute the pair probabilities matrix
//  * @param models References to two models with the same conditions but with
//  *     single- and double-precision datatypes
//  *
//  * @return A pair with the pair probabilities matrix (Tensor) and the
//  *     logarithm of the partition function.
//  */
// std::pair<PairMatrix<real>, real> pair_probability(Env const &env, ::nupack::Complex const &seqs, Model<> const &model, EngineObserver &engobs) {
//     return maybe_repeat(engobs, "pair probability", nt(seqs), false, [&]{
//         auto r = thermo::compute(thermo::Job{seqs, thermo::Job::Pairs()}, model);
//         return std::make_pair(std::move(r.pairs->matrix), r.pfunc->logq);
//         // return thermo::pair_probability<3, 0, 0, 1, 1>(env, seqs, models);
//     });
// }

std::pair<PairMatrix<real>, real> raw_pair_probability(Env const &env, ::nupack::Complex const &seqs, ThermoEnviron &t_env, Sparsity s, EngineObserver &engobs) {
    try {
        return maybe_repeat(engobs, "pair probability", nt(seqs), true, [&]{
            std::unique_lock lock(t_env.mut.mut);
            auto r = thermo::compute(thermo::Job{seqs, thermo::Job::Pairs{s}}, t_env.computers, t_env.compute_options(env));
            NUPACK_ASSERT(r.pairs && r.pfunc, r.pairs, r.pfunc, seqs);
            return std::make_pair(std::move(r.pairs->matrix), r.pfunc->raw_logq);
            // return thermo::pair_probability<3, 0, 0, 1, 1>(env, seqs, t_env.doubled(), t_env.cache, obs);
        });
    } catch (Error const &e) {
        // auto it = t_env.cache.find(seqs);
        // if (it != end_of(t_env.cache)) {
        //     BEEP(it->first);
        //     BEEP(fork(it->second, [](auto const &x) {return std::get<1>(x.contents).size();}));
        // }

        BEEP(seqs);
        throw;
    }
}

/**
 * @brief Adapts thermo::pair_probability() for a sequence where fixed_pairs
 *     are forced to to pair by adding a bonus energy.
 * @details Essential function for computing conditional pair probabilities
 *     and partition functions for nodes in the decomposition tree. This
 *     function removes the extra bonus contributions from the partition
 *     function and pair probabilities before returning.
 *
 * @param seqs The sequence for which to compute the clamped pair
 *     probabilities matrix
 * @param models References to two CachedModels with the same conditions but
 *     with single- and double-precision datatypes
 * @param fixed_pairs pairs which are enforced by multiplying their \f$Q^b\f$
 *     element by bonus
 * @param bonus the energy whose boltzmann factor is multiplied into the
 *     \f$Q^b\f$ elements of the enforced pairs to enforce them
 * @return A pair with the pair probabilities matrix (Tensor) and the
 *     logarithm of the partition function. Both have bonuses removed already.
 */
std::pair<PairMatrix<real>, real> bonus_pair_probability(Env const &env, ::nupack::Complex seqs,
    ThermoEnviron &t_env, Sparsity sparsity, vec<SplitPoint> const &fixed, real bonus, EngineObserver &engobs) {
    izip(seqs, [](auto i, auto &s) {s.id = i;}); // Make strands identifiable based on id
    auto const &mod = t_env.model;
    std::set<SplitPoint> const fixed_pairs(fixed.begin(), fixed.end());     /* assumes unique pairs */
    auto const starts = prefixes(true, indirect_view(seqs, len));
    auto pairing = [&](auto const &v, auto i, auto j) -> std::pair<real, thermo::Action::kind> {
        i += starts[front(v).id]; // go to absolute indices
        j += starts[back(v).id];
        return {fixed_pairs.count(std::minmax(i, j)) ? bonus : 0, thermo::Action::bonus};
    };

    auto [P, logq] = maybe_repeat(engobs, "bonused pair probability", nt(seqs), false, [&]{
        auto options = t_env.compute_options(env);
        options.cache_complexes = false;
        options.action = thermo::Action{pairing};
        std::unique_lock lock(t_env.mut.mut);
        auto r = thermo::compute(thermo::Job{seqs, thermo::Job::Pairs{{.clamp=false}}}, t_env.computers, options); // Get all dense elements
        NUPACK_ASSERT(r.pairs && r.pfunc, r.pairs, r.pfunc, seqs);
        return std::make_pair(r.pairs->matrix.full(), r.pfunc->raw_logq);
    });

    real const exp_bonus = mod.boltz(bonus);

    // remove extra bonuses from fixed base pairs
    for (auto const &[i, j] : fixed_pairs) {
        P(i, j) = P(j, i) = P(i, j) / exp_bonus;
        NUPACK_QUICK_REQUIRE(P(i, j), >=, 0, i, j, bonus, exp_bonus, mod.beta, seqs, fixed_pairs);
        NUPACK_QUICK_REQUIRE(P(i, j), <=, 1.01, i, j, bonus, exp_bonus, mod.beta, seqs, fixed_pairs);
    }

    /* fix the diagonal unpaired probabilities */
    clamp_pair_matrix(P);

    // if (Debug) for (auto const &[i, j] : fixed_pairs) {NUPACK_REQUIRE(P(i, j), >=, 0.99, i, j, P);}

    if (mod.has_terminal_penalty) {
        Sequence s;
        for (auto const &seq : seqs) s.insert(s.end(), seq.begin(), seq.end());
        for (auto const &[i, j] : fixed_pairs) {
            logq += mod.terminal_penalty(s[i ? i : j], s[i ? j : i]) * mod.beta; // add back positive number usually
        }
    }

    // Before this code used QB but that was annoying to do so I took it out.
    logq -= len(fixed_pairs) * std::log(exp_bonus);
    // Fix salt bonus, I don't think the code was handling this properly before
    // Have to add half of the bonus for ever bonused base pair
    logq += -0.5 * mod.beta * len(fixed_pairs) * mod.parameters.info.loop_bias;

    if (std::isnan(logq)) NUPACK_ERROR("bonused DPA generated NaN", seqs, logq, len(fixed_pairs), std::log(exp_bonus));
    return std::make_pair(PairMatrix<real>(P, sparsity), logq);
}

}
