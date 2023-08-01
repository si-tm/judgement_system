#pragma once
#include "StackJoiner.h"
#include <nupack/model/Move.h>
#include <nupack/model/ModelVariants.h>
#include <nupack/iteration/Search.h>
#include <nupack/algorithms/Utility.h>

namespace nupack::kmc {

/******************************************************************************************/

struct NoLoopJoiner {
    Mat<real> const &join_propensities() const {NUPACK_ERROR("joining not enabled");}

    void update(Ignore, Ignore, Ignore, Ignore, Ignore) {}

    vec<JoinLoc> all_locs(Ignore, Ignore, Ignore, Ignore, Ignore) const {return {};}

    JoinLoc choose_loc(Ignore, Ignore, Ignore, Ignore, Ignore, Ignore, Ignore) const {NUPACK_ERROR("joining not enabled");}
};

/******************************************************************************************/

struct QuadraticLoopJoiner {
    Mat<real> const &join_propensities() const {NUPACK_ERROR("join propensities not enabled");}

    void update(Ignore, Ignore, Ignore, Ignore, Ignore) {}

    vec<JoinLoc> all_locs(SubsequenceList const &v, Base b1, Ignore, Ignore, Ignore) const {
        vec<JoinLoc> locs;
        for (auto s : iterators(v))  // iterate through the sequences
            for (auto b : iterators(*s).offset(1, -1)) if (*b == b1) locs.emplace_back(s-v.begin(), b);
        return locs;
    }

    JoinLoc choose_loc(Ignore, Ignore, Ignore, Ignore, Ignore, Ignore, Ignore) const {NUPACK_ERROR("join choice not enabled based on only one loop");}
};

/******************************************************************************************/

struct ProductLoopJoiner {
    using Active = True;

    NUPACK_REFLECT(ProductLoopJoiner, propensities);

    Mat<real> propensities;
    ProductLoopJoiner() = default;

    /**************************************************************************************/

    Mat<real> const &join_propensities() const {return propensities;};

    template <class State, class EM, class RF>
    void update(State const &, SubsequenceList const &, int, EM const &, RF const &);

    /// Iterate over all locations in the loop with their partial dG and partial rate constant
    template <class EM, class RF, class F>
    void for_locs(SubsequenceList const &, Base, Base, EM const &, RF const &, F &&) const;

    template <class State, class EM, class RF>
    JoinLoc choose_loc(SubsequenceList const &, double &, uint index, bool which, State const &, EM const &, RF const &rf) const;

    /// Return all locations that can form the specified base pair
    template <class EM, class RF>
    vec<JoinLoc> all_locs(SubsequenceList const &, Base, Base, EM const &, RF const &) const;
};

/******************************************************************************************/

struct StackAccum {
    real *data;
    friend void operator<<(StackAccum &a, real o) {*a.data++ += o;}

    void operator++() {++data;}

    template <class I>
    void operator+=(I i) {data += i;}
};

/******************************************************************************************/

// Remember that ProductLoopJoiner updates are replicated in BOTH for_locs and update
template <class State, class EM, class RF>
void ProductLoopJoiner::update(State const &state, SubsequenceList const &v, int n, EM const &em, RF const &rf) {
    auto const &joiner = state.joiner.product_joiner();
    propensities.zeros(joiner.length(em), 2);

    NUPACK_ASSERT(rf.has_bimolecular_propensity(), "ProductJoiner not compatible with both stacking and this rate function");

    if (has_subensemble(em.ensemble)) {
        auto const nick = find_nick(v);
        real * const l = propensities.memptr(), * const r = l + propensities.n_rows;
        StackAccumulator<real> acc(v, nick, em);
        for (auto s : indices(v)) for (auto i : indices(v[s]).offset(+1, -1)) {
            accumulate_left_stacks(StackAccum{l}, acc, em, joiner.base_pairs, v, nick, s, i);
            accumulate_right_stacks(StackAccum{r}, acc, em, joiner.base_pairs, v, nick, s, i);
        }
        propensities *= em.boltz(0.5 * em.join_penalty() - em.loop_energy(v, find_nick(v)));
    } else {
        em.dangle_switch([&](auto const &dangle) {
            for (auto s : iterators(v)) {
                auto const old5 = safe_dangle5(dangle, v, s); // at beginning
                auto const old3 = safe_dangle3(dangle, v, s); // at end

                for (auto b : iterators(*s).offset(1, -1)) {
                    izip(joiner.base_pairs, [&](auto index, auto const &bp) {
                        if (*b != bp.first) return;
                        Base const c = bp.second;
                        real new3 = (b > begin_of(*s) + 1) ? dangle.energy3(b[-1], b[0], c) : real(); // at beginning
                        real new5 = (b < end_of(*s) - 2) ? dangle.energy5(c, b[0], b[1]) : real();    // at end

                        // Beginning of sequence
                        real dG = 0.5 * (em.join_penalty() + em.terminal_penalty(*b, c));
                        if (b == begin_of(*s) + 1) dG += new3 - old5;
                        else dG += dangle.combine(new3, old5, 2 == b - begin_of(*s)) - old5;

                        // End of sequence
                        if (b == end_of(*s) - 2) dG += new5 - old3;
                        else dG += dangle.combine(new5, old3, 3 == end_of(*s) - b) - old3;

                        propensities(index, 0) += rf.bimolecular_propensity(-em.beta * dG);
                    });
                }
            }
        });
        // Fill in second column by looking for reversed base pair 
        izip(joiner.base_pairs, [&](auto i, auto const &p1) {
            izip(joiner.base_pairs, [&](auto j, auto const &p2) {
                if (p1.first == p2.second && p1.second == p2.first)
                    propensities(j, 1) = propensities(i, 0);
            });
        });
    }
}

/******************************************************************************************/

// ProductLoopJoiner::choose_loc could be more efficient to just calculate at that index, right now just programmed to calculate the whole vector
template <class State, class EM, class RF>
JoinLoc ProductLoopJoiner::choose_loc(SubsequenceList const &v, double &r, uint index, bool which, State const &w, EM const &em, RF const &rf) const {
    auto const &joiner = w.joiner.product_joiner();
    if (has_subensemble(em.ensemble)) {
        auto const nick = find_nick(v);

        Col<real> search(joiner.length(em), la::fill::zeros);
        search(index) -= r * em.boltz(em.loop_energy(v, find_nick(v)) - 0.5 * em.join_penalty()); // take off constant scale factors

        StackAccumulator<real> acc(v, nick, em);
        for (auto s : indices(v)) for (auto i : indices(v[s]).offset(+1, -1)) {
            real const last = search(index);

            if (!which) accumulate_left_stacks(StackAccum{search.memptr()}, acc, em, joiner.base_pairs, v, nick, s, i);
            else accumulate_right_stacks(StackAccum{search.memptr()}, acc, em, joiner.base_pairs, v, nick, s, i);
            
            if (search(index) >= 0) {
                real const factor = search(index) - last;
                real const dlogp = std::log(factor);
                r = -last / factor;
                return JoinLoc(s, v[s].begin() + i, dlogp, rf.bimolecular_propensity(dlogp));
            }
        }
        NUPACK_ERROR("Join choice failure (stacking)", v, r, index, which, w);
    } else {
        JoinLoc ret;
        auto [x, y] = joiner.base_pairs[index];
        if (which) swap(x, y);
        for_locs(v, x, y, em, rf, [&](auto i, auto b, auto dlogp, auto hrate){
            if (!minus_divide_if(r, hrate)) return false;
            else ret = JoinLoc(i, b, dlogp, hrate); return true;
        });
        NUPACK_REQUIRE(ret.dlogp_propensity.value().second, >, 0, "Join choice failure (not stacking)", v, r, index, which, w);
        return ret;
    }
}

/******************************************************************************************/

// Iterate over all locations in the loop that can form base pair b1 to b2
// Remember that ProductLoopJoiner updates are replicated in BOTH for_locs and update
template <class EM, class RF, class F>
void ProductLoopJoiner::for_locs(SubsequenceList const &v, Base b1, Base b2, EM const &em, RF const &rf, F &&f) const {
    NUPACK_ASSERT(!has_subensemble(em.ensemble), "for_locs should not be used in stacking ensemble");
    
    em.dangle_switch([&](auto const &dangle) { // switch on the specified type of dangle
        for (auto s : iterators(v)) { // iterate through the sequences
            auto const old5 = safe_dangle5(dangle, v, s); // old dangle contribution at beginning
            auto const old3 = safe_dangle3(dangle, v, s); // old dangle contribution at end

            for (auto b : iterators(*s).offset(1, -1)) if (*b == b1) {
                real new3 = (b > begin_of(*s) + 1) ? dangle.energy3(b[-1], b[0], b2): real(); // at beginning
                real new5 = (b < end_of(*s) - 2) ? dangle.energy5(b2, b[0], b[1]) : real();   // at end

                // Split join penalty in 2 so it is properly counted
                real dG = 0.5 * (em.join_penalty() + em.terminal_penalty(b1, b2));
                // Beginning of sequence
                if (b == begin_of(*s) + 1) dG += new3 - old5;
                else dG += dangle.combine(new3, old5, 2 == b - begin_of(*s)) - old5;

                // End of sequence -- use dangle.combine to delegate none, min, all dangles
                if (b == end_of(*s) - 2) dG += new5 - old3;
                else dG += dangle.combine(new5, old3, 3 == end_of(*s) - b) - old3;
                NUPACK_ASSERT(contains_iter(*s, b));

                // Call the callback with the base location information, free energy, and partial rate constant
                // if the callback returns true, short-circuit and return early.
                real const dlogp = -em.beta * dG;
                if (f(s - begin_of(v), b, dlogp, rf.bimolecular_propensity(dlogp))) return;
            }
        }
    });
}

/******************************************************************************************/

template <class EM, class RF>
vec<JoinLoc> ProductLoopJoiner::all_locs(SubsequenceList const &v, Base b1, Base b2, EM const &em, RF const &rf) const {
    vec<JoinLoc> ret;
    /// Iterate over all locations in the loop with their partial dG and partial rate constant
    for_locs(v, b1, b2, em, rf, [&](auto i, auto b, auto dlogp, auto hrate){
        ret.emplace_back(i, b, dlogp, hrate);
        return false; // dont short circuit out of the iteration over locations
    });
    return ret;
}

/******************************************************************************************/

using LoopJoiner = Variant<NoLoopJoiner, ProductLoopJoiner, QuadraticLoopJoiner>;

}
