#pragma once

#include "ModelVariants.h"
#include "ParameterSet.h"
#include "StackEnumeration.h"
#include "../algorithms/Utility.h"
#include "../algorithms/Numeric.h"
#include "../standard/Ptr.h"
#include "../standard/Optional.h"
#include "../standard/Map.h"
#include "../types/Complex.h"
#include "../types/Alphabet.h"
#include "nupack/common/Constants.h"

namespace nupack {

/******************************************************************************************/

/// Condition descriptors, these are fine as reals since they aren't often used
struct ModelConditions {
    real temperature = DefaultTemperature;
    real na_molarity = 1.0;
    real mg_molarity = 0.0;

    real water_molarity() const {return ::nupack::water_molarity(temperature);}

    NUPACK_REFLECT(ModelConditions, temperature, na_molarity, mg_molarity);
    using is_member_ordered = True;
};

/******************************************************************************************/

template <class T>
struct Model : MemberOrdered {
    using value_type = T;

    /**************************************************************************************/

    NUPACK_REFLECT(Model, parameters, beta, conditions, has_terminal_penalty, ensemble);

    ParameterSet<T> parameters;
    ModelConditions conditions;
    T beta;
    Ensemble ensemble;
    bool has_terminal_penalty = false;

    /**************************************************************************************/

    Model() = default;

    explicit Model(Ensemble, ParameterFile const &p={}, ModelConditions const &cs={}, Optional<BasePairing> pairs={});

    explicit Model(Ensemble, ParameterSet<T> p, ModelConditions const &cs={}, Optional<BasePairing> pairs={});

    template <class U>
    explicit Model(Model<U> const &o) : parameters(o.parameters), conditions(o.conditions),
        beta(o.beta), ensemble(o.ensemble), has_terminal_penalty(o.has_terminal_penalty) {
        static_assert(!is_same<U, T>, "Should use normal copy constructor");
    }

    /**************************************************************************************/

    template <class ...Is>
    decltype(auto) dG(Is const &...is) const {return parameters(is...);}

    auto boltz(T e) const {return boltzmann_factor(beta, e);};

    /**************************************************************************************/

    template <class S> T hairpin_energy(const S &) const;
    template <class S1, class S2> T interior_energy(const S1 &, const S2 &) const;
    template <class V> T linear_multi_energy(V const &) const;
    template <class V> T multi_energy(V const &) const;
    template <class V> T multi_energy_nolinear(V const &) const;
    template <class V> T exterior_energy(V const &, int const) const;
    template <class V> T loop_energy(V const &, int) const;
    template <class V> T minimum_stack_energy(V const &, int) const;
    template <class V> T terminal_penalty_sum(V const &v) const;

    /**************************************************************************************/

    T interior_size_energy(int) const;
    T interior_asymmetry(int, int) const;
    T interior_mismatch(Base, Base, Base, Base) const;
    T interior_mismatch_1(Base, Base) const;
    T terminal_mismatch(Base, Base, Base, Base) const;

    T join_penalty() const {return dG(parameters.join_penalty());}
    T multi_init() const {return dG(parameters.multi_init());}
    T multi_base() const {return dG(parameters.multi_base());}
    T multi_pair() const {return dG(parameters.multi_pair());}
    auto pairs(Base i) const {return pairing().pairs(i);}
    T terminal_penalty(Base i, Base j) const {return dG(parameters.terminal_penalty(), +i, +j);}
    T coaxial_stack_energy(Base, Base, Base, Base) const;

    BasePairing const &pairing() const {return parameters.pairing;}
    Alphabet const &alphabet() const {return parameters.alphabet;}

    /**************************************************************************************/

    template <class F>
    auto dangle_switch(F &&f) const {
        switch (ensemble) {
            case Ensemble::nostacking: return fw<F>(f)(dangle_function(NoStacking(), parameters));
            case Ensemble::stacking: return fw<F>(f)(dangle_function(Stacking(), parameters));
            case Ensemble::dangle: return fw<F>(f)(dangle_function(Stacking(), parameters));
            case Ensemble::coaxial: return fw<F>(f)(dangle_function(Stacking(), parameters));
            case Ensemble::none: return fw<F>(f)(dangle_function(NoStacking(), parameters));
            case Ensemble::min: return fw<F>(f)(dangle_function(MinDangles(), parameters));
            case Ensemble::all: return fw<F>(f)(dangle_function(AllDangles(), parameters));
        }
    }

    auto ensemble_type() const {return ensemble_variant(ensemble);}

    T dangle5(Base i, Base j, Base k) const {return dG(parameters.dangle5(), +i, +j, +k);}
    T dangle3(Base i, Base j, Base k) const {return dG(parameters.dangle3(), +i, +j, +k);}

    T structure_energy(Complex const &, PairList p, bool distinguishable=true) const;
};

NUPACK_DEFINE_TEMPLATE(is_model, Model, class);

/******************************************************************************************/

template <class T>
Model<T>::Model(Ensemble e, ParameterSet<T> p, ModelConditions const &cs, Optional<BasePairing> pairs)
    : parameters(std::move(p)), conditions(cs),
      beta(1.0 / (Kb * cs.temperature)), ensemble(e) {

    // Use parameter wobble setting if not specified
    if (pairs) parameters.pairing = std::move(*pairs);

    if (e != Ensemble::nostacking && e != Ensemble::stacking) {
        if (e == Ensemble::none || e == Ensemble::min || e == Ensemble::all) {
            disable_noncomplement_closing(parameters.pairing, parameters.alphabet);
        } else {
            // Modify parameters. Safer to undo the shared shared_ptr in this case.
            parameters.modify([&](auto &ps) {
                if (e == Ensemble::dangle) {
                    fill(parameter_view(ps, ps.coaxial_stack), 1e100);
                } else if (e == Ensemble::coaxial) {
                    fill(parameter_view(ps, ps.dangle5), 1e100);
                    fill(parameter_view(ps, ps.dangle3), 1e100);
                    fill(parameter_view(ps, ps.terminal_mismatch), 1e100);
                }
            });
        }
    }

    for (auto i : alphabet().all()) for (auto j : alphabet().all())
        if (pairing().can_pair(i, j)) has_terminal_penalty |= terminal_penalty(i, j) != 0;
}

template <class T>
Model<T>::Model(Ensemble e, ParameterFile const &p, ModelConditions const &cs, Optional<BasePairing> pairs)
    : Model(e, ParameterSet<T>(load_parameter_set(ParameterInfo{p, "dG", dna_salt_correction(cs.temperature, cs.na_molarity, cs.mg_molarity), cs.temperature})), cs, std::move(pairs)) {}

/******************************************************************************************/

template <class T> template <class V>
T Model<T>::loop_energy(V const &v, int nick) const {
    if (nick != -1) return exterior_energy(v, nick);
    else if (len(v) == 1) return hairpin_energy(v[0]);
    else if (len(v) == 2) return interior_energy(v[0], v[1]);
    else return multi_energy(v);
}

/******************************************************************************************/

template <class T>
T Model<T>::interior_size_energy(int s) const {
    NUPACK_REQUIRE(s, >, 0, "interior size must be at least 1");
    if (s <= 30) return dG(parameters.interior_size(), s - 1); // Interior with >2, <30 total
    return static_cast<T>(dG(parameters.interior_size(), 29) + log((s) / 30.0) * dG(parameters.log_loop_penalty())); // Big interior loop
}

/******************************************************************************************/

template <class T>
T Model<T>::interior_asymmetry(int n1, int n2) const {
    auto ninio_number = std::min((decltype(n2)) 4, std::min(n2, n1)) - 1;
    auto asymmetry = std::abs(n1 - n2);
    return std::min<T>(asymmetry * dG(parameters.ninio(), ninio_number), dG(parameters.ninio(), 4));
}

/******************************************************************************************/

template <class T> /// Interior mismatch energy for (b1, b2, b3, b4) where b2 and b3 are paired, b1 left of b2, b4 right of b3
T Model<T>::interior_mismatch(Base b1, Base b2, Base b3, Base b4) const {
    return dG(parameters.interior_mismatch(), +b1, +b2, +b3, +b4);
}

/******************************************************************************************/

template <class T> /// Interior mismatch energy
T Model<T>::interior_mismatch_1(Base b1, Base b2) const {
    return dG(parameters.interior_mismatch_1(), +b1, +b2);
}

/******************************************************************************************/

template <class T> /// Terminal mismatch energy for (b1, b2, b3, b4) where b2 and b3 are paired, b1 left of b2, b4 right of b3
T Model<T>::terminal_mismatch(Base b1, Base b2, Base b3, Base b4) const {
    return dG(parameters.terminal_mismatch(), +b1, +b2, +b3, +b4);
}

/******************************************************************************************/

/// Coaxial stack energy for (b1, b2, b3, b4) where b1 and b2 are paired, b3 and b4 are paired,
template <class T> /// and b2 and b3 are on the same strand
T Model<T>::coaxial_stack_energy(Base b1, Base b2, Base b3, Base b4) const {
    return dG(parameters.coaxial_stack(), +b2, +b3, +b4, +b1);
}

/******************************************************************************************/

// Input 5' ---seq1--> 3'
//       3' <--seq2--- 5'
template <class T> template <class S1, class S2>
T Model<T>::interior_energy(S1 const &a, S2 const &b) const {
    T en = T();
    int const n1 = len(a) - 2, n2 = len(b) - 2;
    if (n1 == 0 && n2 == 0) { // Stack loop
        return en + dG(parameters.stack(), +a[0], +a[1], +b[0], +b[1]);
    } else if (n1 == 0 || n2 == 0) { // Bulge loop
        auto sz = std::max(len(a), len(b)) - 2;

        if (sz <= 30) en += dG(parameters.bulge_size(), sz - 1);
        else en += dG(parameters.bulge_size(), 29) + log(sz / 30.0) * dG(parameters.log_loop_penalty());

        if (sz == 1) {
            // add stacking term for single-base bulges. No terminal penalty here
            return en + dG(parameters.stack(), +front(a), +back(a), +front(b), +back(b)) - parameters.info.loop_bias;
        } else {
            // Terminal penalty applies otherwise
            if (has_terminal_penalty) en += terminal_penalty(front(a), back(b))
                                          + terminal_penalty(front(b), back(a));
            return en;
        }
    } else if (!parameters.is_condensed && n1 == 1 && n2 == 1) { // Interior 1x1
        return en + dG(parameters.interior_1_1(), +a[0], +a[1], +a[2], +b[0], +b[1], +b[2]);
    } else if (!parameters.is_condensed && n1 == 1 && n2 == 2) { // Interior 1x2
        return en + dG(parameters.interior_1_2(), +a[0], +a[1], +a[2], +b[0], +b[1], +b[2], +b[3]);
    } else if (!parameters.is_condensed && n1 == 2 && n2 == 1) { // Interior 2x1
        return en + dG(parameters.interior_1_2(), +b[0], +b[1], +b[2], +a[0], +a[1], +a[2], +a[3]);
    } else if (!parameters.is_condensed && n1 == 2 && n2 == 2) { // Interior 2x2
        return en + dG(parameters.interior_2_2(), +a[0], +a[1], +a[2], +a[3], +b[0], +b[1], +b[2], +b[3]);
    } else { // Big interior loop
        en += interior_size_energy(n1 + n2);
    }
    // interior loops n1 > 4 x n2 > 4 are size + asymmetry + mismatch
    en += interior_asymmetry(n1, n2); // after n1 > 4 && n2 > 4 this just depends on |n1-n2|

    if ((n1 == 1 && (parameters.is_condensed || n2 > 2))
     || (n2 == 1 && (parameters.is_condensed || n1 > 2))) {
        en += interior_mismatch_1(back(b), front(a));
        en += interior_mismatch_1(back(a), front(b));
    } else {
        en += interior_mismatch(back_index(b, 1), back(b), front(a), front(a, 1));
        en += interior_mismatch(back_index(a, 1), back(a), front(b), front(b, 1));
    }
    return en;
}

/******************************************************************************************/

template <class T> template <class S>
T Model<T>::hairpin_energy(S const &s) const {
    T en = T();
    if (len(s) < 5) return inf<T>();
    if (len(s) <= 32) en += dG(parameters.hairpin_size(), len(s) - 3);
    else en += dG(parameters.hairpin_size(), 29) + log((len(s) - 2) / 30.0) * dG(parameters.log_loop_penalty());

    if (len(s) == 5) { // triloop
        // terminal penalty is unintuitive but correct
        if (has_terminal_penalty) en += terminal_penalty(back(s), front(s));
        return en + dG(parameters.hairpin_tri(), +s[0], +s[1], +s[2], +s[3], +s[4]);
    } else if (len(s) == 6) { // tetraloop
        en += dG(parameters.hairpin_tetra(), +s[0], +s[1], +s[2], +s[3], +s[4], +s[5]);
    }
    return en + dG(parameters.hairpin_mismatch(), +back_index(s, 1), +back(s), +front(s), +front(s, 1));
}

/******************************************************************************************/

template <class T> template <class V>
T Model<T>::linear_multi_energy(V const &v) const {
    auto const n_unpaired = sum(v, len) - 2 * len(v);
    return len(v) * multi_pair() + multi_init() + multi_base() * n_unpaired;
}

/******************************************************************************************/

// Total of terminal penalty and impossible closing base pair energies for multi and exterior loop
template <class T> template <class V>
T Model<T>::terminal_penalty_sum(V const &v) const {
    if (ensemble == Ensemble::none || ensemble == Ensemble::min || ensemble == Ensemble::all) {
         if (front(front(v)) != Base::null() && front(front(v)) != alphabet().complement(back(back(v)))) return inf<T>();
         for (auto it : iterators(v).offset(0, -1))
             if (front(it[1]) != Base::null() && front(it[1]) != alphabet().complement(back(it[0]))) return inf<T>();
    }

    T t{};
    if (has_terminal_penalty) {
        for (auto it = begin_of(v); it != end_of(v) - 1; ++it)
            if (front(it[1]) != Base::null()) t += terminal_penalty(front(it[1]), back(it[0]));
        if (front(front(v)) != Base::null()) t += terminal_penalty(front(front(v)), back(back(v)));
    }
    return t;
}

/******************************************************************************************/

template <class T> template <class V>
T Model<T>::multi_energy(V const &v) const {
    return terminal_penalty_sum(v)
         + linear_multi_energy(v)
         + fork(ensemble_type(), [&](auto d) -> T {return stacking_energy(d, *this, v, -1);});
}

/******************************************************************************************/

template <class T> template <class V>
T Model<T>::multi_energy_nolinear(V const &v) const {
    return terminal_penalty_sum(v)
         + fork(ensemble_type(), [&](auto d) -> T {return stacking_energy(d, *this, v, -1);});
}

/******************************************************************************************/

// A list of sequences. nick indicates which one the strand break is before
// For example, the starting loop always has a nick of 0
template <class T> template <class V>
T Model<T>::exterior_energy(V const &v, int const nick) const {
    if (std::size(v) <= 1) return T();
    return terminal_penalty_sum(v)
         + fork(ensemble_type(), [&](auto d) -> T {return stacking_energy(d, *this, v, nick);});
}


/******************************************************************************************/

template <class T> template <class V>
T Model<T>::minimum_stack_energy(V const &v, int nick) const {
    if (!has_subensemble(ensemble) || (nick == Ether && len(v) < 3) || len(v) < 2) return loop_energy(v, nick);
    T e = terminal_penalty_sum(v) + minimum_stacking_energy<T>(v, nick, *this);
    if (nick == Ether) e += linear_multi_energy(v);
    return e;
}

/******************************************************************************************/

int find_loop_structure_nick(SequenceList const &, PairList const &);

/******************************************************************************************/

template <class Map, class Model>
Map loop_stacking_energies(Model const &m, Complex const &c, int nick=Ether) {
    Map out;
    auto const v = complex_to_loop(c, nick);
    if (nick != Ether || n_strands(c) > 2) {
        auto const t = m.terminal_penalty_sum(v) + (nick != Ether ? 0 : m.linear_multi_energy(v));
        enumerate_stacking_state_energies(v, nick, m, [&](auto const &p, auto e) {
            out.emplace(loop_stack_sequence_string(p), t + e);
        });
    } else {
        out.emplace(string(n_strands(c), 'n'), m.loop_energy(v, nick));
    }
    return out;
}

/******************************************************************************************/

}
