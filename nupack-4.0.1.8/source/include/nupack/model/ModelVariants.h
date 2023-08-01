#pragma once
#include <cmath>
#include "../types/Sequence.h"
#include "../standard/Variant.h"
#include "StackProgram.h"
#include "ParameterSet.h"

namespace nupack {

/******************************************************************************************/

/// simpler Rigs for use outside dynamic programs
struct RigPF {
    template <class ...Ts>
    static auto plus(Ts const &...ts) {return fold(::nupack::plus, ts...);}

    template <class ...Ts>
    static auto times(Ts const &...ts) {return fold(::nupack::times, ts...);}

    static real zero() {return 0;}
    static real one() {return 1;}

    static auto boltz(real beta, real t) {return boltzmann_factor(beta, t);}
};


struct RigMFE {
    template <class ...Ts>
    static auto plus(Ts const &...ts) {return fold(::nupack::min, ts...);}

    template <class ...Ts>
    static auto times(Ts const &...ts) {return fold(::nupack::plus, ts...);}

    static float zero() {return inf<float>();}
    static float one() {return 0;}

    static auto boltz(real, real t) {return t;}
};

/******************************************************************************************/

enum class Ensemble : std::uint_fast8_t {nostacking, stacking, dangle, coaxial, min, all, none};

inline constexpr bool has_subensemble(Ensemble e) {
    switch (e) {
        case Ensemble::nostacking: return false;
        case Ensemble::stacking: return true;
        case Ensemble::dangle: return true;
        case Ensemble::coaxial: return true;
        case Ensemble::min: return false;
        case Ensemble::all: return false;
        case Ensemble::none: return false;
    }
}

extern std::array<char const *, 7> EnsembleNames;


static std::array<Ensemble, 7> const AllEnsembles
    = {Ensemble::nostacking, Ensemble::stacking, Ensemble::dangle, Ensemble::coaxial, Ensemble::min, Ensemble::all, Ensemble::none};

Ensemble as_ensemble(string_view s);

inline std::ostream& operator<<(std::ostream &os, Ensemble e) {return os << (uint(e) < 7 ? EnsembleNames[uint(e)] : "invalid");}

/******************************************************************************************/

/// Historical dangles=some
struct MinDangles {
    NUPACK_REFLECT_EMPTY(MinDangles);
    static auto can_backtrack() {return True();}
    template <class T>
    constexpr T combine(T e1, T e2, bool same) const {return same ? std::min(e1, e2) : e1 + e2;}
};

/// Historical dangles=all
struct AllDangles {
    NUPACK_REFLECT_EMPTY(AllDangles);
    static auto can_backtrack() {return True();}
    template <class T>
    constexpr T combine(T e1, T e2, Ignore) const {return e1 + e2;}
};

/// New dangles=none
/// Same as historical dangles=none except closing dangle issue
struct NoStacking {
    NUPACK_REFLECT_EMPTY(NoStacking);
    static auto can_backtrack() {return True();}
    template <class T>
    constexpr T combine(T, T, Ignore) const {return T();}
};

/// New dangles=coax
struct Stacking {
    NUPACK_REFLECT_EMPTY(Stacking);
    static auto can_backtrack() {return True();}
    template <class T>
    constexpr T combine(T, T, Ignore) const {NUPACK_ERROR("Incorrect ensemble setting (stacking)."); return T();}
};

// Mixed Materials
struct MixedNoStacking {static auto can_backtrack() {return True();}};
struct MixedStacking {static auto can_backtrack() {return True();}};

using EnsembleType = Variant<NoStacking, Stacking, MinDangles, AllDangles>;

inline EnsembleType ensemble_variant(Ensemble e) {
    switch (e) {
        case Ensemble::nostacking: return NoStacking();
        case Ensemble::stacking: return Stacking();
        case Ensemble::dangle: return Stacking();
        case Ensemble::coaxial: return Stacking();
        case Ensemble::none: return NoStacking();
        case Ensemble::min: return MinDangles();
        case Ensemble::all: return AllDangles();
    }
}

/******************************************************************************************/

template <class Dangles, class P>
struct DangleFunction : Dangles {
    P const *p;
    DangleFunction(P const &pset) : p(&pset) {}

    unsigned alphabet_length() const {return p->alphabet.length();}
    auto energy5(Base i, Base j, Base k) const {return (*p)(p->dangle5(), +i, +j, +k);}
    auto energy3(Base i, Base j, Base k) const {return (*p)(p->dangle3(), +i, +j, +k);}
};

template <class P>
struct DangleFunction<NoStacking, P> : NoStacking {
    DangleFunction(P const &pset) {}
    unsigned alphabet_length() const {return 0;}
    value_type_of<P> energy5(Base, Base, Base) const {return value_type_of<P>();}
    value_type_of<P> energy3(Base, Base, Base) const {return value_type_of<P>();}
};

template <class D, class P>
DangleFunction<D, P> dangle_function(D, P const &p) {return {p};}

/******************************************************************************************/

template <class M, class V>
value_type_of<M> stacking_energy(NoStacking, M const &, V const &, int const) {return value_type_of<M>();}

template <class M, class V>
value_type_of<M> stacking_energy(Stacking, M const &model, V const &v, int const nick) {
    auto pf = stacking_sum<real>(v, nick, model);
    return inverse_boltzmann(model.beta, pf);
}

// A list of sequences. nick indicates which one the strand break is before
// For example, the starting loop always has a nick of 0
template <class D, class M, class V, NUPACK_IF(is_same<D, MinDangles, AllDangles>)>
value_type_of<M> stacking_energy(D dangle, M const &model, V const &v, int const nick) {
    using T = value_type_of<M>;
    if (len(v) == 1) return value_type_of<M>();
    auto const si = begin_of(v), sf = end_of(v) - 1;

    auto get = [&, sn=begin_of(v) + nick](auto s, auto t, auto u) {
        T e5, e3;
        if (sn != t) e5 = model.dangle5(back(*s), front(*t), front(*t, 1));
        if (sn != u) e3 = model.dangle3(back_index(*t, 1), back(*t), front(*u));
        if (sn == t) return e3;
        if (sn == u) return e5;
        return dangle.combine(e5, e3, 3 == len(*t));
    };

    T en = T();

    if (len(*si) != 2) en += get(sf, si, si + 1);
    for (auto s = begin_of(v) + 1; s != end_of(v) - 1; ++s)
        if (len(*s) != 2) en += get(s - 1, s, s + 1);
    if (len(*sf) != 2) en += get(sf - 1, sf, si);

    return en;
}

template <class State, class Model, class F>
void for_stacking_states(State const &w, Model const &m, F &&f) {
    vec<int> stack(w.sys.n_bases());
    if (has_subensemble(m.ensemble)) for_stacking_states_impl(w, m, f, stack, 0, 0);
    else f(stack, w.calculate_energy(m));
}

/******************************************************************************************/


/// Dangle at the beginning of an edge, return 0 if s == sn
template <class D, class V, class It>
real safe_dangle5(D const &dangle, V const &v, It s) {
    if (len(*s) <= 2 || front(*s) == Base::null()) return 0;
    auto const lo = cyclic_prev(v, s);
    return dangle.energy5(back(*lo), front(*s), front(*s, 1));
}

/******************************************************************************************/

/// Dangle at the end_of of an edge, return 0 if s == sn - 1
template <class D, class V, class It>
real safe_dangle3(D const &dangle, V const &v, It s) {
    if (len(*s) <= 2 || back(*s) == Base::null()) return 0;
    auto const up = cyclic_next(v, s);
    return dangle.energy3(back_index(*s, 1), back(*s), front(*up));
}

/******************************************************************************************/

}
