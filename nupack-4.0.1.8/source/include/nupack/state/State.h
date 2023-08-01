#pragma once
#include "StaticState.h"
#include "../model/Model.h"

namespace nupack {

/******************************************************************************************/

template <class Model>
auto structure_energy(System s, PairList p, Model const &em, bool distinguishable=true, bool stack_energy=false) {
    StaticState<> w(std::move(s), std::move(p));
    w.check_structure(em.pairing());
    auto out = w.calculate_energy(em);

    if (!distinguishable) out += std::log(real(w.symmetry())) / em.beta;

    return out;
}

/******************************************************************************************/

template <class Model>
auto minimum_stack_energy(System s, PairList p, Model const &em) {
    StaticState<> w(std::move(s), std::move(p));
    w.check_structure(em.pairing());

    auto e = sum(w.loops, [&](auto const &o) {return em.minimum_stack_energy(o.sequences(), o.nick());});
    e += (len(w.sys.strands()) - len(w.complexes)) * em.join_penalty();
    return e;
}


}
