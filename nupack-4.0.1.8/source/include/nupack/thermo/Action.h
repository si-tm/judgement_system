/**
 * @brief Recursion adapter for modifying QB recursion, used by design mostly
 *
 * @file Action.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "../algorithms/TypeSupport.h"
#include "../types/Complex.h"
#include <functional>

namespace nupack::thermo {

/******************************************************************************************/

// Action is called and returns a log factor
// the pair partition function is yielded times exp(the log factor)
// short-circuiting takes place if the log factor is -inf.

/**
 * Actions are called once for each base i, j during the calculation.
 * The boolean can_pair whether the thermodynamic model permits a pair between i and j
 * The f is the recursion yielding the value of QB
 * The third argument is the slice of sequences being calculated
 * The fourth and fifth arguments are the indices of i and j
 * The last argument is the partition function matrices block Q
 */
struct Action {
    enum kind {forbid, flat, bonus};
    std::function<std::pair<real, kind>(ComplexView, int, int)> rule; // treated as if it returns true if empty

    template <class Block, class Algebra, class F, class Model>
    auto operator()(int i, int j, bool can_pair, Algebra A, Block const &, ComplexView strands, Model const &mod, F &&recursion) const noexcept {
        auto [factor, k] = rule ? std::pair<value_type_of<Model>, kind>(rule(strands, i, j))
                                : std::pair<value_type_of<Model>, kind>(0, bonus);

        return A.sum(k == bonus && can_pair ? A.nullable(A.product(recursion(), mod.boltz(factor))) : A.zero(), // do recursion
                                  k == flat ? mod.boltz(factor) : mod.zero()); // apply flat term or zero
    }
};

/******************************************************************************************/

}
