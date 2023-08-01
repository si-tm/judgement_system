/**
 * @brief Equilibrium toolkit functionality
 *
 * @file Equilibrate.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#define ARMA_DONT_PRINT_ERRORS
#ifndef ARMA_WARN_LEVEL
#   define ARMA_WARN_LEVEL 1
#endif
#include <armadillo>
#include <vector>

namespace nupack { namespace concentration {

using real = double;
using uword = arma::uword;
namespace la = arma;

enum class Method : unsigned int {fit, cd, dogleg};

struct Options {
    la::Col<real> initial_dual_guess;
    uword max_iters = 1e4;
    real tolerance = 1.e-8;
    real delta_min = 1.e-12;
    real delta_max = 1000.0;
    bool orthogonalize = true; //< takes care of cases where # strands > # complexes
    Method method = Method::cd; //< usually the non-uniform method is better
};

la::Col<real> initial_dual_guess(la::Mat<real> const &A, la::Col<real> const &logq, la::Col<real> const &x0);

template <class T>
struct Output {
    la::Col<T> solution, dual_solution;
    real objective, error;
    uword iters = 0;
    bool converged = false;
};

/**
 * @brief Solve equilibrium concentrations
 * @param A Coefficient matrix (complexes, strands)
 * @param logb initial LOG of STRAND concentrations
 * @param q log partition functions
 * @param ops solving options
 * @return V equilibrated concentrations
 */
Output<real> equilibrate(la::Mat<real> const &A, la::Col<real> logb, la::Col<real> const &q, Options const &ops={});

/**
 * @brief Solve equilibrium concentrations for complexes
 * indices: list of ordered indices of strands
 * logq: list of log partition functions (distinguishable)
 * s0: strand concentrations
 */

Output<real> solve_complexes(std::vector<std::vector<uword>> const &indices, la::Col<real> logq, la::Col<real> s0, Options const &ops, bool rotational_correction=true, bool as_strands=true);

}}
