#include <nupack/concentration/Solve.h>
#include <stdexcept>

namespace nupack { namespace concentration {

la::Col<real> initial_dual_guess(la::Mat<real> const &A, la::Col<real> const &logq, la::Col<real> const &x0) {
    return initial_guess(A, logq, x0);
}

/******************************************************************************************/

Output<real> equilibrate(la::Mat<real> const &A, la::Col<real> logb, la::Col<real> const &logq, Options const &ops) {
    if (logq.empty()) throw std::invalid_argument("No complexes in equilibrium input");

    constexpr real const inf = std::numeric_limits<real>::infinity();
    if (logq.has_nan()) throw std::domain_error("Input log Q contains NaN");
    if (logq.max() == inf) throw std::domain_error("Input log Q contains +inf");

    // Reduce complexes of the same composition
    std::map<la::Col<real>, real, LexicographicalCompare> unique;
    izip(logq, [&](auto i, auto logq) {
        if (logq == -inf) return;
        auto &p = unique.emplace(A.row(i).t(), -inf).first->second;
        p = log_sum_exp(p, logq);
    });

    la::Mat<real> A2(unique.size(), A.n_cols);
    la::Col<real> logq2(unique.size());
    izip(unique, [&](auto i, auto const &p) {
        A2.row(i) = p.first.t();
        logq2(i) = p.second;
    });

    auto out = equilibrate_finite<real>(A2, logb, logq2, ops);
    out.solution = la::exp(A * out.dual_solution + logq);
    return out;
}

/******************************************************************************************/

Output<real> solve_complexes(std::vector<std::vector<uword>> const &indices, la::Col<real> logq, la::Col<real> x0, Options const &ops, bool rotational_correction, bool as_strands) {
    if (indices.size() != logq.size()) throw std::invalid_argument("Inconsistent number of complexes");
    uword const n_strands = maxmap(indices, 0, [](auto const &p) {return maxmap(p, 0, [](auto i) {return i;}) + 1;});
    la::Mat<real> A(logq.size(), n_strands, arma::fill::zeros);

    izip(indices, [&](auto c, auto const &x) {
        for (auto s : x) A(c, s) += 1;
        if (rotational_correction) logq(c) -= std::log(real(rotational_symmetry(x)));
    });

    if (!(x0.min() > 0)) throw std::invalid_argument("All concentrations should be positive");

    if (as_strands) {
        if (x0.size() != A.n_cols) throw std::invalid_argument("Incorrect number of concentrations given");
        x0 = la::log(std::move(x0));
    } else {
        if (x0.size() != A.n_rows) throw std::invalid_argument("Incorrect number of concentrations given");
        x0 = la::log(A.t() * x0);
    }

    return equilibrate(std::move(A), std::move(x0), std::move(logq), ops);
}

/******************************************************************************************/

}}
