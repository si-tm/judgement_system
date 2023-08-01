/**
 * @brief Equilibrium toolkit algorithms - usually it will suffice to include Equilibrate.h
 *
 * @file Concentration.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "Equilibrate.h"
#include "Utility.h"

// #include "../common/Random.h"
// #include "../math/BoundSolve.h"
// #include "../iteration/Search.h"
// #include "../reflect/SerializeMatrix.h"
// #include "../reflect/Serialize.h"
// #include <spdlog/spdlog.h>
#include <fmt/ostream.h>
#include <stdexcept>
#include <sstream>

namespace nupack { namespace concentration {

/******************************************************************************************/

// where x is a vector and H is symmetric, return A.T * diag(x) * A
// this is probably the rate limiting step, complexity is (# strands)^2 (# complexes)
template <class M, class V>
void symmetric_mdm(M &H, M const &A, V const &x) {
    H.set_size(A.n_cols, A.n_cols);
    for (uword j = 0; j != H.n_cols; ++j)
        H.col(j).head(j+1) = A.head_cols(j+1).t() * (x % A.col(j));

    for (uword j = 0; j != H.n_cols; ++j)
        for (uword i = 0; i != j; ++i)
            H.at(j, i) = H.at(i, j);
}


/******************************************************************************************/

// x_{c, I} = \sum_j (A_{j, c} - A_{j, c_i}) \lambda_j + q_c - q_{c_i}  \\
// y_{c, I} = A_{c, I} e^{x_c, I}  \\
// \log b - \log a - q = \hat{A}^T \lambda + \log ( 1^T Y ) + \frac{ (\lambda - \lambda^0)^T (A^T Y - \hat{A} d(1^T Y)) } { 1^T Y  }

template <class T, class M>
Output<T> equilibrate_lse(M const &A, la::Col<T> const &logb, la::Col<T> const &logq, Options ops) {
    if (!logb.is_finite()) throw std::invalid_argument("Strand concentrations should be finite");
    if (!logq.is_finite()) throw std::invalid_argument("Partition functions should be finite");
    if (A.min() < 0) throw std::invalid_argument("Coefficient matrix should be non-negative");

    Output<T> out;
    auto &y = out.dual_solution = la::solve(A, la::log(A * la::exp(logb)) - logq);
    auto &x = out.solution;
    la::uvec c(logb.n_cols);
    M X(A.n_rows, A.n_cols), Y, Alog = la::log(A);
    la::Mat<T> G;
    la::Col<T> rhs, error, xc, shift;

    for (out.iters = 0; out.iters != ops.max_iters; ++out.iters) {
        x = A * y + logq;
        c = la::index_max(Alog.each_col() + x, 0).t();
        xc = -(A.rows(c) * y + logq(c));
        X.each_col() = x;
        X.each_row() += xc.t();
        Y = A % la::exp(X);

        rhs = la::sum(Y, 0).t();
        error = logb + xc - la::log(rhs);
        // print(out.iters, la::max(la::abs(error)), json(c));
        if (la::max(la::abs(error)) < ops.tolerance) {out.converged = true; break;}

        G = Y.t() * A;
        rhs %= error;
        y += la::solve(G, rhs); // la::solve_opts::refine
        x = A * y + logq;
    }
    x = la::exp(x);
    return out;
}

/******************************************************************************************/

// This one converges nicely but isn't very guarded against exponentiation overflow
template <class T, class M>
Output<T> equilibrate_cd(M const &A, la::Col<T> const &logb, la::Col<T> const &logq, Options ops) {
    if (!A.is_finite()) throw std::runtime_error("non-finite equilibration coefficients");
    Output<T> out;
    // auto &y = out.dual_solution = la::solve(A, la::log(A * la::exp(logb)) - logq);
    auto &y = out.dual_solution = ops.initial_dual_guess;
    la::Mat<T> H, V, AV;
    la::Col<T> e, logx1, x1, logx = A * y + logq;
    auto &x = out.solution = la::exp(logx);
    if (!x.is_finite()) throw std::runtime_error("non-finite values in initial guess");// , x, logx, vy, shift, s, s0, h, AV.col(i));
    la::Col<T> err1, err = la::log(A.t() * x) - logb;
    la::Col<T> const x0 = la::solve(A.t(), la::exp(logb));

    for (out.iters = 0; out.iters != ops.max_iters; ++out.iters) {
        if (!x.is_finite()) throw std::runtime_error("non-finite values encountered in equilibration");// , x, logx, vy, shift, s, s0, h, AV.col(i));
        out.error = la::max(la::abs(err));
        if (out.error < ops.tolerance) {out.converged = true; break;}

        if (out.iters % logb.n_rows == 0) { // O(m n^2), done every n is O(m n)
            symmetric_mdm(H, A, x);
            if (!la::eig_sym(e, V, H)) {
                if (!H.is_finite()) throw std::runtime_error("non-finite Hessian in equilibrium concentration solver");
                throw std::runtime_error("eigendecomposition failed in equilibrium concentration solver");
            }
            if (la::abs(V.t() * V - la::eye(V.n_cols, V.n_cols)).max() > 1e-10) {
                la::Mat<std::complex<real>> V0;
                la::Col<std::complex<real>> e0;
                if (!la::eig_gen(e0, V0, H))
                    throw std::runtime_error("general eigendecomposition failed in equilibrium concentration solver");
                e = la::real(e0);
                V = la::real(V0);
            }
            AV = A * V;
        }

        for (la::uword i = 0; i != logb.n_rows; ++i) { // O(m n) overall
            T const s = la::dot(AV.col(i), x); // O(m)
            T const s0 = la::dot(AV.col(i), x0); // O(m)
            T h = la::accu(AV.col(i) % AV.col(i) % x); // O(m)
            if (h <= 0) continue;

            // Newton step in this DOF, assuming the Hessian is non-singular
            T shift = (s0 - s) / h;
            if (!std::isfinite(shift) || shift == 0) continue;

            // A mini-bisection method that usually terminates in one iteration (just until a downhill step is found)
            while (true) {
                // Try new solution
                logx1 = logx + shift * AV.col(i); // O(m)
                x1 = la::exp(logx1); // O(m)
                err1 = la::log(A.t() * x1) - logb; // O(m)
                // If it improves, take the step.
                if (la::norm(err1) <= la::norm(err)) {
                    logx1.swap(logx);
                    x1.swap(x);
                    err1.swap(err);
                    y += V.col(i) * shift;
                    break;
                }
                // Otherwise shrink the step and try again
                shift /= 2;
            }
        }
    }
    return out;
}

/******************************************************************************************/

/// Using linear solver, gradient, Hessian, delta, and minimum delta return the dogleg direction to go to
template <class V, class M, class T>
V find_direction(V const &grad, M const &Hess, T const delta, T const min_delta) {
    // Probably could optimize the order of Newton, Cauchy evaluations below
    // Calculate Newton step with a SPD, singular accepting solver
    V newt = -grad;
    bool const newton = la::solve(newt, Hess, newt,
#if (ARMA_VERSION_MAJOR >= 9) && (ARMA_VERSION_MINOR >= 500)
        la::solve_opts::fast + la::solve_opts::allow_ugly + la::solve_opts::likely_sympd
#else
// #       warning "Armadillo 9.500 or newer is recommended."
        la::solve_opts::fast + la::solve_opts::allow_ugly
#endif
    );
    T const newt_norm = norm(newt);

    // Take Newton if we are inside the minimum trust region
    if (newton && newt_norm > 0 && (delta < min_delta || newt_norm < delta))
        return newt;

    // Calculate Cauchy step
    V cauchy = grad / norm(grad);
    cauchy = grad * (-1 / dot(cauchy, Hess * cauchy));
    T const cauchy_norm = norm(cauchy);

    // Take Cauchy if Newton failed or we are outside the trust region
    if (!newton || newt_norm == 0 || !std::isfinite(sq(newt_norm)) || cauchy_norm > delta)
        return sqrt(delta / cauchy_norm) * cauchy;

    // Dogleg - take this if we are in intermediate region
    T const newt_cauchy = dot(newt, cauchy);

    auto q = quadratic_solve<T>(sq(newt_norm) + sq(cauchy_norm) - 2 * newt_cauchy,
                                2 * (newt_cauchy - sq(cauchy_norm)),
                                sq(cauchy_norm) - sq(delta));

    auto beta = std::abs(q.first) <= std::abs(q.second) ? q.first : q.second; // choose correct root for mixing coefficient
    if (!(std::abs(beta) <= 1)) throw std::runtime_error("invalid dogleg step");

    if (beta < 0) return (beta + 1) * cauchy;
    else return (1 - beta) * cauchy + beta * newt;
}

/******************************************************************************************/

template <class V, class P, class O, class D>
struct DualSystem {
    V grad;
    la::Mat<typename V::elem_type> hess;

    DualSystem(V v, P p, O o, D d) : dual(std::move(v)), to_primal(std::move(p)), objective_function(std::move(o)), derivatives(std::move(d)) {
        to_primal(primal, dual);
        obj = objective_function(primal, dual);
        derivatives(*this);
    }

    V primal, dual;
    real obj;

    P to_primal; // (primal, dual) -> primal
    O objective_function;  // (primal, dual) -> float
    D derivatives;

    void calculate_gradient_and_hessian() {derivatives(*this);}

    auto const &hessian() const {return hess;}

    auto const &gradient() const {return grad;}

    auto objective() const {return obj;}

    // Take another system and set self to its dual value plus a shift
    template <class T>
    void update(DualSystem const &s, T &&shift) {
        dual = s.dual + std::forward<T>(shift);
        to_primal(primal, dual);
        obj = objective_function(primal, dual);
    }

    void swap(DualSystem &x) {
        std::swap(primal, x.primal);
        std::swap(dual, x.dual);
        std::swap(obj, x.obj);
    }

    bool valid() const {
        return primal.is_finite() && dual.is_finite() && std::isfinite(obj);
    }
};

template <class ...Ts>
auto dual_system(Ts &&...ts) {return DualSystem<std::decay_t<Ts>...>(std::forward<Ts>(ts)...);}

/******************************************************************************************/

/// Converge a system using a trust-region method
/// Provide initial system and delta, functions to compute gradient, offset, direction, convergence, and radius adjustment
template <class Sys, class D, class C, class A>
Sys trust_region(Sys s, D &&direction, C &&condition, A &&adjust_delta, real delta) {
    auto s2 = s;
    while (true) {
        if (condition(s, s.gradient(), s.hessian())) break;
        auto const p = direction(s.gradient(), s.hessian(), delta);
        s2.update(s, p);

        // rho measures actual improvement divided by expected improvement (if exact, = 1)
        auto const expected = dot(s.gradient(), p) + dot(p, s.hessian() * p) / 2; // should almost always be < 0
        auto const rho = (s2.objective() - s.objective()) / expected; // better for rho to be >0

        if (std::isnan(rho)) throw std::runtime_error("trust region solver encountered NaN value");

        // Possibly accept the new values. Also, adjust the trust region size
        delta = adjust_delta(delta, rho);
        if (s2.objective() <= s.objective()) {
            s.swap(s2);
            s.calculate_gradient_and_hessian();
        }

        // Check validity of current solution
        if (!s.valid()) throw std::runtime_error("trust region solver encountered non-finite value");
    }
    return s;
}

/******************************************************************************************/

/**
 * @brief Solve equilibrium concentrations
 * The objective is O(y) = 1^T exp(A y + q) + y^T A^T x0
 * @param A Coefficient matrix (complexes, strands)
 * @param x0 initial concentrations of complexes
 * @param q log partition functions
 * @param ops solving options
 * @return V equilibrated concentrations
 */
template <class T>
Output<T> equilibrate_gradient(la::Mat<T> A, la::Col<T> const &s0, la::Col<T> const &q, Options const &ops) {
    la::Mat<T> orth_A;
    if (ops.orthogonalize) {
        orth_A = orth(A.t());
        A = A * orth_A;
    }

    // Direction finder with minimum radius delta
    auto direction = [dmin=ops.delta_min](auto const &...ts) {
        auto p = find_direction(ts..., dmin);
        // It is unknown why this can happen, but it does in some irreproducible design cases
        for (auto &x : p) if (!std::isfinite(x)) x = 0;
        return p;
    };

    uword n = 0; // number of condition checks (= 1 + number of iterations)
    real error;
    bool good = false;
    // Negative total concentrations of each strand type in the chosen basis
    la::Col<T> const rhs = -(ops.orthogonalize ? orth_A.t() * s0 : s0);
    // Inverse total concentrations of each strand
    la::Col<T> const normalization = 1 / s0;
    if (!normalization.is_finite() || normalization.min() <= 0) throw std::runtime_error("invalid normalization");
    // Convergence criterion based on gradient norm
    auto condition = [&] (auto const &sys, auto const &grad, auto const &...) {
        if (!ops.orthogonalize) error = la::max(la::abs(grad) % normalization);
        else error = la::max(la::abs(orth_A * grad) % normalization);
        good = error < ops.tolerance;
        return (ops.max_iters < ++n) || good;
    };
    // Function to calculate primal vector from previous primal and dual
    auto primal = [&](auto &x, auto const &y) {
        // In the future, 1e100 is a bit arbitrary.
        // This is done to prevent infinity from occurring during intermediate solution
        // It would be better to constrain the dual vector y to be in a reasonable regime.
        x = la::clamp(la::exp(A * y + q), std::numeric_limits<real>::min(), 1e100);
        if (!x.is_finite()) throw std::runtime_error("primal transformation failed"); //NUPACK_ASSERT(x.is_finite(), x, A, y, q, A * y + q, n);
    };
    // Function to calculate objective from primal and dual
    auto objective = [&](auto const &x, auto const &y) {return accu(x) + dot(y, rhs);};

    // Function to calculate gradient in dual space from system
    auto derivatives = [&](auto &s) {
        // Function to calculate gradient in dual space from system
        s.grad = A.t() * s.primal + rhs;
        if (!s.grad.is_finite()) throw std::runtime_error("gradient is non-finite");
        // this is probably the rate limiting step, complexity is (# strands)^2 (# complexes)
        symmetric_mdm(s.hess, A, s.primal); // / same as H = A.t() * la::diagmat(s.primal) * A;
    };

    auto sys = dual_system(ops.orthogonalize ? orth_A.t() * ops.initial_dual_guess : ops.initial_dual_guess,
                           primal, objective, derivatives);

    // Function to adjust trust region radius and return if the new values should be accepted
    auto adjust_delta = [=](auto const delta, auto const rho) {
        if (delta <= ops.delta_min) return delta;
        if (rho > 0.75) return std::min(2 * delta, ops.delta_max);
        if (rho < 0.25) return delta / 4;
        return delta;
    };
    // Run trust region with initial radius 1000
    auto ret = trust_region(std::move(sys), direction, condition, adjust_delta, ops.delta_max);

    // Undo the threshold that was put on in primal()
    ret.primal.replace(std::numeric_limits<real>::min(), 0);

    if (ops.orthogonalize) ret.dual = orth_A * ret.dual;
    return {std::move(ret.primal), std::move(ret.dual), ret.objective(), error, n - 1, good};
}

/******************************************************************************************/

template <class T>
Output<T> equilibrate_finite(la::Mat<T> const &A, la::Col<T> const &logb, la::Col<T> const &logq, Options ops) {
    if (logq.empty()) return {{}, {}, 0, 0, 0, true};

    if (ops.initial_dual_guess.empty())
        ops.initial_dual_guess = initial_dual_guess(A, logq, la::Col<T>(A * la::exp(logb)));

    if (logq.n_rows != A.n_rows)
        throw std::runtime_error("equilibrate_finite(): inconsistent number of complexes");
    if (A.n_cols != ops.initial_dual_guess.n_rows || A.n_cols != logb.n_rows)
        throw std::runtime_error("equilibrate_finite(): inconsistent number of strands");

    switch (ops.method) {
        case Method::cd: return equilibrate_cd<T>(A, logb, logq, ops);
        case Method::fit: return equilibrate_lse<T>(A, logb, logq, ops);
        case Method::dogleg: return equilibrate_gradient(A, la::Col<T>(la::exp(logb)), logq, ops);
    }
}

/******************************************************************************************/

template <class T>
la::Col<T> initial_guess(la::Mat<T> const &A, la::Col<T> const &logq, la::Col<T> const &x0) {
    if (!(x0.min() > 0)) throw std::runtime_error("initial complex concentrations must be positive");
    la::Col<T> weight(logq.n_rows, la::fill::ones);
    la::Mat<T> const logA = la::log(A);
    la::Col<T> y, logx;

    // Refocus least squares until all are finite
    while (true) {
        la::Mat<T> const AwA = A.t() * la::diagmat(weight) * A;
        if (!la::solve(y, AwA, A.t() * (weight % (la::log(x0) - logq)))) {
            std::stringstream ss;
            ss << "initial solution failed\n";
            ss << "y\n" << y.t() << "\n";
            ss << "w\n" << weight.t() << "\n";
            ss << "AwA\n" << AwA << "\n";
            ss << "logq\n" << logq.t() << "\n";
            ss << "x0\n" << x0.t() << "\n";
            throw std::runtime_error(ss.str());
        }
        logx = A * y + logq;
        bool bad = false;
        for (la::uword i = 0; i != A.n_cols; ++i) {
            la::Col<T> logc = logA.col(i) + logx;
            auto j = logc.index_max();
            if (logc(j) < -200 || logc(j) > 200) {bad = true; weight(j) *= 2;}
        }
        if (!bad) return y;
    }
}

}}
