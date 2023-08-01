#pragma once
#include "../types/Matrix.h"

namespace nupack {

/******************************************************************************************/

template <class T>
real sparsity(T const &t) {return 1 - std::count(begin_of(t), end_of(t), 0) / real(len(t));}

/******************************************************************************************/

/// Make a matrix from a list of tuples where it is nonzero (x, y, value)
template <class Mat, class Tuples, NUPACK_IF(la::is_dense<Mat>)>
Mat matrix_from_tuples(usize n_rows, usize n_cols, Tuples const &tups) {
    Mat ret(n_rows, n_cols); la::fill_zero(ret);
    for (auto const &t : tups) ret(first_of(t), second_of(t)) += third_of(t);
    return ret;
};

/******************************************************************************************/

template <class Mat, class Tuples, NUPACK_IF(la::is_sparse<Mat>)>
Mat matrix_from_tuples(usize const n_rows, usize const n_cols, Tuples const &tups) {
    la::vec values(len(tups));
    la::umat locs(2, len(tups));
    for (auto i : indices(tups)) {
        locs(0, i) = first_of(tups[i]);
        locs(1, i) = second_of(tups[i]);
        values(i) = third_of(tups[i]);
    }
    return {locs, values, n_rows, n_cols};
}

/******************************************************************************************/

struct Sparsity {
    real threshold = 0;
    std::size_t row_size = 0;
    bool clamp = true, diagonal_only=false;
    NUPACK_REFLECT(Sparsity, threshold, row_size, clamp, diagonal_only);
};

/******************************************************************************************/

template <class T>
struct PairMatrix {
    Col<T> values, diagonal;
    Col<std::uint32_t> rows, cols;

    PairMatrix() = default;

    PairMatrix(Col<T> d, Col<T> v, Col<std::uint32_t> r, Col<std::uint32_t> c)
        : values(std::move(v)), diagonal(std::move(d)), rows(std::move(r)), cols(std::move(c)) {
            NUPACK_REQUIRE(len(values), ==, len(rows));
            NUPACK_REQUIRE(len(values), ==, len(cols));
        }

    explicit PairMatrix(Mat<T> m, Sparsity const &sp={});

    Mat<T> full() const {
        Mat<T> out = la::diagmat(diagonal);
        zip(rows, cols, values, [&](auto i, auto j, auto v) {out(j, i) = out(i, j) = v;});
        return out;
    }

    T operator()(std::uint32_t i, std::uint32_t j) const {
        if (i == j) return diagonal(i);
        using P = std::pair<std::uint32_t, std::uint32_t>;
        P const p  = std::minmax(i, j);
        auto const idx = indices(values);
        auto const it = std::lower_bound(idx.begin(), idx.end(), p, 
            [&](auto k, P const &p) {return P(rows.at(k), cols.at(k)) < p;});
        return it != idx.end() && rows.at(*it) == p.first && cols.at(*it) == p.second ? values.at(*it) : T();
    }

    template <class V>
    T defect(V const &v) const {
        NUPACK_REQUIRE(len(v), ==, len(diagonal), "Dimension mismatch");
        T out = 0;
        izip(diagonal, [&](auto i, auto x) {
            if (v[i] == i) out += x;
        });
        zip(rows, cols, values, [&](auto i, auto j, auto x) {
            if (v[i] == j) out += 2 * x;
        });
        return len(diagonal) - out;
    }

    // void normalize() {
    //     1 P = 1
    //     P 1 = 1
    //     diagonal.fill(0);
    //     zip(rows, cols, values, [](auto i, auto j, auto x) {
    //         diagonal(i) += x;
    //         diagonal(j) += x;
    //     });
    //     diagonal = 1 - diagonal;
    //     // calc totals of BPs for each row.
    //     // if smaller than 1, add to unpaired
    //     // if bigger than 1, ...
    //     // diagonal will be negative ... need to scale the pairs accordingly
    // }

    NUPACK_REFLECT(PairMatrix, values, diagonal, rows, cols);
};

template <class P, class F>
void visit_pair_matrix(P &&p, bool sym, F &&f) {
    izip(p.diagonal, [&](la::uword i, auto &&x) {f(i, i, x);});
    zip(p.rows, p.cols, p.values, [&](la::uword i, la::uword j, auto &&x) {f(i, j, x);});
    if (sym) zip(p.rows, p.cols, p.values, [&](la::uword i, la::uword j, auto &&x) {f(j, i, x);});
}

template <class T>
void clamp_pair_matrix(Mat<T> &m) {
    m.diag().zeros();
    for (auto i : range(m.n_rows)) {
        real const s = la::accu(m.col(i));
        NUPACK_QUICK_REQUIRE(s, <, 1.01);
        if (s <= 1) {
            m.at(i, i) = 1 - s;
        } else { // do some clamping to avoid tiny disobeying of the [0, 1] domain
            m.col(i) *= 1/s;
            m.row(i) *= 1/s;
            m.at(i, i) = 0;
        }
    }
}

/// Take square matrix of probabilities
/// Return diagonal entries, off-diagonal entries, rows, columns
/// if n is 0, return all entries, else return the n highest non-diagonal entries in each row
// complexity is N^2 log(n) + N n log(n N) \approx N^2 log(n) + N n log(N)
template <class T>
PairMatrix<T>::PairMatrix(Mat<T> m, Sparsity const &sp) {
    NUPACK_REQUIRE(m.n_rows, ==, m.n_cols);
    NUPACK_ASSERT(m.is_symmetric(), "sparse_pair_matrix(): pair matrix should be symmetric", m);

    if (sp.clamp) clamp_pair_matrix(m);
    diagonal = m.diag();
    if (sp.diagonal_only) return;

    if (m.n_rows <= 1) return;

    auto const filter = [t0=sp.threshold](T const &t) {return t > t0;};
    auto const nnz = (count_if(m, filter) - count_if(diagonal, filter)) / 2;
    bool const simple = (sp.row_size == 0) || (sp.row_size >= m.n_rows/2) || (nnz <= sp.row_size * m.n_rows);

    if (simple) {
        for_each(std::tie(values, rows, cols), [nnz](auto &x) {x.set_size(nnz);});
        la::uword p = 0;
        for (auto i : range(m.n_rows)) for (auto j : range(i+1, m.n_rows)) {
            if (filter(m.at(j, i))) {
                rows(p) = i;
                cols(p) = j;
                values(p) = m.at(j, i);
                ++p;
            }
        }
        NUPACK_REQUIRE(p, ==, nnz, "Sparse pair matrix failure", count_if(m, filter), count_if(diagonal, filter));
    } else {
        std::vector<la::uword> idx;
        std::vector<std::pair<la::uword, la::uword>> v;
        v.reserve(m.n_rows * sp.row_size);

        for (auto const j : range(m.n_rows)) {
            // all possible column indices
            idx.assign(range(m.n_rows).begin(), range(m.n_rows).end());
            // partially sort to get the n biggest column indices
            std::nth_element(idx.begin(), idx.begin() + sp.row_size, idx.end(), [j, c=m.col(j)](auto i1, auto i2) {
                if (i1 == j) return false; // don't include on-diagonal elements
                if (i2 == j) return true;
                return c(i1) > c(i2);
            });
            // only include nonzero elements
            for (auto it = idx.begin(); it != idx.begin() + sp.row_size; ++it)
                if (*it != j && filter(m.at(*it, j))) v.emplace_back(std::minmax(*it, j));
        }
        v = unique_sorted(std::move(v));

        for_each(std::tie(values, rows, cols), [s=v.size()](auto &x) {x.set_size(s);});
        zip(v, rows, cols, values, [&m](auto const &p, auto &r, auto &c, auto &v) {
            r = p.first;
            c = p.second;
            v = m.at(p.second, p.first);
        });
    }
}


}
