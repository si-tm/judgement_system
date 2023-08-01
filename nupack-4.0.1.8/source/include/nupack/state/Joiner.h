#pragma once
#include "../loop/LoopJoiner.h"
#include "StateComplexes.h"

namespace nupack::kmc {

/******************************************************************************************/

struct ProductJoiner : TotallyOrdered {
    using Index = iseq;
    using Matrix = Mat<real>;
    using Column = Col<real>;

    /// Loop join data
    vec<Matrix> join_rates; // join rates for each exterior loop
    /// Complex join data
    Fenwick<Matrix> complex_rates; // all join propensities
    Fenwick<Column> complex_rates_squared; // all self join propensities
    BasePairList base_pairs;

    NUPACK_REFLECT(ProductJoiner, join_rates, complex_rates, complex_rates_squared, base_pairs);

    /**************************************************************************************/

    // Comparing is hard ... just make all equal to be simpler.
    bool operator<(ProductJoiner const &) const {return false;}
    bool operator==(ProductJoiner const &) const {return true;}

    template <class State>
    void initialize_joiner(State const &w);

    /// Total join rate for state
    template <class RF> 
    real join_rate(RF const &) const;

    /// Update join rates for strand "i"
    void update_join_rates_matrix(StateComplexes const &complexes, Index i, Matrix const &);
    
    template <class State, class O>
    void update_join_rates(State const &state, Index i, O const &loop) {update_join_rates_matrix(state.complexes, i, loop.join_propensities());}

    /// Choose a join move
    template <class State, class RF>
    ComplexJoinMove get_join_move(State const &complexes, real, RF const &) const;
    /// Choose a join move
    JoinMove get_join_move_nondimensional(StateComplexes const &, real const) const;
    /// Register that strands i and j have joined
    void register_join(StateComplexes const &, Index i, Index j);
    /// Register that strands i and j have split
    template <class State>
    void register_split(State &, Index i, Index j);

    auto loop_joiner() const {return ProductLoopJoiner();}

    template <class T, NUPACK_IF(std::is_same_v<T, Matrix>)>
    Col<real> calculate_rate_product(T const &t) const {return t.col(0) % t.col(1);}

    template <class M>
    std::size_t length(M const &model) const {
        if (has_subensemble(model.ensemble)) {
            std::size_t n_bases = model.alphabet().length(), n_dangles = n_bases + 1;
            return n_bases * n_dangles * n_dangles + 2 * n_bases + 4 * n_bases * n_dangles + 2 * len(base_pairs);
        } else return len(base_pairs);
    }
};

/******************************************************************************************/

struct QuadraticJoiner : MemberOrdered {
    using Index = iseq;

    Fenwick<real> join_rates; // triangular array for every ordered pair of strands
    vec<std::pair<Index, Index>> strand_pairs; // pair of strand indices for each join rate

    NUPACK_REFLECT(QuadraticJoiner, join_rates, strand_pairs);

    template <class State>
    void initialize_joiner(State const &w);

    /// Register that strands i and j have joined
    void register_join(StateComplexes const &, Index i, Index j);
    /// Register that strands i and j have split
    template <class State>
    void register_split(State &, Index i, Index j);
    // Index of a pair of strand indices inside the triangular array
    auto join_index(Index i, Index j) const {NUPACK_QUICK_REQUIRE(i, !=, j); if (j < i) std::swap(i, j); return j * (j - 1) / 2 + i;}
    // Update join rates between revised loop and every other exterior loop in a different complex
    template <class State, class O, class Predicate=AlwaysTrue>
    void update_join_rates(State const &state, Index i, O const &loop, Predicate const &pred={});

    template <class State, class RF, class Predicate=AlwaysTrue>
    ComplexJoinMove get_join_move(State const &state, real, RF const &, Predicate const &pred={}) const;

    template <class RF> 
    real join_rate(RF const &) const {return join_rates.total();}

    auto loop_joiner() const {return QuadraticLoopJoiner();}
};

/******************************************************************************************/

struct NoJoiner {
    NUPACK_REFLECT_EMPTY(NoJoiner);
    
    NoLoopJoiner loop_joiner() const {return {};}

    void register_split(Ignore, Ignore, Ignore) {}
    void register_join(Ignore, Ignore, Ignore) {NUPACK_ERROR("NoJoiner should not be used this way");}

    void initialize_joiner(Ignore) {}
    constexpr real join_rate(Ignore) const {return 0;}

    void update_join_rates(Ignore, Ignore, Ignore) {}

    ComplexJoinMove get_join_move(Ignore, Ignore, Ignore) const {NUPACK_ERROR("No joins possible");}
};

/******************************************************************************************/

// In the future the LoopJoiner should just be moved into Joiner, that was sort of a mistake.
struct Joiner : MemberOrdered {
    using Index = iseq;
    Variant<ProductJoiner, QuadraticJoiner, NoJoiner> impl;

    Joiner() = default;
    Joiner(ProductJoiner j) : impl(std::move(j)) {}
    Joiner(QuadraticJoiner j) : impl(std::move(j)) {}
    Joiner(NoJoiner j) : impl(std::move(j)) {}
    
    NUPACK_REFLECT(Joiner, impl);

    LoopJoiner loop_joiner() const {return fork(impl, [](auto const &j) -> LoopJoiner {return j.loop_joiner();});}

    ProductJoiner const &product_joiner() const {return std::get<ProductJoiner>(impl);}
    QuadraticJoiner const &quadratic_joiner() const {return std::get<QuadraticJoiner>(impl);}

    template <class State, class O>
    void update_join_rates(State const &state, Index i, O const &loop) {
        fork(impl, [&](auto &j) {j.update_join_rates(state, i, loop);});
    }

    /// Register that strands i and j have split
    template <class State>
    void register_split(State &w, Index i, Index j) {
        fork(impl, [&](auto &joiner) {joiner.register_split(w, i, j);});
    }

    // Register join move. This happens before StateComplexes::register_join() is called.
    void register_join(StateComplexes &complexes, Index i, Index j) {
        fork(impl, [&](auto &joiner) {joiner.register_join(complexes, i, j);});
    }

    template <class State>
    void initialize_joiner(State const &w) {
        fork(impl, [&](auto &joiner) {joiner.initialize_joiner(w);});
    }

    template <class RF> 
    real join_rate(RF const &rf) const {
        return fork(impl, [&](auto const &j) {return j.join_rate(rf);});
    }

    template <class State, class RF>
    ComplexJoinMove get_join_move(State const &state, real value, RF const &rf) const {
        return fork(impl, [&](auto const &j) {return j.get_join_move(state, value, rf);});
    }
};

/******************************************************************************************/

template <class RF>
real ProductJoiner::join_rate(RF const &rf) const {
    if (complex_rates.size() < 2) return 0.0;

    return 0.5 * rf.bimolecular_propensity_prefactor() 
        * (la::accu(calculate_rate_product(complex_rates.total())) - la::accu(complex_rates_squared.total()));
}

/******************************************************************************************/

template <class State, class RF>
ComplexJoinMove ProductJoiner::get_join_move(State const &state, double r, RF const &rf) const {
    NUPACK_REQUIRE(r, <, join_rate(rf));
    NUPACK_REQUIRE(r, >=, 0);

    // Choose the loops
    auto const m = get_join_move_nondimensional(state.complexes, 2 * r / rf.bimolecular_propensity_prefactor());
    auto const &p = state.loops[m.o1];
    auto const &k = state.loops[m.o2];
    
    // Choose the sequence/base within the loops
    auto margin = m.margin / m.scale;
    auto c1 = p.seqs.choose_join(margin, m.index, false, state, state.model, rf); 
    NUPACK_REQUIRE(margin, <=, 1, c1, m.index, state);
    margin *= m.scale;
    auto c2 = k.seqs.choose_join(margin, m.index, true, state, state.model, rf); 
    NUPACK_REQUIRE(margin, <=, 1, c2, m.index, state);

    return ComplexJoinMove{
        .loop1=m.o1, .loop2=m.o2,
        .seq1=uint(c1.s), .seq2=uint(c2.s), .base1=c1.b, .base2=c2.b, 
        .dlogp=rf.log_mole_fraction + c1.dlogp_propensity->first + c2.dlogp_propensity->first,
        .rate=c1.dlogp_propensity->second * c2.dlogp_propensity->second * rf.bimolecular_propensity_prefactor()
    };
}
/******************************************************************************************/

/// Call a callback each join location in complex "x" of state "w" of base "b"
// that can pair to an external base "c"
template <class State, class X, class F>
void for_join_sites_in_complex(State &&w, X const &x, Base b, Base c, F &&f) {
    // for each exterior loop in the complex
    fork(w.rate_function, [&](auto const &rf) {
        for_exterior_loops_in_complex(w, x, [&](auto &&o) {
            // for each matching base in the exterior loop
            for_join_locs_in_loop(o.seqs, b, c, w.model, rf, [&](auto &&m) {
                // call back with the loop index and base location
                f(o.index(), m);
            });
        });
    });
}

/******************************************************************************************/

/// Call a callback on each ComplexJoin Move possible between complexes x and y in the given State
template <class State, class X, class F>
void for_joins_between(State const &w, X const &x, X const &y, F &&f) {
    // For every base pair (e.g. AT GC CG TA)
    for (auto b : w.model.alphabet().all()) for (auto c : w.model.alphabet().all()) if (w.model.pairing().can_close(b, c)) {
        // For every base in the first complex of nucleotide "b" which can pair to "c"
        for_join_sites_in_complex(w, x, b, c, [&](auto const &o1, auto const &m1){
            // For every base in the second complex of nucleotide "c" which can pair to "b"
            for_join_sites_in_complex(w, y, c, b, [&](auto const &o2, auto const &m2){
                // Call back with the join move between these bases
                f(ComplexJoinMove{.loop1=static_cast<iseq>(o1), .loop2=static_cast<iseq>(o2), .seq1=static_cast<uint>(m1.s), .seq2=static_cast<uint>(m2.s), .base1=m1.b, .base2=m2.b});
            });
        });
    }
}

/******************************************************************************************/

/// Call a callback on each ComplexJoin Move possible between complexes in the State
template <class State, class F>
void for_all_joins(State const &w, F &&f) {
    // iterate over all complexes (x, y) where x < y
    for_each_ordered_pair(begin_of(w.complexes), end_of(w.complexes),
        [&](auto const &x, auto const &y){
            // callback for each join move between the two complexes
             for_joins_between(w, *x, *y, f);
        }
    );
}

/******************************************************************************************/

template <class State, class O, class Predicate>
void QuadraticJoiner::update_join_rates(State const &state, Index i, O const &loop1, Predicate const &pred) {
    auto const &model = state.model;
    uint const strand = state.sys.strand_of(loop1.strand_begin());
    uint const complex_index = state.complexes.strand_map[strand].x;
    izip(state.complexes.strand_map, [&](auto is, auto const &s) {
        NUPACK_QUICK_REQUIRE(state.sys.strand_of(state.loops[s.loop].seqs.strand_begin()), ==, is, state.loops[s.loop]);

        if (is == strand) return;
        real rate = 0;
        if (s.x != complex_index) { // no joins between loops in same complex
            auto const &loop2 = state.loops[s.loop].seqs;
            auto const e0 = loop1.energy + loop2.energy;
            fork(state.rate_function, [&](auto const &rf) {
                for (auto s1 : iterators(loop1)) for (auto b1 : iterators(*s1).offset(+1, -1)) 
                    for (auto s2 : iterators(loop2)) for (auto b2 : iterators(*s2).offset(+1, -1)) 
                        if (pred(b1, b2) && model.pairing().can_close(*b1, *b2)) {
                            auto seqs1 = loop1, seqs2 = loop2;
                            auto [n1, n2] = seqs1.associate(seqs2, s1-loop1.begin(), s2-loop2.begin(), b1, b2, state, model, rf);
                            auto dlogp = -model.beta * (seqs1.energy + seqs2.energy + model.join_penalty() - e0) + rf.log_mole_fraction;
                            rate += rf.bimolecular_addition_weight(dlogp)
                                 * rf.loop_prefactor(true, seqs1, n1 ? len(seqs1) - n1 : 0)
                                 * rf.loop_prefactor(true, seqs2, n2 ? len(seqs2) - n2 : 0);
                        }
            });
        }
        join_rates.update(join_index(is, strand), rate);
    });
}

/******************************************************************************************/

template <class State, class RF, class Predicate>
ComplexJoinMove QuadraticJoiner::get_join_move(State const &state, real rate0, RF const &rf, Predicate const &pred) const {
    auto const &model = state.model;
    
    auto const [entry, margin] = join_rates.find(rate0);
    NUPACK_QUICK_REQUIRE(state.complexes.strand_map[strand_pairs[entry].first].x, 
                     !=, state.complexes.strand_map[strand_pairs[entry].second].x);
    auto const iloop1 = state.complexes.strand_map[strand_pairs[entry].first].loop;
    auto const iloop2 = state.complexes.strand_map[strand_pairs[entry].second].loop;
    auto const &loop1 = state.loops[iloop1].seqs;
    auto const &loop2 = state.loops[iloop2].seqs;

    auto const e0 = loop1.energy + loop2.energy;

    real rate_sum = 0;
    for (auto s1 : iterators(loop1)) for (auto b1 : iterators(*s1).offset(+1, -1))
        for (auto s2 : iterators(loop2)) for (auto b2 : iterators(*s2).offset(+1, -1)) 
            if (pred(b1, b2) && model.pairing().can_close(*b1, *b2)) {
                auto seqs1 = loop1, seqs2 = loop2;
                auto const [n1, n2] = seqs1.associate(seqs2, s1-loop1.begin(), s2-loop2.begin(), b1, b2, state, model, rf);
                auto dlogp = -model.beta * (seqs1.energy + seqs2.energy + model.join_penalty() - e0) + rf.log_mole_fraction;
                real rate = rf.bimolecular_addition_weight(dlogp)
                        * rf.loop_prefactor(true, seqs1, n1 ? len(seqs1) - n1 : 0)
                        * rf.loop_prefactor(true, seqs2, n2 ? len(seqs2) - n2 : 0);
                rate_sum += rate;

                if (rate_sum >= margin) {
                    return ComplexJoinMove{.loop1=iloop1, .loop2=iloop2, 
                        .seq1=uint(s1-loop1.begin()), .seq2=uint(s2-loop2.begin()), .base1=b1, .base2=b2, 
                        .dlogp=dlogp, .rate=rate};
                }
            }
    NUPACK_ERROR("Join move not found", rate_sum, margin, rate0, join_rates.total(), state, loop1, loop2);
}

/******************************************************************************************/

template <class State>
void QuadraticJoiner::initialize_joiner(State const &w) {
    strand_pairs.clear();
    std::size_t const n_strands = len(w.sys.strands());
    join_rates.reset(0, n_strands * (n_strands - 1) / 2);
    for (auto j : range(n_strands)) for (auto i : range(j)) strand_pairs.emplace_back(i, j);
        for (auto i : indices(strand_pairs)) NUPACK_REQUIRE(i, ==, join_index(strand_pairs[i].first, strand_pairs[i].second));
}

/******************************************************************************************/

template <class Model>
BasePairList base_pair_list(Model const &model) {
    BasePairList base_pairs;
    for (auto b : model.alphabet().all()) for (auto c : model.alphabet().all()) 
        if (model.pairing().can_close(b, c)) base_pairs.emplace_back(b, c);
    return base_pairs;
}

/******************************************************************************************/

template <class State>
void ProductJoiner::initialize_joiner(State const &w) {
    std::size_t const n_strands = len(w.sys.strands()), n_complexes = len(w.complexes), n_bases = w.model.alphabet().length();
    base_pairs = base_pair_list(w.model);

    std::size_t const n = length(w.model);

    complex_rates.reset(Matrix(n, 2, la::fill::zeros), n_complexes);
    join_rates = vec<Matrix>(n_strands, Matrix(n, 2, la::fill::zeros));
    complex_rates_squared.reset(Column(n, la::fill::zeros), n_complexes);
}


/******************************************************************************************/

//// Given two strands indexed at "p" and "k", split the complex they belong to in 2.
/// If k.pos > p.pos, make a new complex from [p:k)
/// If p.pos > k.pos, make a new complex from [k:p)
template <class State>
void ProductJoiner::register_split(State &w, Index p, Index k) {
    auto [newx, old_ix] = w.complexes.register_split(p, k);
    auto diff = complex_rates.zero_value;
    for (auto j : newx) diff += join_rates[j];
    complex_rates.emplace_back(diff); 
    complex_rates.update(old_ix, complex_rates[old_ix] - diff);

    complex_rates_squared.emplace_back(calculate_rate_product(complex_rates.back()));
    complex_rates_squared.update(old_ix, calculate_rate_product(complex_rates[old_ix]));
}

/******************************************************************************************/

template <class State>
void QuadraticJoiner::register_split(State &w, Index p, Index k) {
    auto [newx, old_ix] = w.complexes.register_split(p, k);

    for (auto strand : newx) // this will revise all new--old contact info. a bit extra work is done... but this is a slower method anyway
        update_join_rates(w, strand, w.loops[w.complexes.strand_map[strand].loop].seqs);
}

/******************************************************************************************/

}