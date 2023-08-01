#pragma once
#include "Models.h"
#include "Complex.h"
#include "Defect.h"
#include "Granularity.h"
#include "Logging.h"

namespace nupack::design {

using ComplexWeights = std::map<uint, vec<real>>;

struct sum_vec_t {
    template <class A, class B>
    auto operator()(A const &a, B const &b) const {
        if (len(a) == 0) return b;
        if (len(b) == 0) return a;

        NUPACK_REQUIRE(len(a), ==, len(b), "vectors must be same length if both not empty.");
        decay<decltype(a)> ret(len(a), 0.0);
        izip(a, b, [&](auto i, auto ai, auto bi) {
                at(ret, i) = ai + bi;
        });
        return ret;
    };
};

static constexpr auto sum_vec = sum_vec_t();


struct TubeTarget : MemberOrdered {
    uint complex_index;
    real target_conc;
    vec<uint> indices;
    NUPACK_REFLECT(TubeTarget, complex_index, target_conc, indices);

    /* Constructors */
    TubeTarget() = default;
    TubeTarget(uint c, real conc=0.0) : complex_index(c), target_conc(conc) {}

    bool is_on_target() const {return target_conc > 0.0;}
};


struct Tube : MemberOrdered {
    /** @brief includes indices to complexes included in the master designer list. */
    vec<TubeTarget> targets;
    string name;
    real_mat stoichiometry;
    real nucleotide_concentration;
    real water_conc;

    NUPACK_REFLECT(Tube, name, targets, water_conc, stoichiometry, nucleotide_concentration);

    Tube() = default;
    Tube(vec<TubeTarget> c, string name, vec<Complex> const &cs, real water) :
        targets(std::move(c)), name(std::move(name)), water_conc(water) {
        compute_invariants(cs);
    }

    /* initialization functions requiring const access to vector of complexes */
    /**
     * @brief compute the total concentration of individual nucleotides in the tube
     *
     * @param complexes the set of all complexes in the design
     * @return the total concentration of nucleotides
     */
    void compute_nucleotide_concentration(vec<Complex> const &cs) {
        nucleotide_concentration = sum(targets, [&](auto const &t) {return t.target_conc * nt(cs[t.complex_index]);});
    }
    void compute_stoichiometry(vec<Complex> const &cs);
    void store_complex_indices(vec<Complex> const &cs);
    void compute_invariants(vec<Complex> const &cs) {
        compute_stoichiometry(cs);
        compute_nucleotide_concentration(cs);
        store_complex_indices(cs);
    }
    vec<StrandView> strand_types(vec<Complex> const &cs) const;
    // real_mat single_strand_matrix(vec<Complex> const &cs) const;

    vec<real> concentrations(vec<real> const &log_pfuncs, EnsemblePartition const &part) const;
    vec<real> fractions(vec<real> const &log_pfuncs, EnsemblePartition const &part) const;

    vec<real> concentrations(vec<real> const &log_pfuncs) const;
    vec<real> fractions(vec<real> const &log_pfuncs) const;



    Defect defect(vec<real> const &log_pfuncs, vec<Defect> const &comp_defects,
        EnsemblePartition const &part={}, ComplexWeights const &weights={}) const;
    Defect normalized_defect(vec<real> const &log_pfuncs, vec<Defect> const &comp_defects,
        EnsemblePartition const &part={}, ComplexWeights const &weights={}) const;

    std::tuple<real_mat, real_col, real_col> deflate(vec<real> const &log_pfuncs, EnsemblePartition const &part) const;
    real_col reinflate(real_col const &x, EnsemblePartition const &part) const;

    /**
     * @brief create a view of the complexes making up the tube from the targets and vector of complexes
     *
     * @param cs the complexes which the tube is indexed into
     * @return view returning const references to each of the complexes in the tube
     */
    auto complexes(vec<Complex> const &cs) const {
        return indirect_view(targets, [&](auto const &t) -> Complex const & {return cs[t.complex_index];});
    }

    // uint add_missing_strands(real_mat &A, real_col &x0, real_col &dG, vec<Complex> const &cs) const;
};

void remove_added_strands(uint num_strands, real_col &x);

Defect structural_defect(TubeTarget const &t, Defect const &comp_defect, real concentration);
Defect concentration_defect(TubeTarget const &t, real concentration);

/* base, shared computation */
real_col solve_concentrations(real_mat const &A, real_col const &x0, real_col const &dG);

}
