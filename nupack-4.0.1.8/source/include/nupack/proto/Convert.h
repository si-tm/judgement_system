#pragma once
#include "Spec.h"
#include "../model/Model.h"
#include "../math/Sparse.h"
#include "../types/Structure.h"
#include "../types/Named.h"
#include "../design/Specification.h"

namespace nupack::jobs {

/******************************************************************************/

Ensemble from_spec(spec::Ensemble e);
ParameterFile from_spec(spec::Parameters p, string const &);

/******************************************************************************/

extern bool UseRestrictedParameters;
extern bool UseModelCache;

Model<> from_spec(spec::Model const &m);
spec::Model to_spec(Model<> const &m);

/******************************************************************************/

Complex from_spec(vec<spec::Strand> const &strands, spec::Complex const &x);
vec<Complex> from_spec(spec::ComplexSet const &job);

/******************************************************************************/

spec::PairProbability to_spec(la::SpMat<double> const &P);
spec::PairProbability to_spec(PairMatrix<double> const &);

PairMatrix<double> pairs_from_spec(spec::PairProbability const &);

/******************************************************************************/


struct UtilityPairsInput {
    Model<> model;
    Complex complex;
    Structure structure;
    NUPACK_REFLECT(UtilityPairsInput, model, complex, structure);
};

struct AnalysisInput {
    vec<Model<>> models;
    vec<Tube> tubes;
    NUPACK_REFLECT(AnalysisInput, models, tubes);
};

struct DesignInput {
    design::Specification spec;
    uint trials;
    NUPACK_REFLECT(DesignInput, spec, trials);
};

using WebInput = std::variant<UtilityPairsInput, AnalysisInput, DesignInput>;

AnalysisInput load_input(spec::FrontendAnalysisJob const &j);
DesignInput load_input(spec::FrontendDesignJob const &j);
UtilityPairsInput load_input(spec::FrontendUtilityJob const &j);

/******************************************************************************/

struct UtilityOutput {
    PairMatrix<real> pair_probabilities;
    Structure structure;
    real log_partition_function, structure_energy;
    NUPACK_REFLECT(UtilityOutput, pair_probabilities, structure, log_partition_function, structure_energy);
};

struct ComplexOutput {
    PairMatrix<real> pair_probabilities;
    Structure mfe_structure;
    real mfe_energy, log_partition_function;
    NUPACK_REFLECT(ComplexOutput, pair_probabilities, mfe_structure, mfe_energy, log_partition_function);
};

struct AnalysisTubeOutput {
    Col<real> concentrations;
    NUPACK_REFLECT(AnalysisTubeOutput, concentrations);
};

struct AnalysisOutput {
    vec<AnalysisTubeOutput> tubes;
    vec<ComplexOutput> complexes;
    NUPACK_REFLECT(AnalysisOutput, tubes, complexes);
};

struct DesignOutput {
    vec<design::DesignResult> trials;
    NUPACK_REFLECT(DesignOutput, trials);
};

using WebOutput = std::variant<UtilityOutput, AnalysisOutput, DesignOutput>;

AnalysisOutput load_output(spec::FrontendAnalysisResults const &r);
DesignOutput load_output(spec::FrontendDesignResults const &r);
UtilityOutput load_output(spec::FrontendUtilityResult const &r);

/******************************************************************************/

}