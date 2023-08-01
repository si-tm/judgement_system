#pragma once
#include "Schedule.h"
#include "Action.h"
#include "Common.h"
#include "Cache.h"
#include "../model/Model.h"
#include "../standard/Map.h"
#include "../standard/Function.h"
#include "../types/PairList.h"
#include "../types/Sequence.h"
#include "../math/Sparse.h"

#include <atomic>

namespace nupack::thermo {

/******************************************************************************************/

struct StructureEnergy {
    PairList structure;
    real energy;
    real stack_energy;

    NUPACK_REFLECT(StructureEnergy, structure, energy, stack_energy);
};

struct Result {
    struct PF     {real logq, raw_logq; NUPACK_REFLECT(PF, logq, raw_logq);};
    struct MFE    {real energy, raw_energy; NUPACK_REFLECT(MFE, energy, raw_energy);};
    struct Pairs  {PairMatrix<real> matrix; NUPACK_REFLECT(Pairs, matrix);};
    struct Costs  {Mat<real> matrix; NUPACK_REFLECT(Costs, matrix);};
    struct Sample {vec<PairList> structures; NUPACK_REFLECT(Sample, structures);};
    struct Subopt {vec<StructureEnergy> structures; NUPACK_REFLECT(Subopt, structures);};

    std::optional<PF>     pfunc;
    std::optional<MFE>    mfe;
    std::optional<Pairs>  pairs;
    std::optional<Subopt> subopt;
    std::optional<Sample> sample;
    std::optional<Costs>  costs;

    NUPACK_REFLECT(Result, pfunc, mfe, pairs, subopt, sample, costs);
};

using Update = std::variant<Result::PF, Result::MFE, Result::Pairs, Result::Subopt, Result::Sample, Result::Costs>;

/******************************************************************************************/

struct Job {
    struct PF {
        NUPACK_REFLECT_EMPTY(PF);
    };
    struct MFE {
        NUPACK_REFLECT_EMPTY(MFE);
    };
    struct Pairs {
        Sparsity sparsity;
        NUPACK_REFLECT(Pairs, sparsity);
    };
    struct Subopt {
        real gap=0;
        std::size_t max_number=100000;
        NUPACK_REFLECT(Subopt, gap, max_number);
    };
    struct Costs {
        NUPACK_REFLECT_EMPTY(Costs);
    };
    struct Sample {
        std::size_t number;
        std::uint64_t seed; // The real seed used for the RNG initialization
        NUPACK_REFLECT(Sample, number, seed);
    };

    using Kind = std::variant<PF, Pairs, Sample, MFE, Costs, Subopt>;
    Complex strands; // also need indices on each strand.
    Kind kind;

    NUPACK_REFLECT(Job, strands, kind);
};

using Sink = std::function<void(Complex const &, Update)>;

/******************************************************************************************/

// Undefined class that should be specialized for each model type to yield a compatible thermo computer generator
template <class U>
struct DefaultFactory;

/******************************************************************************************/

using ComputerList = vec<Computer>;

struct Factory {
    ComputerList computers;

    Factory() = default;

    template <class Model>
    Factory(Model mod) : computers(DefaultFactory<Model>::create(std::move(mod))) {}
};

/******************************************************************************************/

// Specialization to make a factory of thermo computers using a Model<>
template <>
struct DefaultFactory<Model<>> {
    static ComputerList create(Model<>);
};

/******************************************************************************************/

using ResultMap = std::map<Complex, Result>;

/******************************************************************************************/

struct Future {
    std::future<void> future;
    SharedExecutor executor;
    std::shared_ptr<tf::Taskflow> flow;
    std::shared_ptr<ResultMap> result;
    SharedError error;

    NUPACK_REFLECT(Future, future, executor, flow, result, error);

    ResultMap get();
};

/******************************************************************************************/

void compute_callback(Sink put, vec<Job> jobs, Factory const &, ComputeOptions const &ops={});

Future submit(vec<Job> jobs, Factory const &, ComputeOptions const &ops={});

tf::Task schedule(tf::Taskflow &flow, Sink put, vec<Job> jobs, Factory const &, ComputeOptions const &ops);

ResultMap computes(vec<Job> jobs, Factory const &, ComputeOptions const &ops={});

Result compute(Job, Factory const &, ComputeOptions const &ops={});

std::pair<PairMatrix<real>, real> sparse_pair_probability(SequenceList seqs, Factory const &, ComputeOptions const &ops={});

std::pair<Mat<real>, real> pair_probability(SequenceList seqs, Factory const &, ComputeOptions const &ops={});

std::pair<Mat<real>, real> mfe_cost(SequenceList seqs, Factory const &, ComputeOptions const &ops={});

std::pair<vec<PairList>, real> sample(std::size_t n, SequenceList seqs, Factory const &, ComputeOptions const &ops={});

vec<StructureEnergy> subopt(real gap, SequenceList seqs, Factory const &, ComputeOptions const &ops={});

real log_partition_function(SequenceList seqs, Factory const &, ComputeOptions const &ops={});

real min_free_energy(SequenceList seqs, Factory const &, ComputeOptions const &ops={});

// Compute unpaired melt fraction for a list of complexes and strand mole fractions; assumes strand IDs have been set to count up from 0!
real fraction_bases_unpaired(vec<Complex> complexes, Factory const &factory, Col<real> const &mole_fractions, ComputeOptions const &ops={});

// Compute unpaired melt fractions for a list of distinguishable strands, temperatures, model generating function, strand concentrations, and max complex size
Col<real> melt_curve(vec<Sequence> strands, Col<real> const &concentrations, uint max_size, Col<real> const &temperatures, std::function<Factory(real)> factory, ComputeOptions const &ops={});

/******************************************************************************************/

}
