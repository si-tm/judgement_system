#pragma once

#include "TypeImports.h"
#include "Logging.h"
#include "Split.h"
#include "../math/Sparse.h"

namespace nupack::design {

struct EngineObserver;
struct ThermoEnviron;

real log_partition_function(Env const &env, ::nupack::Complex const &, ThermoEnviron &, EngineObserver &obs=NullEngineObserver);

std::pair<PairMatrix<real>, real> raw_pair_probability(Env const &env, ::nupack::Complex const &, ThermoEnviron &, Sparsity s, EngineObserver &obs=NullEngineObserver);

std::pair<PairMatrix<real>, real> bonus_pair_probability(Env const &env, ::nupack::Complex, ThermoEnviron &, Sparsity s, vec<SplitPoint> const &fixed_pairs, real bonus, EngineObserver &obs=NullEngineObserver);

}
