#pragma once
#include "../Forward.h"
#include "../types/Complex.h"
#include "../model/Model.h"
#include "../execution/Local.h"
#include "../execution/Executor.h"

namespace nupack::design {


struct Env {
    Local pool;
    SharedExecutor executor;

    explicit Env(uint n) : pool(n), executor(n) {}
    Env(Local l, SharedExecutor e) : pool(std::move(l)), executor(std::move(e)) {}

    NUPACK_REFLECT(Env, pool, executor);
};


using ProbabilityMatrix = arma::sp_mat;

}
