#pragma once
#include "../reflect/Reflection.h"
#include "../reflect/Serialize.h"

namespace nupack::design {

// Triages data to multiple loggers
struct Logs {
    // "time", "type", "depth", "psi_active", "psi_passive", "sequence", "defect"
    std::function<void(real, string_view, int, int, int, json const &, real)> basic;
    // "index", "name", "decomposition"
    std::function<void(uint, string_view, json const &)> decomposition;

    NUPACK_REFLECT(Logs, basic, decomposition);
};

// static_assert(std::is_convertible_v<Logs, json>);
// static_assert(!std::is_convertible_v<Logs, json>);
static_assert(std::is_constructible_v<json, Logs>);
// static_assert(!std::is_constructible_v<json, Logs>);


struct EngineObserver {
    // "type", "length", "time", "cache possible"
    std::function<void(string_view, uint, real, bool)> log;
    uint slowdown=0;

    NUPACK_REFLECT(EngineObserver, log, slowdown);
};

extern EngineObserver NullEngineObserver;

}
