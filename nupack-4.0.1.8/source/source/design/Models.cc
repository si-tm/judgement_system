#include <nupack/design/Models.h>

namespace nupack::design {

thermo::ComputeOptions ThermoEnviron::compute_options(Env const &env) const {
    thermo::ComputeOptions ops;
    ops.executor = env.executor;
    ops.max_bytes = max_bytes;
    return ops;
}

ThermoEnviron & ModelMap::get(Model<real> const &key) const {
    //std::unique_lock lock(mut);
    return mod_map.try_emplace(key, key).first->second;
}


void ThermoEnviron::add_pfunc(::nupack::Complex const &s, real log_pfunc) {
    std::unique_lock lock(mut.mut);
    log_pfuncs.try_emplace(s, log_pfunc);
}


std::optional<real> ThermoEnviron::get_pfunc(::nupack::Complex const &s) const {
    std::shared_lock lock(mut.mut);
    if (auto it = log_pfuncs.find(s); it != end_of(log_pfuncs)) return it->second;
    else return {};
}


void ThermoEnviron::clear_cache() {for (auto &c : computers.computers) c.cache->clear();}

void ModelMap::clear_caches() {
    for (auto &i : item_view(mod_map)) i.clear_cache();
}


void ModelMap::create_caches(std::size_t ram) {
    auto num_models = len(mod_map);
    auto ram_per_model = ram / num_models;
    for (auto &i : item_view(mod_map)) i.max_bytes = ram_per_model;
}


}
