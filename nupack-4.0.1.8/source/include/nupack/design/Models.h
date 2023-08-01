#pragma once

#include "TypeImports.h"
#include "../thermo/Engine.h"
#include "../model/Model.h"
#include "../common/Threading.h"

namespace nupack::design {

struct CopyableMutex {
    CopyableMutex() = default;
    CopyableMutex(CopyableMutex const &) : CopyableMutex() {}
    CopyableMutex(CopyableMutex &&) : CopyableMutex() {}
    CopyableMutex & operator=(CopyableMutex const &) {return *this;}
    CopyableMutex & operator=(CopyableMutex &&) {return *this;}

    std::shared_mutex mutable mut;
};

/**
 * @brief Wrapper for all arguments to thermo computation (model, cache, etc, cache of log pfs for shortcuts)
 *
 */
struct ThermoEnviron {
    Model<> model;
    thermo::Factory computers;
    std::unordered_map<::nupack::Complex, real> log_pfuncs;
    CopyableMutex mut;
    std::size_t max_bytes = 0;

    ThermoEnviron() = default;
    /* create cache of correct type immediately when environment is made */
    // ThermoEnviron(ModelsTuple mods) : models(std::move(mods)) {initialize_cache(0);};
    ThermoEnviron(Model<real> mod)
        : model(std::move(mod)), computers(model) {}

    thermo::ComputeOptions compute_options(Env const &) const;

    void add_pfunc(::nupack::Complex const &s, real log_pfunc);
    std::optional<real> get_pfunc(::nupack::Complex const &s) const;
    void clear_cache();

    NUPACK_REFLECT(ThermoEnviron, model, computers, log_pfuncs, max_bytes);
};


/**
 * @brief Maintains cache of all models needed during design to avoid constantly
 * recreating the model. Furthermore, keeps 32-bit and 64-bit versions of the
 * same model together to support model fallback operations in a relatively
 * seamless manner when calling thermo code.
 *
 */
class ModelMap : MemberOrdered {
    mutable std::map<Model<real>, ThermoEnviron> mod_map;

public:
    NUPACK_REFLECT(ModelMap, mod_map);

    ModelMap() = default;

    Model<> const & first_model() const {
        NUPACK_ASSERT(!mod_map.empty());
        return mod_map.begin()->first;
    }

    /**
     * @brief Create new model if key has not been seen before and return
     * pointer to model.
     *
     * @param key specification of the requested model
     * @return reference to the requested pair of models
     */
    ThermoEnviron & get(Model<real> const & key) const;
    // auto const &cached_models(Model<real> const & key) const {return get(key).models;}

    void create_caches(std::size_t ram);
    void clear_caches();
    auto size() const {return len(mod_map);}
};



}
