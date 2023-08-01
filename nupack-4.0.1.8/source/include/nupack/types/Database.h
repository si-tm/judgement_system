#pragma once
#include "../reflect/Serialize.h"
#include "Matrix.h"

namespace nupack {

/**
 * @brief Database of key-value pairs serializable to/from JSON - generally for testing
 *
 * @tparam V container type for the data (probably map or vector)
 */
template <class V>
struct JsonDatabase : ConstIterable<JsonDatabase<V>> {
    using const_iterator = const_iterator_of<V>;

    JsonDatabase() = default;

    JsonDatabase(string const &p) {
        auto s = path_join(DefaultDataPath, p);
        NUPACK_ASSERT(path_exists(s), "Failed to find data", DefaultDataPath, p);
        std::ifstream ifs(s);
        json j;
        ifs >> j;
        *this = std::move(j).get<JsonDatabase>();
    }

    auto const & iter() const {return data;}

    string modified, revision;
    V data;
    NUPACK_REFLECT(JsonDatabase, modified, revision, data);

    auto const & operator[](usize s) const {return *next(data, s);}

    void write(string const &p, bool prepend=false, uint indent=0) {
        modified = timestamp();
        revision = GitRevision;
        std::ofstream ofs(prepend ? path_join(DefaultDataPath, p) : p);
        if (indent) ofs << std::setw(indent);
        ofs << json{*this};
    };
};

template <bool variable_params=false>
struct EnergyDatum {
    string structure;
    real energy;
    NUPACK_REFLECT(EnergyDatum, structure, energy);
};

template <>
struct EnergyDatum<true> {
    string structure;
    real energy;
    string material;
    real T;
    real sodium;
    real magnesium;
    NUPACK_REFLECT(EnergyDatum, structure, energy, material, T, sodium, magnesium);
};

using EnergyDatabase = JsonDatabase<std::map<string, EnergyDatum<false>>>;

struct ConcentrationDatum {
    Col<real> x0; Col<real> x; Mat<real> A; Col<real> g;
    NUPACK_REFLECT(ConcentrationDatum, x0, x, A, g);
};

struct PairsDatum {
    vec<int> i; vec<int> j; vec<real> prob;
    NUPACK_REFLECT(PairsDatum, i, j, prob);
};

struct PfuncDatum {
    real result;
    string material;
    real T;
    real sodium;
    real magnesium;
    NUPACK_REFLECT(PfuncDatum, result, material, T, sodium, magnesium);
};

using PairsDatabase = JsonDatabase<std::map<string, PairsDatum>>;

}
