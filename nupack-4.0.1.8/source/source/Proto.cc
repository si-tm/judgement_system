#include <nupack/proto/Convert.h>
#include <nupack/types/Named.h>
#include <nupack/proto/Spec.h>
#include <nupack/proto/Design.h>

namespace nupack::jobs {

/******************************************************************************/

Ensemble from_spec(spec::Ensemble e) {
    switch (e) {
        case spec::Ensemble::NOSTACKING: return Ensemble::nostacking;
        case spec::Ensemble::STACKING: return Ensemble::stacking;
        case spec::Ensemble::DANGLE: return Ensemble::dangle;
        case spec::Ensemble::COAXIAL: return Ensemble::coaxial;
        case spec::Ensemble::NONE: return Ensemble::none;
        case spec::Ensemble::SOME: return Ensemble::min;
        case spec::Ensemble::ALL: return Ensemble::all;
    }
    NUPACK_ERROR("Bad Ensemble value");
}

spec::Ensemble to_spec(Ensemble e) {
    switch (e) {
        case Ensemble::nostacking: return spec::Ensemble::NOSTACKING;
        case Ensemble::stacking: return spec::Ensemble::STACKING;
        case Ensemble::dangle: return spec::Ensemble::DANGLE;
        case Ensemble::coaxial: return spec::Ensemble::COAXIAL;
        case Ensemble::none: return spec::Ensemble::NONE;
        case Ensemble::min: return spec::Ensemble::SOME;
        case Ensemble::all: return spec::Ensemble::ALL;
    }
    NUPACK_ERROR("Bad Ensemble value");
}

bool UseRestrictedParameters = false;

ParameterFile from_spec(spec::Parameters p, string const &filename) {
    if (!UseRestrictedParameters && !filename.empty()) return ParameterFile(filename);
    switch (p) {
        case spec::Parameters::RNA: {return ParameterFile("rna06");}
        case spec::Parameters::DNA: {return ParameterFile("dna04");}
        case spec::Parameters::RNA95: {return ParameterFile("rna95");}
        case spec::Parameters::RNA95_NUPACK3: {return ParameterFile("rna95-nupack3");}
        case spec::Parameters::RNA99_NUPACK3: {return ParameterFile("rna99-nupack3");}
        case spec::Parameters::DNA04_NUPACK3: {return ParameterFile("dna04-nupack3");}
    }
    NUPACK_ERROR("Bad parameters value", static_cast<int>(p));
}

std::pair<string, spec::Parameters> get_parameter_name(string_view n) {
    string s(n);
    for (auto &c : s) c = std::toupper(c);
    std::replace(s.begin(), s.end(), '-', '_');
    if (s == "RNA95_NUPACK3")       return {{}, spec::Parameters::RNA95_NUPACK3};
    if (s == "RNA99_NUPACK3")       return {{}, spec::Parameters::RNA99_NUPACK3};
    if (s == "DNA04_NUPACK3")       return {{}, spec::Parameters::DNA04_NUPACK3};
    if (s == "RNA95")               return {{}, spec::Parameters::RNA95};
    if (s == "DNA" || s == "DNA04") return {{}, spec::Parameters::DNA};
    if (s == "RNA" || s == "RNA06") return {{}, spec::Parameters::RNA};
    return {string(n), spec::Parameters::RNA};
}

/******************************************************************************/

bool UseModelCache = false;

Model<> from_spec(spec::Model const &m) {
    ModelConditions const c{.temperature=m.kelvin, .na_molarity=m.sodium, .mg_molarity=m.magnesium};
    ParameterFile const p = from_spec(m.parameters, m.filename);
    Ensemble const e = from_spec(m.ensemble);

    if (!UseModelCache) return Model(e, p, c);
    static std::list<Model<>> cache;
    static std::mutex mut;
    {
        std::lock_guard<std::mutex> lock(mut);
        auto it = std::find_if(cache.begin(), cache.end(), [&](auto const &m) {
            return m.conditions == c && m.ensemble == e && m.parameters.info.file == p;
        });
        if (it == cache.end()) {
            cache.emplace_front(e, p, c);
            if (cache.size() > 10) cache.pop_back();
        } else {
            cache.splice(cache.begin(), cache, it);
        }
        return cache.front();
    }
}

spec::Model to_spec(Model<> const &m) {
    auto [s, p] = get_parameter_name(m.parameters.info.file.path);
    return {
        .kelvin=m.conditions.temperature,
        .sodium=m.conditions.na_molarity,
        .magnesium=m.conditions.mg_molarity,
        .parameters=p,
        .ensemble=to_spec(m.ensemble),
        .filename=std::move(s)
    };
}

/******************************************************************************/

Complex from_spec(std::vector<spec::Strand> const &strands, spec::Complex const &x) {
    std::unordered_map<std::size_t, Sequence> lookup;
    for (auto const &s : strands)
        lookup.emplace(s.index, DNA.sequence(s.sequence));

    return vmap<Complex>(x.strand_index, [&](auto i) {
        return lookup.at(i);
    });
}

vec<Complex> from_spec(spec::ComplexSet const &set) {
    std::unordered_map<std::size_t, Sequence> strands;
    for (auto const &s : set.strands)
        strands.emplace(s.index, DNA.sequence(s.sequence));

    return vmap(set.complex_indices, [&](auto const &x) {
        return vmap<Complex>(x.strand_index, [&](auto i) {
            return strands.at(i);
        });
    });
}

/******************************************************************************/

spec::PairProbability to_spec(PairMatrix<double> const &pairs) {
    return {
        .diagonal=view(pairs.diagonal),
        .rows=view(pairs.rows),
        .cols=view(pairs.cols),
        .values=view(pairs.values),
    };
}

spec::PairProbability to_spec(la::SpMat<double> const &P) {
    la::Col<double> diag(P.diag());
    return {
        .diagonal=view(diag),
        .rows=indirect_view(iterators(P), [](auto const &it) {return it.row();}),
        .cols=indirect_view(iterators(P), [](auto const &it) {return it.col();}),
        .values=view(P)
    };
}


PairMatrix<double> pairs_from_spec(spec::PairProbability const &P) {
    return PairMatrix<double>(
        vmap<Col<double>>(P.diagonal), 
        vmap<Col<double>>(P.values), 
        vmap<Col<std::uint32_t>>(P.rows), 
        vmap<Col<std::uint32_t>>(P.cols));
}

/******************************************************************************/

AnalysisInput load_input(spec::FrontendAnalysisJob const &j) {
    return {
        .models=indirect_view(j.models, [](auto const &m) {return from_spec(m);}),
        .tubes=indirect_view(j.tubes, [&](auto const &t) {
            vec<NamedSequence> strands = indirect_view(t.strands, [](auto const &s) {return NamedSequence(DNA.sequence(s.sequence), s.name, DNA);});
            return Tube(
               ComplexSet(
                   strands, 
                   indirect_view(t.complexes, [&](auto const &x) {return NamedComplex(indirect_view(x.strand_index, 
                       [&](auto i) {return strands.at(i);}), {});})
               ),
               indirect_view(t.log_concentrations, exp), 
               t.name
            );
        })
    };
}

DesignInput load_input(spec::FrontendDesignJob const &j) {
    return {.spec=from_spec(j.spec), .trials=j.trials};
}

UtilityPairsInput load_input(spec::FrontendUtilityPairsJob const &j) {
    return {
        .model=from_spec(j.model),
        .complex=indirect_view(j.complex.strand_index, [&](auto i) {return DNA.sequence(j.strands.at(i).sequence);}),
        .structure=Structure(j.structure)
    };
}

UtilityPairsInput load_input(spec::FrontendUtilityDesignJob const &j) {
    NUPACK_ERROR("Utility design jobs cannot be loaded");
}

UtilityPairsInput load_input(spec::FrontendUtilityJob const &j) {
    return fork(j, [](auto const &j) -> UtilityPairsInput {return load_input(j);});
}

AnalysisInput load_input(spec::FrontendConcentrationJob const &) {NUPACK_ERROR("Concentration jobs cannot be loaded");}

WebInput load_input(std::string_view fn) {
    std::ifstream ifs{std::string(fn)};
    auto o = spec::load_stream<proto::FrontendJob>(ifs);
    return fork(o, [](auto const &o) -> WebInput {return load_input(o);});
}

/******************************************************************************/

ComplexOutput load_complex_output(spec::PairsResult const &r, spec::MFEResult const &m) {
    return {
        .pair_probabilities=pairs_from_spec(r.pairs),
        .mfe_structure=Structure(m.mfe_structure),
        .mfe_energy=m.mfe,
        .log_partition_function=r.log_partition_function
    };
}

AnalysisOutput load_output(spec::FrontendAnalysisResults const &r) {
    AnalysisOutput out;
    std::map<Complex, ComplexOutput> complexes;
    for (auto const &t : r.tube_result) {
        out.tubes.emplace_back();
        zip(t.complexes, t.pairs.results, t.mfe.results, [&](auto const &c, auto const &p, auto const &m) {
            NUPACK_ERROR("not impl");
            auto [it, inserted] = complexes.try_emplace(Complex());
            if (inserted) it->second = load_complex_output(p, m);
        });
    }
    return out;
}

DesignOutput load_output(spec::FrontendDesignResults const &r) {
    return {vmap(r.results, [](auto const &r) {return from_spec(r);})};
}

AnalysisOutput load_output(spec::ExceptionResult const &) {
    NUPACK_ERROR("Exception results cannot be loaded");
}

AnalysisOutput load_output(spec::FrontendCancellationResult const &) {
    NUPACK_ERROR("Cancelled results cannot be loaded");
}

UtilityOutput load_output(spec::FrontendUtilityResult const &r) {
    return {
        .pair_probabilities=pairs_from_spec(r.pairs.pairs),
        .structure=Structure(r.structure),
        .log_partition_function=r.pairs.log_partition_function,
        .structure_energy=r.energy
    };
}

WebOutput load_output(std::string_view fn) {
    std::ifstream ifs{std::string(fn)};
    auto o = spec::load_stream<proto::FrontendResult>(ifs);
    return fork(o, [](auto const &o) -> WebOutput {return load_output(o);});
}

/******************************************************************************/

}