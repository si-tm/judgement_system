#include <nupack/proto/Design.h>
#include <nupack/proto/YAML.h>
#include "Design.h"

namespace nupack::design {

template <class S, class ...Ts>
struct SafeConstruct {
    static S call(Ts ...ts) {
        S s{std::forward<Ts>(ts)...};
        s.check_constraint();
        return s;
    }
};

template <class T, class ...Ts>
auto make_vector(T &&t, Ts &&...ts) {
    std::vector<std::decay_t<T>> out;
    out.reserve(1 + sizeof...(Ts));
    (out.emplace_back(std::forward<Ts>(ts)), ...);
    return out;
}

template <class ...Ts, class S>
auto safe_construct(Type<S>) {return &SafeConstruct<S, Ts...>::call;}

void render(Document &doc, Type<Env> t) {
    doc.type(t, "design.core.Env");
    doc.method(t, "new", [](SharedExecutor ex) {
        NUPACK_ASSERT(ex);
        return Env(Local(1), std::move(ex));
    });
}

void render(Document &doc, Type<Timer> t) {
    doc.type(t, "design.components.Timer");
    doc.method(t, "new", rebind::construct(t));
    doc.method(t, "start", &Timer::start);
    doc.method(t, "elapsed", &Timer::elapsed);
    doc.method(t, "stop", &Timer::stop);
}

void render(Document &doc, Type<DesignStats> t) {
    doc.type(t, "design.results.Stats");

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<EnsemblePartition> t) {
    doc.type(t, "design.results.Partition");

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}

void render(Document &doc, Type<SSMObjective> t) {
    doc.type(t, "design.constraints.SSM");
    doc.method(t, "new", rebind::construct<vec<TargetComplex>, uint, real>(t));
    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}

void render(Document &doc, Type<EnergyEqualizationObjective> t) {
    doc.type(t, "design.constraints.EnergyMatch");
    doc.method(t, "new", rebind::construct<vec<NamedDomain>, Optional<real>, real>(t));
    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}

void render(Document &doc, Type<SingleResult> t) {
    doc.type(t, "design.results.Single");

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<ComplexResult> t) {
    doc.type(t, "design.results.ComplexResult");

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<TubeComplex> t) {
    doc.type(t, "design.results.TubeComplex");

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<TubeResult> t) {
    doc.type(t, "design.results.TubeResult");

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<DesignResult> t) {
    doc.type(t, "design.results.RawResult");
    doc.method(t, "new", rebind::construct(t));
    render_json(doc, t);
    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}

void render(Document &doc, Type<CompareConstraint> t) {
    doc.type(t, "design.constraints.Compare");
    doc.method(t, "new", rebind::construct<vec<NamedDomain>, vec<NamedDomain>>(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}

void render(Document &doc, Type<MatchConstraint> t) {
    doc.type(t, "design.constraints.Match");
    doc.method(t, "new", rebind::construct<vec<NamedDomain>, vec<NamedDomain>>(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<ComplementarityConstraint> t) {
    doc.type(t, "design.constraints.Complementarity");
    doc.method(t, "new", rebind::construct<vec<NamedDomain>, vec<NamedDomain>, bool>(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}

void render(Document &doc, Type<PairingConstraint> t) {
    doc.type(t, "design.constraints.Pairing");
    doc.method(t, "new", rebind::construct<vec<NamedDomain>, vec<NamedDomain>>(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<PatternConstraint> t) {
    doc.type(t, "design.constraints.Pattern");
    doc.method(t, "new", safe_construct<vec<NamedDomain>, DomainList, real>(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}

void render(Document &doc, Type<DiversityConstraint> t) {
    doc.type(t, "design.constraints.Diversity");
    doc.method(t, "new", safe_construct<vec<NamedDomain>, uint, uint>(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}

void render(Document &doc, Type<LibraryConstraint> t) {
    doc.type(t, "design.constraints.Library");
    doc.method(t, "new", safe_construct<vec<NamedDomain>, vec<DomainList>>(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}

void render(Document &doc, Type<WindowConstraint> t) {
    doc.type(t, "design.constraints.Window");
    doc.method(t, "new", safe_construct<vec<NamedDomain>, DomainList>(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<SimilarityConstraint> t) {
    doc.type(t, "design.constraints.Similarity");
    doc.method(t, "new", safe_construct<vec<NamedDomain>, Domain, std::pair<real, real>, real>(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


using Condition = rebind::Callback<rebind::Integer>; // This has to be so for now, bindings are a bit broken if using a non-primitive type.
using Handler = rebind::Callback<void>;

struct DesignRunner {
    auto operator()(Specification const &spec, Env env, Condition condition, Handler handler, std::optional<DesignResult> restart_={}) const {
        double max_time = spec.options.max_time;

        SignalRuntime l;

        NUPACK_ASSERT(!handler.function || condition.function, "If using checkpointing with designer, you must supply a checkpoint condition");

        Designer d = spec.create_designer();

        // Use unlinked serial executor if the design seems pretty small
        real approx_nt = std::pow(sum(d.design.complexes, [](auto const &c) {
            return cube(sum(c.strands, [](auto const &s) {return real(s.size());}));
        }), 1.0/3.0);
        if (approx_nt < 1000) env = Env(1);

        d.initialize(env);

        if (bool(restart_)) {
            try {
                auto const &restart = restart_.value();
                auto const &res = restart.results[0];
                auto mapping = spec.ensure_compatibility(res); // throw useful errors where there are mismatches
                auto &seqs = d.design.sequences;

                for (auto const &domain: res.domains) seqs.set_domain(domain.first, domain.second);
                d.Psi.mask = vmap<decltype(d.Psi.mask)>(mapping, [&](auto i) {return restart.stats.final_Psi.active(i);});
                d.stats = restart.stats;
                d.redecompose_active(env, 0);
            } catch (...) {
                print("nupack: Failure in loading Design specification from intermediate result. Does it correspond to the correct design?");
                throw;
            }
        }
        auto end = std::chrono::steady_clock::now() + std::chrono::duration<double>(max_time);
        using TP = decltype(end);
        auto const deadline = max_time ? std::make_shared<TP>(end) : std::shared_ptr<TP>();
        if (condition.function) {
            // whether to emit a checkpoint
            auto real_condition = [c=std::move(condition), deadline](Designer const &des, bool done) {
                if (deadline && std::chrono::steady_clock::now() > *deadline) {return 2;}
                switch (c(des.stats, des.timer, done)) {
                    case +2: return 2; // finish the design.
                    case +1: return 1; // make a checkpoint
                    case -1: throw SignalError::sigint(); // throw KeyboardInterrupt
                    default: return 0; // don't make a checkpoint
                }
            };
            if (handler.function) {
                auto real_handler = [env, h=std::move(handler)](Designer &des) {
                    des.stats.design_time += des.timer.stop();
                    des.stats.final_Psi = des.Psi;
                    auto sequence = des.best_sequence(env);
                    auto result = DesignResult(env, des);
                    h(std::move(result));
                    des.timer.start();
                };
                d.checkpoint = [handler=std::move(real_handler), condition=std::move(real_condition)] (Designer &des, bool done) -> bool {
                    auto ret = condition(des, done);
                    if (ret == 1) handler(des);
                    return ret == 2;
                };
            } else {
                d.checkpoint = [condition=std::move(real_condition)] (Designer &des, bool done) -> bool {
                    return condition(des, done) == 2;
                };
            }
        } else {
            d.checkpoint = [deadline](Designer &des, bool done) -> bool {
                return deadline && std::chrono::steady_clock::now() > *deadline;
            };
        }

        d.optimize_tubes(env);

        // d.design.set_sequence(d.best.full.sequence);
        return DesignResult(env, d);
    }
};



void render(Document &doc, Type<Specification> t) {
    doc.type(t, "design.core.Specification");
    doc.method(t, "new", rebind::construct<vec<TargetTube>, Model<>, DesignParameters, vec<HardConstraint>, vec<SoftConstraint>, vec<Weight>, real>(t));
    render_json(doc, t);

    doc.method<5>(t, "run_one", DesignRunner());
    doc.method(t, "replace", [](Specification &spec, vec<NamedDomain> const &domains) {
        DomainMap map;
        for (auto const &d : domains) map.map.emplace(d.name, d);
        spec.replace(map);
    });
    doc.method(t, "evaluate", [](Specification const &spec, Env const &env) {
        for (auto const &t : spec.tubes)
            for (auto const &x : t.complexes)
                for (auto const &s : x.strands)
                    for (auto const &domain : s.domains)
            NUPACK_ASSERT(all_of(domain, is_determined), "Variable nucleotides in design evaluation", domain);

        Designer d = spec.create_designer();

        d.initialize(env);
        d.time_analysis(env);
        d.best.full = d.evaluate_objectives(env, 0, {}, d.weights);
        d.best.full.full_evaluation(env, d);

        return DesignResult(env, d);
    });

    doc.method(t, "to_yaml", [](Specification const &s) {
        return to_yaml_string(spec::to_proto<proto::DesignJob>(jobs::to_spec(s)));
    });

    doc.function("design.core.Specification.from_yaml", [](string const &s) {
        return jobs::from_spec(spec::from_proto(proto::from_yaml_string<proto::DesignJob>(s)));
    });

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<Weight> t) {
    doc.type(t, "design.weights.Weight");
    doc.method(t, "new", rebind::construct<string, string, string, string, real>(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<ReversedComplex> t) {
    doc.type(t, "design.weights.ReversedComplex");
    doc.method(t, "new", rebind::construct(t));
    doc.method(t, "domains", &ReversedComplex::domains);
    doc.method(t, "strands", &ReversedComplex::strands);

    doc.method(t, "{}", dumpable(t));
}


void render(Document &doc, Type<Weights> t) {
    doc.type(t, "design.weights.WeightFactors");
    doc.method(t, "new", rebind::construct(t));
    doc.method(t, "add", &Weights::add);
    doc.method(t, "add_objective_weight", &Weights::add_objective_weight);

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}

void render(Document &doc, Type<DesignParameters> t) {
    // using base_type = DesignParameters;
    doc.type(t, "design.core.Options");
    doc.method(t, "new", rebind::construct(t));
    doc.method(t, "new", rebind::construct<DesignParameters>(t));

    doc.method(t, "{}", dumpable(t));
    NUPACK_PUBLIC(t, seed, f_stop, f_passive, H_split, N_split, f_split,
            f_stringent, dG_clamp, M_bad, M_reseed, M_reopt, f_redecomp, f_refocus,
            cache_bytes_of_RAM, f_sparse, slowdown, log, decomposition_log, thermo_log,
            time_analysis, max_time, wobble_mutations);
}

}

namespace nupack {

void render_design(Document &doc) {
    doc.render<design::MultitubeObjective>();
    doc.render<design::TubeObjective>();
    doc.render<design::ComplexObjective>();
    doc.render<design::SSMObjective>();
    doc.render<design::SimilarityObjective>();
    doc.render<design::EnergyEqualizationObjective>();
    doc.render<design::PatternObjective>();
    doc.render<design::Objective>();

    doc.render<design::Specification>();
    doc.render<design::Timer>();

    doc.render<design::DesignResult>();

    doc.render<design::EnsemblePartition>();
}

}
