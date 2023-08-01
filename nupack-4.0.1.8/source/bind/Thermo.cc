#include "Thermo.h"
// #include <nupack/thermo/CachedModel.h>
// #include <nupack/thermo/Adapters.h>
#include <nupack/Forward.h>
#include <nupack/model/Model.h>
#include <nupack/thermo/Engine.h>


namespace nupack {

namespace thermo {

void render(Document &doc, Type<ComputeOptions> t) {
    doc.type(t, "thermo.ComputeOptions");
    doc.method(t, "new", [](SharedExecutor exec, std::size_t memory) {
        ComputeOptions ops;
        ops.max_bytes = memory;
        ops.executor = std::move(exec);
        return ops;
    });
    NUPACK_PUBLIC(t, max_bytes);
}

void render(Document &doc, Type<Cache> t) {doc.type(t, "thermo.Cache"); render_public(doc, t);}

void render(Document &doc, Type<StructureEnergy> t) {
    doc.type(t, "thermo.PairsEnergy");
    render_public(doc, t);
}

void render(Document &doc, Type<Job> t) {
    doc.type(t, "thermo.Job");
    doc.method(t, "new", rebind::construct<Complex, Job::Kind>(t));
    render_public(doc, t);
}

void render(Document &doc, Type<Job::PF> t) {
    doc.type(t, "thermo.PFJob");
    doc.method(t, "new", rebind::construct(t));
    render_public(doc, t);
}

void render(Document &doc, Type<Job::MFE> t) {
    doc.type(t, "thermo.MFEJob");
    doc.method(t, "new", rebind::construct(t));
    render_public(doc, t);
}

void render(Document &doc, Type<Job::Pairs> t) {
    doc.type(t, "thermo.PairsJob");
    doc.method(t, "new", rebind::construct<Sparsity>(t));
    render_public(doc, t);
}

void render(Document &doc, Type<Job::Costs> t) {
    doc.type(t, "thermo.CostsJob");
    doc.method(t, "new", rebind::construct(t));
    render_public(doc, t);
}

void render(Document &doc, Type<Job::Sample> t) {
    doc.type(t, "thermo.SampleJob");
    doc.method(t, "new", rebind::construct<std::size_t, std::uint64_t>(t));
    render_public(doc, t);
}

void render(Document &doc, Type<Job::Subopt> t) {
    doc.type(t, "thermo.SuboptJob");
    doc.method(t, "new", rebind::construct<real, std::size_t>(t));
    render_public(doc, t);
}

void render(Document &doc, Type<Result> t) {
    doc.type(t, "thermo.Result");
    render_public(doc, t);
}

void render(Document &doc, Type<Result::PF> t) {
    doc.type(t, "thermo.PFResult");
    render_public(doc, t);
}

void render(Document &doc, Type<Result::MFE> t) {
    doc.type(t, "thermo.MFEResult");
    render_public(doc, t);
}

void render(Document &doc, Type<Result::Pairs> t) {
    doc.type(t, "thermo.PairsResult");
    render_public(doc, t);
}

void render(Document &doc, Type<Result::Costs> t) {
    doc.type(t, "thermo.CostsResult");
    render_public(doc, t);
}

void render(Document &doc, Type<Result::Sample> t) {
    doc.type(t, "thermo.SampleResult");
    render_public(doc, t);
}

void render(Document &doc, Type<Result::Subopt> t) {
    doc.type(t, "thermo.SuboptResult");
    render_public(doc, t);
}

void render(Document &doc, Type<Future> t) {
    doc.type(t, "thermo.Future");
    doc.method(t, "get", [](Future &f) {return f.get();});
}

}

// void render(Document &doc, Type<ComplexSampler> t) {
//     doc.type(t, "thermo.ComplexSampler");
//     doc.method(t, "new", rebind::construct<SequenceList const &, vec<std::pair<Complex, real>> const &, real>(t));
//     doc.method(t, "()", [](ComplexSampler &s, Local &env, CachedModel<PF, Model<>> &mod, usize n) {return s(env, mod, n);});
// }

void render_thermo(Document &doc) {
    doc.function("thermo.submit", [](vec<thermo::Job> jobs, Model<> const &model, thermo::ComputeOptions const &ops) {
        // for (auto &j : jobs) for (auto &s : j.strands) s.id = 0; // don't distinguish strands.
        return submit(std::move(jobs), model, ops);
    });
}

/******************************************************************************************/

}


