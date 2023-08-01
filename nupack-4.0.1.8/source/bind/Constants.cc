#include "Bind.h"
#include <nupack/execution/Local.h>
#include <nupack/types/IO.h>
#include <nupack/types/Sequence.h>
#include <nupack/types/Named.h>
#include <nupack/types/Structure.h>
#include <nupack/common/Costs.h>

namespace nupack {

/******************************************************************************************/

void render(Document &doc, Type<Local> t) {
    doc.type(t, "core.Local");
    doc.method<0>(t, "new", rebind::construct<usize>(t));
    doc.method(t, "n_workers", &Local::n_workers);
}

/******************************************************************************************/

void render_constants(Document &doc) {
    doc.render<Base>();
    doc.render<Wildcard>();

    doc.render<Sequence>();
    doc.render<Domain>();

    doc.render<NamedSequence>();
    doc.render<Complex>();
    doc.render<NamedComplex>();
    doc.render<Tube>();

    doc.render<NamedDomain>();
    doc.render<TargetStrand>();
    doc.render<TargetComplex>();
    doc.render<TargetTube>();

    doc.render<PairList>();
    doc.render<Structure>();

    doc.function("core.open_log", [](std::string s) {
        io::set_out(std::make_shared<std::ofstream>(s));
    });

    doc.function("core.close_log", [] {
        io::set_out(std::shared_ptr<std::ostream>(&std::cout, NoOp()));
    });

    doc.function("constants.ldexp", [](float f, std::int32_t i) {
        print(bool(i));
    });

    doc.function("constants.read_lines", [](std::string_view path) {
        std::ifstream ifs{path.data()};
        if (!ifs.good()) NUPACK_ERROR("File does not exist", path);
        vec<std::string> lines;
        for (std::string s; std::getline(ifs, s);) lines.emplace_back(s);
        std::cout << lines.front() << std::endl;
        return lines;
    });

    doc.function("core.test_matrix", [](Mat<double> const &x) {return x;});

    doc.function("constants.run_length_encoding", io::run_length_encoding);

    doc.render<PairMatrix<real>>();

    /// Utility functions
    doc.function("constants.dp_to_pairs", [](std::string_view s) {return io::to_pairs(s);});

    doc.function("constants.unit_evaluation_cost_table", unit_evaluation_cost_table);
    doc.function("constants.unit_evaluation_costs", unit_evaluation_costs);
    doc.function("constants.unit_subblock_cost", unit_subblock_cost);
    doc.function("constants.subblock_cost", subblock_cost<vec<std::size_t>>);

    doc.function("constants.trim_cxx", [](std::string s) {return trim_type_name(std::move(s), 10000);});
    doc.function("constants.rotational_symmetry", &rotational_symmetry<vec<uint>>);
    doc.function("core.compute_necklaces", [](rebind::AnnotatedCallback<void, vec<uint>> f, uint size, uint n) {
        return compute_necklaces(vec<uint>(size, 0), n, std::move(f));
    });

    // doc.function("numeric.nnls", nnls<float, std::uint16_t, float>);
    // doc.function("numeric.nnls", nnls<float, float, float>);
    // doc.function("numeric.nnls", nnls<double, double, double>);
    doc.function("constants.water_molarity", water_molarity);
    doc.function("constants.dna_salt_correction", dna_salt_correction);
    doc.object("constants.ZeroCinK", ZeroCinK);
    doc.object("constants.default_temperature", DefaultTemperature);
    doc.object("constants.boltzmann_constant", Kb);
    doc.object("constants.git_branch", GitBranch);
    doc.object("constants.git_revision", GitRevision);
    doc.object("constants.version", Version);
    doc.object("constants.build_type", BuildType);

#   define NUPACK_TMP(scope, name, val)         \
        doc.function(scope name, [] {return val;}); \
        doc.function(scope "set_" name, [](decltype(val) const &v) {val = v;})
        NUPACK_TMP("constants.", "default_parameters_path", DefaultParametersPath);
        NUPACK_TMP("constants.", "total_ram", TotalRAM);
        NUPACK_TMP("constants.", "total_cpu", TotalCPU);
#   undef NUPACK_TMP
}

void render(Document &doc, Type<AlwaysTrue> t) {doc.type(t, "constants.AlwaysTrue");}

void render(Document &doc, Type<AlwaysFalse> t) {doc.type(t, "constants.AlwaysFalse");}

void render(Document &doc, Type<True> t) {doc.type(t, "constants.TrueType");} // TODO convert to bool?

void render(Document &doc, Type<False> t) {doc.type(t, "constants.FalseType");}

}

/******************************************************************************************/

namespace rebind {

void Renderer<nupack::json>::operator()(Document &doc) const {
    using J = nupack::json;
    Type<J> t;
    doc.type(t, "core.JSON");
    doc.method(t, "new", construct<>(t));
    doc.method(t, "new", [](std::string_view s) {return J::parse(s);});
    doc.method(t, "load", [](J &j, std::string_view s) {j = J::parse(s);});
    doc.method(t, "dump", [](J const &j, unsigned indent) {return indent ? j.dump(indent) : j.dump();});

    doc.method(t, "load_binary", [](J &j, std::vector<std::uint8_t> const &s) {j = J::from_msgpack(s);});
    doc.method(t, "dump_binary", [](J const &j) {return J::to_msgpack(j);});

    doc.method(t, "load_file", [](J &j, std::string_view path) {
        std::ifstream ifs{std::string(path)};
        if (!ifs.good()) NUPACK_ERROR("invalid file", path);
        ifs >> j;
    });
    nupack::render_comparisons(doc, t);
}

}
