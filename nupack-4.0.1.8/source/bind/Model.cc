#include "Model.h"
#include <nupack/state/State.h>

namespace nupack {

/******************************************************************************************/

void render(Document &doc, Type<ParameterFile> t) {
    doc.type(t, "model.ParameterFile");
    doc.method(t, "new", rebind::construct<string>(t));
    NUPACK_PUBLIC(t, path);
}

void render(Document &doc, Type<ModelConditions> t) {
    doc.type(t, "model.Conditions");
    doc.method(t, "new", rebind::construct<>(t));
    render_public(doc, t);
}

void render(Document &doc, Type<ParameterInfo> t) {
    doc.type(t, "model.ParameterInfo");
    render_public(doc, t);
    doc.method(t, "new", rebind::construct<ParameterFile, string, real, real>(t));
}

/******************************************************************************************/

void render(Document &doc, Type<ParameterBase> t) {
    doc.type(t, "model.ParameterBase");
    render_public(doc, t);

    doc.method(t, "load_json", [](ParameterIndex i, Alphabet const &a, json const &j) {
        return load_parameter_data(a, i, j);
    });

    doc.method(t, "save_json", [](ParameterIndex const &i, Alphabet const &a, ParameterArray<real> const &p) {
        return save_parameter_data(a, i, p);
    });
}

/******************************************************************************************/

void render(Document &doc, Type<ParameterArray<real32>> t) {
    render(doc, t, 0);
}

void render(Document &doc, Type<ParameterArray<real64>> t) {
    render(doc, t, 0);
//     doc.method(t, "check", [](ParameterData<real64> const &x, string const &s) {
//         auto y = set_parameters(x, json::parse(s)["dG"]);
//         zip(x.as_array(), y.as_array(), [](auto x, auto y) {
//             NUPACK_REQUIRE(x, ==, y);
//         });
//     });
}

void render(Document &doc, Type<ParameterSet<real32>> t) {render(doc, t, 0);}
void render(Document &doc, Type<ParameterSet<real64>> t) {render(doc, t, 0);}
void render(Document &doc, Type<Model<real64>> t) {render(doc, t, 0);}
void render(Document &doc, Type<Model<real32>> t) {render(doc, t, 0);}

void render_model(Document &doc) {
    doc.render<ParameterSet<real64>>();
    doc.render<Model<real64>>();
    doc.render<ParameterSet<real32>>();
    doc.render<Model<real32>>();

    doc.function("model.loop_structure", find_loop_structure_nick);
}

/******************************************************************************************/

}
