#include "Bind.h"
#include <nupack/types/Matrix.h>
#include <nupack/concentration/Equilibrate.h>

namespace nupack {

namespace concentration {

/**************************************************************************************/

void render(Document &doc, Type<Options> t);

/**************************************************************************************/

template <class V>
void render(Document &doc, Type<Output<V>> t) {
    doc.type(t, "concentration.Output");
    NUPACK_PUBLIC(t, solution, dual_solution, objective, error, iters, converged);
}

/**************************************************************************************/

void render(Document &doc, Type<Options> t) {
    doc.type(t, "concentration.Options");
    doc.method(t, "new", rebind::construct(t));
    doc.method(t, "set_method", [](Options &o, uint n) {NUPACK_REQUIRE(n, <, 6); o.method = static_cast<Method>(n);});
    NUPACK_PUBLIC(t, max_iters, tolerance, delta_min, delta_max, orthogonalize);
}

}

/**************************************************************************************/

void render_math(Document &doc) {
    doc.function("concentration.solve", &concentration::equilibrate);
    doc.function("concentration.solve_complexes", &concentration::solve_complexes);
}

/**************************************************************************************/

}
