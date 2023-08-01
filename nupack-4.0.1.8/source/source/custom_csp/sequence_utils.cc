#include <nupack/design/custom_csp/sequence_utils.h>
#include <nupack/design/custom_csp/pathway_utils.h>
#include <nupack/design/custom_csp/design_debug.h>
#include <nupack/design/custom_csp/nupack_invariants.h>

namespace nupack {
namespace custom_csp {
namespace SequenceUtils {

vec<trinary> nuc_to_bool(Alphabet const &a, Wildcard w) {
    return vmap<vec<trinary>>(a.all(), [w](Base b) {return w.includes(b);});
}

vec<vec<trinary>> nucs_to_bools(Alphabet const &a, Domain const &d) {
    return indirect_view(d, [&](Wildcard w) {return nuc_to_bool(a, w);});
}

Wildcard bool_to_nuc(vec<trinary> const &v) {return Wildcard::from_mask(v);}

vec<Wildcard> bools_to_nucs(const vec<vec<trinary>> &in) {
    return indirect_view(in, [](auto const &v) {return bool_to_nuc(v);});
}

}
}
}
