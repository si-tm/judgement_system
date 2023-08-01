#pragma once

#include "types.h"

#include "adapter.h"

#include <vector>
#include <string>

namespace nupack {
namespace custom_csp {

using namespace design;


struct NupackInvariants;

namespace SequenceUtils {

vec<trinary> nuc_to_bool(Alphabet const &, Wildcard in);
AllowTable nucs_to_bools(Alphabet const &, const vec<Wildcard> & in);
AllowTable nucs_to_bools(Alphabet const &, const Domain & in);

Wildcard bool_to_nuc(vec<trinary> const & in);
vec<Wildcard> bools_to_nucs(vec<vec<trinary> > const & in);


};

}
}
