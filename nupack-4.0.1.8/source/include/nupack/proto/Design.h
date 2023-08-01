#pragma once
#include "Spec.h"
#include "nupack/design/DesignParameters.h"
#include <nupack/design/Specification.h>

namespace nupack::jobs {

design::Specification from_spec(spec::DesignJob const &in);
design::DesignParameters from_spec(spec::DesignOptions const &in);
design::DesignStats from_spec(spec::DesignStats const &d);

/******************************************************************************/

spec::DesignStats to_spec(design::DesignStats const &d);
spec::DesignResult to_spec(design::SingleResult const &r, design::DesignStats const &stats, spec::DefectHistory history, Model<> const &model, std::uint32_t trial);


design::DesignResult from_spec(spec::DesignResult const &r);
//  {
//     design::DesignResult out;
//     // out.model = 
//     // out.parameters = 
//     // out.stats = 
//     // out.objectives = 
//     // out.results = 
//     // out.weights = 
//     return out;
// }

/******************************************************************************/

spec::DesignJob to_spec(design::Specification const &d);

// spec::Domain to_spec(design::DomainSpec const &);
// spec::Strand to_spec(design::StrandSpec const &);
// spec::Complex to_spec(design::ComplexSpec const &);
// spec::Tube to_spec(design::TubeSpec const &);
// spec::DesignOptions to_spec(design::Options const &);

// spec::DesignJob to_spec(design::Specification const &spec) {
//     return {
//         .model=to_spec(spec.model),
//         .domains=to_spec(spec.domains),
//         .strands=to_spec(spec.strands),
//         .complexes=to_spec(spec.complexes),
//         .tubes=to_spec(spec.tubes),
//         .constraints={},
//         .options=to_spec(spec.options),
//         .weights=to_spec(spec.weights)
//     };
// }

/******************************************************************************/

}