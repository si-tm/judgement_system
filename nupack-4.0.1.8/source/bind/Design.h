#pragma once
#include "Bind.h"
#include <nupack/design/Objectives.h>
#include <nupack/design/OutputResult.h>
#include <nupack/design/Specification.h>
#include <nupack/design/Design.h>
#include <nupack/design/DesignComponents.h>

namespace nupack::design {

/******************************************************************************/

void render(Document &doc, Type<Env> t);
void render(Document &doc, Type<ComplexResult> t);
// void render(Document &doc, Type<ComplexSpec> t);
void render(Document &doc, Type<ComplementarityConstraint> t);
void render(Document &doc, Type<PairingConstraint> t);
void render(Document &doc, Type<MatchConstraint> t);
void render(Document &doc, Type<PatternConstraint> t);
void render(Document &doc, Type<DiversityConstraint> t);
void render(Document &doc, Type<WindowConstraint> t);
void render(Document &doc, Type<LibraryConstraint> t);
void render(Document &doc, Type<SimilarityConstraint> t);

void render(Document &doc, Type<DesignParameters> t);
void render(Document &doc, Type<DesignResult> t);
void render(Document &doc, Type<DesignStats> t);

void render(Document &doc, Type<SingleResult> t);
void render(Document &doc, Type<Specification> t);

void render(Document &doc, Type<Timer> t);
void render(Document &doc, Type<TubeComplex> t);
void render(Document &doc, Type<TubeResult> t);

void render(Document &doc, Type<Weight> t);
void render(Document &doc, Type<Weights> t);

/******************************************************************************/

}
