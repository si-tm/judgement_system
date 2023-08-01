#pragma once
#include "Bind.h"
#include <nupack/thermo/Engine.h>
#include <nupack/thermo/CachedModel.h>
#include <nupack/thermo/ComplexSampler.h>
#include <nupack/types/Structure.h>
#include <nupack/execution/Local.h>

namespace nupack::thermo {

void render(Document &, Type<ComputeOptions>);

void render(Document &, Type<Cache>);

void render(Document &, Type<Job>);
void render(Document &, Type<Job::PF>);
void render(Document &, Type<Job::MFE>);
void render(Document &, Type<Job::Pairs>);
void render(Document &, Type<Job::Costs>);
void render(Document &, Type<Job::Sample>);
void render(Document &, Type<Job::Subopt>);

void render(Document &, Type<StructureEnergy>);

void render(Document &, Type<Result>);
void render(Document &, Type<Result::PF>);
void render(Document &, Type<Result::MFE>);
void render(Document &, Type<Result::Pairs>);
void render(Document &, Type<Result::Sample>);
void render(Document &, Type<Result::Subopt>);

void render(Document &, Type<Future>);

}
