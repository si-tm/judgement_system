#pragma once
#include "common/Config.h"
#include "algorithms/Traits.h"
#include <nlohmann/json_fwd.hpp>

/**************************************************************************************/

namespace nupack {

using json = nlohmann::json;

/**************************************************************************************/

struct ParameterFile;

struct ParameterInfo;

template <class T>
struct ParameterSet;

template <class T>
struct ParameterArray;

template <class T=real>
struct Model;

struct System;

struct PairList;

struct Structure;

struct StateBase;

class SequenceSet;

template <class SS=SequenceSet>
struct StaticLoop;

template <class Base_=StateBase, class Loop_=StaticLoop<>>
struct StaticState;

using State = StaticState<>;

struct Local;

/**************************************************************************************/

NUPACK_NAMESPACE(kmc);
namespace kmc {
    struct JumpSequenceSet;

    template <class SS=JumpSequenceSet>
    struct JumpLoop;

    template <class Loop_=JumpLoop<>, class EM=Model<real>>
    struct JumpState;

    struct EnumeratedObserver;
    struct Stopwatch;
    struct Timer;
    struct FirstCollision;
    template <bool B=true> struct TimeIntegrator;
    template <bool B=true> struct ScaledIntegrator;
    struct HammingObserver;
    struct HittingTimeObserver;
    template <class> struct CovarianceIntegrator;
    template <class> class PairIntegrator;
    template <class> class PairProbabilityBoltzmann;
}

/**************************************************************************************/

NUPACK_NAMESPACE(design);
NUPACK_NAMESPACE(traits);
NUPACK_NAMESPACE(concentration);
NUPACK_NAMESPACE(lapack);
NUPACK_NAMESPACE(mma);
NUPACK_NAMESPACE(matlab);
NUPACK_NAMESPACE(memory);
NUPACK_NAMESPACE(simd);
NUPACK_NAMESPACE(inference);
NUPACK_NAMESPACE(schedule);
NUPACK_NAMESPACE(math);
NUPACK_NAMESPACE(tensor);

/**************************************************************************************/

NUPACK_NAMESPACE(thermo);
namespace thermo {

NUPACK_NAMESPACE(coax);

template <class T> struct PF;
template <class T> struct MFE;

template <class Rig, class Model>
class CachedModel;

struct Storage;
struct Cache;

}

/**************************************************************************************/


}
