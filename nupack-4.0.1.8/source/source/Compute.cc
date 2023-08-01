#include <nupack/thermo/Strategy.h>
#include <nupack/thermo/Big.h>

namespace nupack::thermo {

ComputerList DefaultFactory<Model<>>::create(Model<> m) {
    ComputerList v;
    fork(m.ensemble_type(), [&](auto d) {
        v.emplace_back(make_computer<DynamicStrategy<Model<real64>, Model<real32>, decltype(d), 3>>(m));
        v.emplace_back(make_computer<MFEStrategy<Model<real32>, decltype(d), 3>>(m));
    });
    return v;
};

}
