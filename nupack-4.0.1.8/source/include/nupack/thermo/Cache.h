#pragma once
#include "../standard/Ptr.h"
#include "../types/Complex.h"
#include "Schedule.h"

namespace nupack::thermo {

/**************************************************************************************/

struct Cache {
    virtual std::size_t shrink_to_limit(std::size_t) = 0;
    virtual std::size_t used() const = 0;
    virtual void clear() = 0;
    virtual ~Cache() {}

    static std::shared_ptr<Cache> create();
};

/**************************************************************************************/

}
