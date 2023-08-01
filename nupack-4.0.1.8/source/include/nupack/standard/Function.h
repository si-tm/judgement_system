#pragma once
//#include "../reflect/Memory.h"
//#include "../reflect/Print.h"
//#include "../reflect/Hash.h"

#include <functional>

namespace nupack {

/******************************************************************************************/

template <class R, class ...Ts>
using Function = std::function<R(Ts...)>;

// NUPACK_DEFINE_VARIADIC(is_function, Function, class);

template <class ...Ts>
struct Observer {
    Function<void, Ts...> function;
    void operator()(Ts ...ts) const {if (function) function(static_cast<Ts>(ts)...);}
};

/******************************************************************************************/

template <class R, class ...Ts>
struct memory::impl<std::function<R(Ts...)>> {
    constexpr std::size_t operator()(std::function<R(Ts...)> const &t) const {
        return sizeof(t);
    }
    void erase(std::function<R(Ts...)> &t) const {t = {};}
};

}

