#pragma once
#include <string_view>
#include <magic_enum.hpp>
#include "../common/Error.h"

namespace nupack {

template <class T>
auto enum_to_string(T const &t) {return magic_enum::enum_name(t);}

template <class T>
T enum_from_string(std::string_view name) {
    auto out = magic_enum::enum_cast<T>(name);
    NUPACK_ASSERT(out, "invalid name for enum", name, magic_enum::enum_type_name<T>(), magic_enum::enum_names<T>());
    return *out;
}

}