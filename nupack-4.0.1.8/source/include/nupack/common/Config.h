/**
 * @brief Defines type aliases and a few global constants
 *
 * @file Config.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once

#include <cstddef>
#include <cstring>
#include <string>
#include <string_view>
#include <iomanip>
#include <complex>
#include <limits>
#include <vector>

#include <boost/container/small_vector.hpp>

/******************************************************************************************/

/// These overloads currently interfere with some printing functionality
#define BOOST_NO_IOSTREAM

/******************************************************************************************/

/// File name macro which excludes directory path
#define NUPACK_FILE (std::strrchr(__FILE__, '/') ? std::strrchr(__FILE__, '/') + 1 : \
    std::strrchr(__FILE__, '\\') ? std::strrchr(__FILE__, '\\') + 1 : __FILE__)

/******************************************************************************************/

namespace nupack {

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    static constexpr bool is_windows = true;
#else
    static constexpr bool is_windows = false;
#endif

/******************************************************************************************/

/// Some simple default types for use across the whole project
using real         = double;

using real32       = float;
static_assert(sizeof(real32) == 4);

using real64       = double;
static_assert(sizeof(real64) == 8);

using string       = std::string;
using string_view  = std::string_view;

using namespace std::string_literals;

using uint         = std::uint32_t;
using complex_real = std::complex<real>;

using iseq         = std::uint_fast32_t;
static_assert(std::numeric_limits<iseq>::max() >= 1e6, "iseq constrains maximum sequence length");

using usize        = std::size_t;
static_assert(sizeof(usize) >= sizeof(std::size_t), "usize should be at least size_t size");

/******************************************************************************************/

template <class T, int N=-1, class Alloc=std::allocator<T>>
using small_vec = boost::container::small_vector<T, N == -1 ? 16/sizeof(T) : N, Alloc>;

// Vector class for type T, suggested length N, allocator Alloc
template <class T, int N=-1, class Alloc=std::allocator<T>>
using vec = std::conditional_t<(N == -1 && sizeof(T) <= 4) || (N * sizeof(T) <= 40),
    small_vec<T, N, Alloc>, std::vector<T, Alloc>>;

template <class T>
using default_vec = vec<T>;

/******************************************************************************************/

/// Some string values that will be set by CMake
extern string const GitRevision, GitBranch, Version, BuildType;
extern string DefaultParametersPath, DefaultDataPath, MatlabCommand, MathematicaCommand;

/// Number of logical CPU cores - you can change this if desired
extern unsigned int TotalCPU;

/// Total RAM in bytes - you can change this if desired
extern std::size_t TotalRAM;

/// Print backtraces in exceptions
extern bool DebugInfo;

/******************************************************************************************/

/// constexpr Debug, should be optimized out by compiler when put in if()
#if NUPACK_DEBUG
#   if NUPACK_DEBUG == 1
        static constexpr bool const Release = false;
        static constexpr bool const Debug = false;
        static constexpr bool const DebugBounds = false;
#   endif
#   if NUPACK_DEBUG == 2
        static constexpr bool const Release = false;
        static constexpr bool const Debug = true;
        static constexpr bool const DebugBounds = true;
#   endif
#else
    static constexpr bool const Release = true;
    static constexpr bool const Debug = false;
    static constexpr bool const DebugBounds = false;
#endif

/******************************************************************************************/

static constexpr std::size_t CharCapacity = std::size_t(std::numeric_limits<unsigned char>::max()) + 1;
static_assert(CharCapacity == 256);

}
