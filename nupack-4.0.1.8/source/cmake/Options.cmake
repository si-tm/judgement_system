option(NUPACK_BUILD_SUPER        "Whether or not a superbuild should be invoked" OFF)
option(NUPACK_BUILD_CXX          "Whether or not to build C++ libraries"         ON)
option(NUPACK_BUILD_TESTS        "Whether or not to build C++ unit tests"        ON)
option(NUPACK_BUILD_PYTHON       "Whether or not to build Python bindings"       ON)
option(NUPACK_BUILD_DOCS         "Whether or not to build website documentation" OFF)
option(NUPACK_CCACHE             "Enable CCache"                                 ON)
option(NUPACK_SERIALIZE          "Enable serialization via JSON and msgpack"     ON)
option(NUPACK_PIC                "Use position independent code"                 ON)
option(NUPACK_DETERMINISTIC      "Use deterministic random number generator"     OFF)
option(NUPACK_FORTRAN            "Compile Fortran sandbox"                       OFF)
option(NUPACK_IWYU               "Enable include-what-you-use"                   OFF)
option(NUPACK_HDF5               "Enable HDF5"                                   OFF)
option(NUPACK_MATHEMATICA        "Enable Mathematica"                            OFF)
option(NUPACK_MATLAB             "Enable MATLAB"                                 OFF)
option(NUPACK_MLPACK             "Enable MLPACK"                                 OFF)
option(NUPACK_MPI                "Enable MPI compilation"                        OFF)
option(NUPACK_OMP                "Enable OpenMP"                                 OFF)
option(NUPACK_ONLY               "Compile only a given file"                     OFF)
option(NUPACK_PGO                "Enable PGO, can be OFF, READ, or WRITE"        OFF)
option(NUPACK_EXTERNAL_ARMADILLO "Use external version of armadillo"             OFF)
option(NUPACK_SIMD               "Use SIMD instructions for dynamic programs"    ON)

################################################################################

message(STATUS "--------------------------------------------------------------------------------")

################################################################################

cmake_host_system_information(RESULT ram_in_mb QUERY TOTAL_PHYSICAL_MEMORY)
set(NUPACK_RAM_IN_MB  ${ram_in_mb} CACHE STRING "Maximum amount of RAM")
message(STATUS "-- Assuming maximum RAM to use is ${NUPACK_RAM_IN_MB} MB")

################################################################################

# for each path in a list, add a directory prefix to it
function(prefix_transform var prefix)
    set(listVar "")
    foreach(f ${ARGN})
        list(APPEND listVar "${prefix}/${f}")
    endforeach(f)
    set(${var} "${listVar}" PARENT_SCOPE)
endfunction(prefix_transform)

################################################################################

if(NUPACK_PIC)
    message("-- Compiling with position-independent code linkage")
else()
    message("-- Not compiling with position-independent code linkage")
endif()

function(nupack_cxx_target var)
    set_target_properties(${var} PROPERTIES POSITION_INDEPENDENT_CODE ${NUPACK_PIC})
endfunction()

################################################################################