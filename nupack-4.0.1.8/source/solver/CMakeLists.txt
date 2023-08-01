cmake_minimum_required(VERSION 3.10)
################################################################################

set(CMAKE_TOOLCHAIN_FILE ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake
  CACHE STRING "Vcpkg toolchain file")

################################################################################

enable_language(CXX)

################################################################################
project(NUPACK-solver VERSION 1.0.0.0 LANGUAGES CXX)
################################################################################

add_library(nupack-solver STATIC ${CMAKE_CURRENT_SOURCE_DIR}/source/Concentration.cc)
target_include_directories(nupack-solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(nupack-solver PUBLIC cxx_std_14)

find_package(Armadillo CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
target_link_libraries(nupack-solver PUBLIC ${ARMADILLO_LIBRARIES} fmt::fmt)

target_compile_features(nupack-solver PUBLIC cxx_std_14)

################################################################################