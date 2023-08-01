README file for NUPACK 4.0
Copyright (c) 2003-2022. California Institute of Technology. All Rights Reserved.

See LICENSE.txt file

See NUPACK User Guide for installation instructions:
https://docs.nupack.org

Technical support: support@nupack.org

This source folder contains two text files:
- LICENSE.txt: NUPACK Software License Agreement
- CMakeLists.txt: The main configuration script for the project (using CMake)

and several sub-directories:

- bind: Code for creating Python bindings from the NUPACK C++ library
- cmake: Configuration scripts to use with cmake
- external: External git repositories used by the C++ code
- include: Include C++ (*.h) directories and files
- package: Miscellaneous files used to build the Python package
- parameters: Free energy model parameter files
- python: Python files for the NUPACK Python module
- solver: Include and source C++ (*.h and *.cc) files for equilibrium concentration solver
- source: Source C++ (*.cc) directories and files
- test: Python test files