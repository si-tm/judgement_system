# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.21.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.21.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build"

# Include any dependencies generated for this target.
include contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/compiler_depend.make

# Include the progress variables for this target.
include contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/progress.make

# Include the compile flags for this target's objects.
include contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/flags.make

contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/src/Observables/TSPAnalysis.cpp.o: contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/flags.make
contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/src/Observables/TSPAnalysis.cpp.o: ../contrib/rovigatti/src/Observables/TSPAnalysis.cpp
contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/src/Observables/TSPAnalysis.cpp.o: contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/src/Observables/TSPAnalysis.cpp.o"
	cd "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/rovigatti" && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/src/Observables/TSPAnalysis.cpp.o -MF CMakeFiles/TSPAnalysis.dir/src/Observables/TSPAnalysis.cpp.o.d -o CMakeFiles/TSPAnalysis.dir/src/Observables/TSPAnalysis.cpp.o -c "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/contrib/rovigatti/src/Observables/TSPAnalysis.cpp"

contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/src/Observables/TSPAnalysis.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TSPAnalysis.dir/src/Observables/TSPAnalysis.cpp.i"
	cd "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/rovigatti" && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/contrib/rovigatti/src/Observables/TSPAnalysis.cpp" > CMakeFiles/TSPAnalysis.dir/src/Observables/TSPAnalysis.cpp.i

contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/src/Observables/TSPAnalysis.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TSPAnalysis.dir/src/Observables/TSPAnalysis.cpp.s"
	cd "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/rovigatti" && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/contrib/rovigatti/src/Observables/TSPAnalysis.cpp" -o CMakeFiles/TSPAnalysis.dir/src/Observables/TSPAnalysis.cpp.s

# Object files for target TSPAnalysis
TSPAnalysis_OBJECTS = \
"CMakeFiles/TSPAnalysis.dir/src/Observables/TSPAnalysis.cpp.o"

# External object files for target TSPAnalysis
TSPAnalysis_EXTERNAL_OBJECTS =

../contrib/rovigatti/TSPAnalysis.so: contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/src/Observables/TSPAnalysis.cpp.o
../contrib/rovigatti/TSPAnalysis.so: contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/build.make
../contrib/rovigatti/TSPAnalysis.so: contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared module ../../../contrib/rovigatti/TSPAnalysis.so"
	cd "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/rovigatti" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TSPAnalysis.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/build: ../contrib/rovigatti/TSPAnalysis.so
.PHONY : contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/build

contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/clean:
	cd "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/rovigatti" && $(CMAKE_COMMAND) -P CMakeFiles/TSPAnalysis.dir/cmake_clean.cmake
.PHONY : contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/clean

contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/depend:
	cd "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script" "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/contrib/rovigatti" "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build" "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/rovigatti" "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : contrib/rovigatti/CMakeFiles/TSPAnalysis.dir/depend

