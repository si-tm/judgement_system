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
include contrib/romano/CMakeFiles/NematicS.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include contrib/romano/CMakeFiles/NematicS.dir/compiler_depend.make

# Include the progress variables for this target.
include contrib/romano/CMakeFiles/NematicS.dir/progress.make

# Include the compile flags for this target's objects.
include contrib/romano/CMakeFiles/NematicS.dir/flags.make

contrib/romano/CMakeFiles/NematicS.dir/src/Observables/NematicS.cpp.o: contrib/romano/CMakeFiles/NematicS.dir/flags.make
contrib/romano/CMakeFiles/NematicS.dir/src/Observables/NematicS.cpp.o: ../contrib/romano/src/Observables/NematicS.cpp
contrib/romano/CMakeFiles/NematicS.dir/src/Observables/NematicS.cpp.o: contrib/romano/CMakeFiles/NematicS.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object contrib/romano/CMakeFiles/NematicS.dir/src/Observables/NematicS.cpp.o"
	cd "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/romano" && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT contrib/romano/CMakeFiles/NematicS.dir/src/Observables/NematicS.cpp.o -MF CMakeFiles/NematicS.dir/src/Observables/NematicS.cpp.o.d -o CMakeFiles/NematicS.dir/src/Observables/NematicS.cpp.o -c "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/contrib/romano/src/Observables/NematicS.cpp"

contrib/romano/CMakeFiles/NematicS.dir/src/Observables/NematicS.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NematicS.dir/src/Observables/NematicS.cpp.i"
	cd "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/romano" && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/contrib/romano/src/Observables/NematicS.cpp" > CMakeFiles/NematicS.dir/src/Observables/NematicS.cpp.i

contrib/romano/CMakeFiles/NematicS.dir/src/Observables/NematicS.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NematicS.dir/src/Observables/NematicS.cpp.s"
	cd "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/romano" && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/contrib/romano/src/Observables/NematicS.cpp" -o CMakeFiles/NematicS.dir/src/Observables/NematicS.cpp.s

# Object files for target NematicS
NematicS_OBJECTS = \
"CMakeFiles/NematicS.dir/src/Observables/NematicS.cpp.o"

# External object files for target NematicS
NematicS_EXTERNAL_OBJECTS =

../contrib/romano/NematicS.dylib: contrib/romano/CMakeFiles/NematicS.dir/src/Observables/NematicS.cpp.o
../contrib/romano/NematicS.dylib: contrib/romano/CMakeFiles/NematicS.dir/build.make
../contrib/romano/NematicS.dylib: contrib/romano/CMakeFiles/NematicS.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library ../../../contrib/romano/NematicS.dylib"
	cd "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/romano" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/NematicS.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
contrib/romano/CMakeFiles/NematicS.dir/build: ../contrib/romano/NematicS.dylib
.PHONY : contrib/romano/CMakeFiles/NematicS.dir/build

contrib/romano/CMakeFiles/NematicS.dir/clean:
	cd "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/romano" && $(CMAKE_COMMAND) -P CMakeFiles/NematicS.dir/cmake_clean.cmake
.PHONY : contrib/romano/CMakeFiles/NematicS.dir/clean

contrib/romano/CMakeFiles/NematicS.dir/depend:
	cd "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script" "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/contrib/romano" "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build" "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/romano" "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/romano/CMakeFiles/NematicS.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : contrib/romano/CMakeFiles/NematicS.dir/depend

