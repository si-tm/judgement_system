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
include contrib/rovigatti/CMakeFiles/AOInteraction.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include contrib/rovigatti/CMakeFiles/AOInteraction.dir/compiler_depend.make

# Include the progress variables for this target.
include contrib/rovigatti/CMakeFiles/AOInteraction.dir/progress.make

# Include the compile flags for this target's objects.
include contrib/rovigatti/CMakeFiles/AOInteraction.dir/flags.make

contrib/rovigatti/CMakeFiles/AOInteraction.dir/src/Interactions/AOInteraction.cpp.o: contrib/rovigatti/CMakeFiles/AOInteraction.dir/flags.make
contrib/rovigatti/CMakeFiles/AOInteraction.dir/src/Interactions/AOInteraction.cpp.o: ../contrib/rovigatti/src/Interactions/AOInteraction.cpp
contrib/rovigatti/CMakeFiles/AOInteraction.dir/src/Interactions/AOInteraction.cpp.o: contrib/rovigatti/CMakeFiles/AOInteraction.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object contrib/rovigatti/CMakeFiles/AOInteraction.dir/src/Interactions/AOInteraction.cpp.o"
	cd "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/rovigatti" && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT contrib/rovigatti/CMakeFiles/AOInteraction.dir/src/Interactions/AOInteraction.cpp.o -MF CMakeFiles/AOInteraction.dir/src/Interactions/AOInteraction.cpp.o.d -o CMakeFiles/AOInteraction.dir/src/Interactions/AOInteraction.cpp.o -c "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/contrib/rovigatti/src/Interactions/AOInteraction.cpp"

contrib/rovigatti/CMakeFiles/AOInteraction.dir/src/Interactions/AOInteraction.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/AOInteraction.dir/src/Interactions/AOInteraction.cpp.i"
	cd "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/rovigatti" && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/contrib/rovigatti/src/Interactions/AOInteraction.cpp" > CMakeFiles/AOInteraction.dir/src/Interactions/AOInteraction.cpp.i

contrib/rovigatti/CMakeFiles/AOInteraction.dir/src/Interactions/AOInteraction.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/AOInteraction.dir/src/Interactions/AOInteraction.cpp.s"
	cd "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/rovigatti" && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/contrib/rovigatti/src/Interactions/AOInteraction.cpp" -o CMakeFiles/AOInteraction.dir/src/Interactions/AOInteraction.cpp.s

# Object files for target AOInteraction
AOInteraction_OBJECTS = \
"CMakeFiles/AOInteraction.dir/src/Interactions/AOInteraction.cpp.o"

# External object files for target AOInteraction
AOInteraction_EXTERNAL_OBJECTS =

../contrib/rovigatti/AOInteraction.so: contrib/rovigatti/CMakeFiles/AOInteraction.dir/src/Interactions/AOInteraction.cpp.o
../contrib/rovigatti/AOInteraction.so: contrib/rovigatti/CMakeFiles/AOInteraction.dir/build.make
../contrib/rovigatti/AOInteraction.so: contrib/rovigatti/CMakeFiles/AOInteraction.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared module ../../../contrib/rovigatti/AOInteraction.so"
	cd "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/rovigatti" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/AOInteraction.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
contrib/rovigatti/CMakeFiles/AOInteraction.dir/build: ../contrib/rovigatti/AOInteraction.so
.PHONY : contrib/rovigatti/CMakeFiles/AOInteraction.dir/build

contrib/rovigatti/CMakeFiles/AOInteraction.dir/clean:
	cd "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/rovigatti" && $(CMAKE_COMMAND) -P CMakeFiles/AOInteraction.dir/cmake_clean.cmake
.PHONY : contrib/rovigatti/CMakeFiles/AOInteraction.dir/clean

contrib/rovigatti/CMakeFiles/AOInteraction.dir/depend:
	cd "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script" "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/contrib/rovigatti" "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build" "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/rovigatti" "/Users/maya/OneDrive - お茶の水女子大学/lab/judgement_system/script/build/contrib/rovigatti/CMakeFiles/AOInteraction.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : contrib/rovigatti/CMakeFiles/AOInteraction.dir/depend

