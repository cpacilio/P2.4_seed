# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/dealii/P2.4_seed/exercises/lab03/step-3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dealii/P2.4_seed/exercises/lab03/step-3/build

# Utility rule file for strip_comments.

# Include the progress variables for this target.
include CMakeFiles/strip_comments.dir/progress.make

CMakeFiles/strip_comments:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/dealii/P2.4_seed/exercises/lab03/step-3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "strip comments"
	/usr/bin/perl -pi -e 's#^[ \t]*//.*\n##g;' step-3_v2.cc

strip_comments: CMakeFiles/strip_comments
strip_comments: CMakeFiles/strip_comments.dir/build.make

.PHONY : strip_comments

# Rule to build all files generated by this target.
CMakeFiles/strip_comments.dir/build: strip_comments

.PHONY : CMakeFiles/strip_comments.dir/build

CMakeFiles/strip_comments.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/strip_comments.dir/cmake_clean.cmake
.PHONY : CMakeFiles/strip_comments.dir/clean

CMakeFiles/strip_comments.dir/depend:
	cd /home/dealii/P2.4_seed/exercises/lab03/step-3/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dealii/P2.4_seed/exercises/lab03/step-3 /home/dealii/P2.4_seed/exercises/lab03/step-3 /home/dealii/P2.4_seed/exercises/lab03/step-3/build /home/dealii/P2.4_seed/exercises/lab03/step-3/build /home/dealii/P2.4_seed/exercises/lab03/step-3/build/CMakeFiles/strip_comments.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/strip_comments.dir/depend

