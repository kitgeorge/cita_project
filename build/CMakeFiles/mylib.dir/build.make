# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/kit/Documents/cita_project

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/kit/Documents/cita_project/build

# Include any dependencies generated for this target.
include CMakeFiles/mylib.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mylib.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mylib.dir/flags.make

CMakeFiles/mylib.dir/src/df/sample.cpp.o: CMakeFiles/mylib.dir/flags.make
CMakeFiles/mylib.dir/src/df/sample.cpp.o: ../src/df/sample.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/kit/Documents/cita_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/mylib.dir/src/df/sample.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mylib.dir/src/df/sample.cpp.o -c /home/kit/Documents/cita_project/src/df/sample.cpp

CMakeFiles/mylib.dir/src/df/sample.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mylib.dir/src/df/sample.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/kit/Documents/cita_project/src/df/sample.cpp > CMakeFiles/mylib.dir/src/df/sample.cpp.i

CMakeFiles/mylib.dir/src/df/sample.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mylib.dir/src/df/sample.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/kit/Documents/cita_project/src/df/sample.cpp -o CMakeFiles/mylib.dir/src/df/sample.cpp.s

CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.cpp.o: CMakeFiles/mylib.dir/flags.make
CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.cpp.o: ../src/multithreading/execute_in_parallel.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/kit/Documents/cita_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.cpp.o -c /home/kit/Documents/cita_project/src/multithreading/execute_in_parallel.cpp

CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/kit/Documents/cita_project/src/multithreading/execute_in_parallel.cpp > CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.cpp.i

CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/kit/Documents/cita_project/src/multithreading/execute_in_parallel.cpp -o CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.cpp.s

CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.test.cpp.o: CMakeFiles/mylib.dir/flags.make
CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.test.cpp.o: ../src/multithreading/execute_in_parallel.test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/kit/Documents/cita_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.test.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.test.cpp.o -c /home/kit/Documents/cita_project/src/multithreading/execute_in_parallel.test.cpp

CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.test.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/kit/Documents/cita_project/src/multithreading/execute_in_parallel.test.cpp > CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.test.cpp.i

CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.test.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/kit/Documents/cita_project/src/multithreading/execute_in_parallel.test.cpp -o CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.test.cpp.s

CMakeFiles/mylib.dir/src/potential/mestel.cpp.o: CMakeFiles/mylib.dir/flags.make
CMakeFiles/mylib.dir/src/potential/mestel.cpp.o: ../src/potential/mestel.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/kit/Documents/cita_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/mylib.dir/src/potential/mestel.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mylib.dir/src/potential/mestel.cpp.o -c /home/kit/Documents/cita_project/src/potential/mestel.cpp

CMakeFiles/mylib.dir/src/potential/mestel.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mylib.dir/src/potential/mestel.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/kit/Documents/cita_project/src/potential/mestel.cpp > CMakeFiles/mylib.dir/src/potential/mestel.cpp.i

CMakeFiles/mylib.dir/src/potential/mestel.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mylib.dir/src/potential/mestel.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/kit/Documents/cita_project/src/potential/mestel.cpp -o CMakeFiles/mylib.dir/src/potential/mestel.cpp.s

CMakeFiles/mylib.dir/src/potential/spiral.cpp.o: CMakeFiles/mylib.dir/flags.make
CMakeFiles/mylib.dir/src/potential/spiral.cpp.o: ../src/potential/spiral.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/kit/Documents/cita_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/mylib.dir/src/potential/spiral.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mylib.dir/src/potential/spiral.cpp.o -c /home/kit/Documents/cita_project/src/potential/spiral.cpp

CMakeFiles/mylib.dir/src/potential/spiral.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mylib.dir/src/potential/spiral.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/kit/Documents/cita_project/src/potential/spiral.cpp > CMakeFiles/mylib.dir/src/potential/spiral.cpp.i

CMakeFiles/mylib.dir/src/potential/spiral.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mylib.dir/src/potential/spiral.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/kit/Documents/cita_project/src/potential/spiral.cpp -o CMakeFiles/mylib.dir/src/potential/spiral.cpp.s

CMakeFiles/mylib.dir/src/potential/spiral.test.cpp.o: CMakeFiles/mylib.dir/flags.make
CMakeFiles/mylib.dir/src/potential/spiral.test.cpp.o: ../src/potential/spiral.test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/kit/Documents/cita_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/mylib.dir/src/potential/spiral.test.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mylib.dir/src/potential/spiral.test.cpp.o -c /home/kit/Documents/cita_project/src/potential/spiral.test.cpp

CMakeFiles/mylib.dir/src/potential/spiral.test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mylib.dir/src/potential/spiral.test.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/kit/Documents/cita_project/src/potential/spiral.test.cpp > CMakeFiles/mylib.dir/src/potential/spiral.test.cpp.i

CMakeFiles/mylib.dir/src/potential/spiral.test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mylib.dir/src/potential/spiral.test.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/kit/Documents/cita_project/src/potential/spiral.test.cpp -o CMakeFiles/mylib.dir/src/potential/spiral.test.cpp.s

CMakeFiles/mylib.dir/src/utility/add_functions.test.cpp.o: CMakeFiles/mylib.dir/flags.make
CMakeFiles/mylib.dir/src/utility/add_functions.test.cpp.o: ../src/utility/add_functions.test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/kit/Documents/cita_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/mylib.dir/src/utility/add_functions.test.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mylib.dir/src/utility/add_functions.test.cpp.o -c /home/kit/Documents/cita_project/src/utility/add_functions.test.cpp

CMakeFiles/mylib.dir/src/utility/add_functions.test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mylib.dir/src/utility/add_functions.test.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/kit/Documents/cita_project/src/utility/add_functions.test.cpp > CMakeFiles/mylib.dir/src/utility/add_functions.test.cpp.i

CMakeFiles/mylib.dir/src/utility/add_functions.test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mylib.dir/src/utility/add_functions.test.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/kit/Documents/cita_project/src/utility/add_functions.test.cpp -o CMakeFiles/mylib.dir/src/utility/add_functions.test.cpp.s

CMakeFiles/mylib.dir/src/utility/flatten.test.cpp.o: CMakeFiles/mylib.dir/flags.make
CMakeFiles/mylib.dir/src/utility/flatten.test.cpp.o: ../src/utility/flatten.test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/kit/Documents/cita_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/mylib.dir/src/utility/flatten.test.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mylib.dir/src/utility/flatten.test.cpp.o -c /home/kit/Documents/cita_project/src/utility/flatten.test.cpp

CMakeFiles/mylib.dir/src/utility/flatten.test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mylib.dir/src/utility/flatten.test.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/kit/Documents/cita_project/src/utility/flatten.test.cpp > CMakeFiles/mylib.dir/src/utility/flatten.test.cpp.i

CMakeFiles/mylib.dir/src/utility/flatten.test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mylib.dir/src/utility/flatten.test.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/kit/Documents/cita_project/src/utility/flatten.test.cpp -o CMakeFiles/mylib.dir/src/utility/flatten.test.cpp.s

CMakeFiles/mylib.dir/src/utility/shape.test.cpp.o: CMakeFiles/mylib.dir/flags.make
CMakeFiles/mylib.dir/src/utility/shape.test.cpp.o: ../src/utility/shape.test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/kit/Documents/cita_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/mylib.dir/src/utility/shape.test.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mylib.dir/src/utility/shape.test.cpp.o -c /home/kit/Documents/cita_project/src/utility/shape.test.cpp

CMakeFiles/mylib.dir/src/utility/shape.test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mylib.dir/src/utility/shape.test.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/kit/Documents/cita_project/src/utility/shape.test.cpp > CMakeFiles/mylib.dir/src/utility/shape.test.cpp.i

CMakeFiles/mylib.dir/src/utility/shape.test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mylib.dir/src/utility/shape.test.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/kit/Documents/cita_project/src/utility/shape.test.cpp -o CMakeFiles/mylib.dir/src/utility/shape.test.cpp.s

CMakeFiles/mylib.dir/src/utility/vector_io.cpp.o: CMakeFiles/mylib.dir/flags.make
CMakeFiles/mylib.dir/src/utility/vector_io.cpp.o: ../src/utility/vector_io.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/kit/Documents/cita_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/mylib.dir/src/utility/vector_io.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mylib.dir/src/utility/vector_io.cpp.o -c /home/kit/Documents/cita_project/src/utility/vector_io.cpp

CMakeFiles/mylib.dir/src/utility/vector_io.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mylib.dir/src/utility/vector_io.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/kit/Documents/cita_project/src/utility/vector_io.cpp > CMakeFiles/mylib.dir/src/utility/vector_io.cpp.i

CMakeFiles/mylib.dir/src/utility/vector_io.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mylib.dir/src/utility/vector_io.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/kit/Documents/cita_project/src/utility/vector_io.cpp -o CMakeFiles/mylib.dir/src/utility/vector_io.cpp.s

CMakeFiles/mylib.dir/src/utility/vector_io.test.cpp.o: CMakeFiles/mylib.dir/flags.make
CMakeFiles/mylib.dir/src/utility/vector_io.test.cpp.o: ../src/utility/vector_io.test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/kit/Documents/cita_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/mylib.dir/src/utility/vector_io.test.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mylib.dir/src/utility/vector_io.test.cpp.o -c /home/kit/Documents/cita_project/src/utility/vector_io.test.cpp

CMakeFiles/mylib.dir/src/utility/vector_io.test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mylib.dir/src/utility/vector_io.test.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/kit/Documents/cita_project/src/utility/vector_io.test.cpp > CMakeFiles/mylib.dir/src/utility/vector_io.test.cpp.i

CMakeFiles/mylib.dir/src/utility/vector_io.test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mylib.dir/src/utility/vector_io.test.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/kit/Documents/cita_project/src/utility/vector_io.test.cpp -o CMakeFiles/mylib.dir/src/utility/vector_io.test.cpp.s

CMakeFiles/mylib.dir/src/vectors/force.cpp.o: CMakeFiles/mylib.dir/flags.make
CMakeFiles/mylib.dir/src/vectors/force.cpp.o: ../src/vectors/force.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/kit/Documents/cita_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/mylib.dir/src/vectors/force.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mylib.dir/src/vectors/force.cpp.o -c /home/kit/Documents/cita_project/src/vectors/force.cpp

CMakeFiles/mylib.dir/src/vectors/force.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mylib.dir/src/vectors/force.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/kit/Documents/cita_project/src/vectors/force.cpp > CMakeFiles/mylib.dir/src/vectors/force.cpp.i

CMakeFiles/mylib.dir/src/vectors/force.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mylib.dir/src/vectors/force.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/kit/Documents/cita_project/src/vectors/force.cpp -o CMakeFiles/mylib.dir/src/vectors/force.cpp.s

CMakeFiles/mylib.dir/src/vectors/force.test.cpp.o: CMakeFiles/mylib.dir/flags.make
CMakeFiles/mylib.dir/src/vectors/force.test.cpp.o: ../src/vectors/force.test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/kit/Documents/cita_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/mylib.dir/src/vectors/force.test.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mylib.dir/src/vectors/force.test.cpp.o -c /home/kit/Documents/cita_project/src/vectors/force.test.cpp

CMakeFiles/mylib.dir/src/vectors/force.test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mylib.dir/src/vectors/force.test.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/kit/Documents/cita_project/src/vectors/force.test.cpp > CMakeFiles/mylib.dir/src/vectors/force.test.cpp.i

CMakeFiles/mylib.dir/src/vectors/force.test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mylib.dir/src/vectors/force.test.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/kit/Documents/cita_project/src/vectors/force.test.cpp -o CMakeFiles/mylib.dir/src/vectors/force.test.cpp.s

# Object files for target mylib
mylib_OBJECTS = \
"CMakeFiles/mylib.dir/src/df/sample.cpp.o" \
"CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.cpp.o" \
"CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.test.cpp.o" \
"CMakeFiles/mylib.dir/src/potential/mestel.cpp.o" \
"CMakeFiles/mylib.dir/src/potential/spiral.cpp.o" \
"CMakeFiles/mylib.dir/src/potential/spiral.test.cpp.o" \
"CMakeFiles/mylib.dir/src/utility/add_functions.test.cpp.o" \
"CMakeFiles/mylib.dir/src/utility/flatten.test.cpp.o" \
"CMakeFiles/mylib.dir/src/utility/shape.test.cpp.o" \
"CMakeFiles/mylib.dir/src/utility/vector_io.cpp.o" \
"CMakeFiles/mylib.dir/src/utility/vector_io.test.cpp.o" \
"CMakeFiles/mylib.dir/src/vectors/force.cpp.o" \
"CMakeFiles/mylib.dir/src/vectors/force.test.cpp.o"

# External object files for target mylib
mylib_EXTERNAL_OBJECTS =

libmylib.so: CMakeFiles/mylib.dir/src/df/sample.cpp.o
libmylib.so: CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.cpp.o
libmylib.so: CMakeFiles/mylib.dir/src/multithreading/execute_in_parallel.test.cpp.o
libmylib.so: CMakeFiles/mylib.dir/src/potential/mestel.cpp.o
libmylib.so: CMakeFiles/mylib.dir/src/potential/spiral.cpp.o
libmylib.so: CMakeFiles/mylib.dir/src/potential/spiral.test.cpp.o
libmylib.so: CMakeFiles/mylib.dir/src/utility/add_functions.test.cpp.o
libmylib.so: CMakeFiles/mylib.dir/src/utility/flatten.test.cpp.o
libmylib.so: CMakeFiles/mylib.dir/src/utility/shape.test.cpp.o
libmylib.so: CMakeFiles/mylib.dir/src/utility/vector_io.cpp.o
libmylib.so: CMakeFiles/mylib.dir/src/utility/vector_io.test.cpp.o
libmylib.so: CMakeFiles/mylib.dir/src/vectors/force.cpp.o
libmylib.so: CMakeFiles/mylib.dir/src/vectors/force.test.cpp.o
libmylib.so: CMakeFiles/mylib.dir/build.make
libmylib.so: CMakeFiles/mylib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/kit/Documents/cita_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Linking CXX shared library libmylib.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mylib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mylib.dir/build: libmylib.so

.PHONY : CMakeFiles/mylib.dir/build

CMakeFiles/mylib.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mylib.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mylib.dir/clean

CMakeFiles/mylib.dir/depend:
	cd /home/kit/Documents/cita_project/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/kit/Documents/cita_project /home/kit/Documents/cita_project /home/kit/Documents/cita_project/build /home/kit/Documents/cita_project/build /home/kit/Documents/cita_project/build/CMakeFiles/mylib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mylib.dir/depend

