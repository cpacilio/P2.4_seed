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
CMAKE_SOURCE_DIR = /home/dealii/P2.4_seed/exercises/lab02/step-2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dealii/P2.4_seed/exercises/lab02/step-2/build

# Include any dependencies generated for this target.
include CMakeFiles/step-2.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/step-2.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/step-2.dir/flags.make

CMakeFiles/step-2.dir/step-2.cc.o: CMakeFiles/step-2.dir/flags.make
CMakeFiles/step-2.dir/step-2.cc.o: ../step-2.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dealii/P2.4_seed/exercises/lab02/step-2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/step-2.dir/step-2.cc.o"
	/usr/bin/mpicxx   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/step-2.dir/step-2.cc.o -c /home/dealii/P2.4_seed/exercises/lab02/step-2/step-2.cc

CMakeFiles/step-2.dir/step-2.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/step-2.dir/step-2.cc.i"
	/usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dealii/P2.4_seed/exercises/lab02/step-2/step-2.cc > CMakeFiles/step-2.dir/step-2.cc.i

CMakeFiles/step-2.dir/step-2.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/step-2.dir/step-2.cc.s"
	/usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dealii/P2.4_seed/exercises/lab02/step-2/step-2.cc -o CMakeFiles/step-2.dir/step-2.cc.s

CMakeFiles/step-2.dir/step-2.cc.o.requires:

.PHONY : CMakeFiles/step-2.dir/step-2.cc.o.requires

CMakeFiles/step-2.dir/step-2.cc.o.provides: CMakeFiles/step-2.dir/step-2.cc.o.requires
	$(MAKE) -f CMakeFiles/step-2.dir/build.make CMakeFiles/step-2.dir/step-2.cc.o.provides.build
.PHONY : CMakeFiles/step-2.dir/step-2.cc.o.provides

CMakeFiles/step-2.dir/step-2.cc.o.provides.build: CMakeFiles/step-2.dir/step-2.cc.o


# Object files for target step-2
step__2_OBJECTS = \
"CMakeFiles/step-2.dir/step-2.cc.o"

# External object files for target step-2
step__2_EXTERNAL_OBJECTS =

step-2: CMakeFiles/step-2.dir/step-2.cc.o
step-2: CMakeFiles/step-2.dir/build.make
step-2: /home/dealii/dealii-v9.0.0/lib/libdeal_II.g.so.9.0.0
step-2: /home/dealii/libs/p4est-2.0/DEBUG/lib/libp4est.so
step-2: /home/dealii/libs/p4est-2.0/DEBUG/lib/libsc.so
step-2: /usr/lib/x86_64-linux-gnu/libz.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libmuelu-adapters.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libmuelu-interface.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libmuelu.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libteko.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libstratimikos.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libstratimikosbelos.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libstratimikosaztecoo.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libstratimikosamesos.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libstratimikosml.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libstratimikosifpack.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libifpack2-adapters.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libifpack2.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libanasazitpetra.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libModeLaplace.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libanasaziepetra.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libanasazi.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libamesos2.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libbelostpetra.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libbelosepetra.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libbelos.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libml.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libifpack.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libzoltan2.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libpamgen_extras.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libpamgen.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libamesos.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libgaleri-xpetra.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libgaleri-epetra.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libaztecoo.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libisorropia.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libxpetra-sup.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libxpetra.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libthyratpetra.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libthyraepetraext.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libthyraepetra.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libthyracore.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libepetraext.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libtpetraext.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libtpetrainout.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libtpetra.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libkokkostsqr.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libtpetrakernels.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libtpetraclassiclinalg.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libtpetraclassicnodeapi.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libtpetraclassic.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libtriutils.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libzoltan.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libepetra.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libsacado.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/librtop.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libteuchoskokkoscomm.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libteuchoskokkoscompat.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libteuchosremainder.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libteuchosnumerics.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libteuchoscomm.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libteuchosparameterlist.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libteuchoscore.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libkokkosalgorithms.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libkokkoscontainers.so
step-2: /home/dealii/libs/trilinos-release-12-10-1/lib/libkokkoscore.so
step-2: /home/dealii/libs/superlu_dist_5.1.2/lib/libsuperlu_dist.so
step-2: /home/dealii/libs/adolc-2.6.4-rc1/lib64/libadolc.so
step-2: /home/dealii/libs/arpack-ng-3.6.2/lib/libarpack.so
step-2: /home/dealii/libs/assimp-3.3.1/lib/libassimp.so
step-2: /usr/lib/x86_64-linux-gnu/libgsl.so
step-2: /usr/lib/x86_64-linux-gnu/libgslcblas.so
step-2: /home/dealii/libs/hdf5-1.10.1/lib/libhdf5_hl.so
step-2: /home/dealii/libs/hdf5-1.10.1/lib/libhdf5.so
step-2: /usr/lib/x86_64-linux-gnu/libnetcdf_c++.so
step-2: /usr/lib/x86_64-linux-gnu/libnetcdf.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKBO.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKBool.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKBRep.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKernel.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKFeat.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKFillet.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKG2d.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKG3d.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKGeomAlgo.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKGeomBase.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKHLR.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKIGES.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKMath.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKMesh.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKOffset.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKPrim.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKShHealing.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKSTEP.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKSTEPAttr.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKSTEPBase.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKSTEP209.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKSTL.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKTopAlgo.so
step-2: /home/dealii/libs/oce-OCE-0.18.2/lib/libTKXSBase.so
step-2: /home/dealii/libs/slepc-3.7.3/lib/libslepc.so
step-2: /home/dealii/libs/petsc-3.7.6/lib/libpetsc.so
step-2: /home/dealii/libs/petsc-3.7.6/lib/libcmumps.a
step-2: /home/dealii/libs/petsc-3.7.6/lib/libdmumps.a
step-2: /home/dealii/libs/petsc-3.7.6/lib/libsmumps.a
step-2: /home/dealii/libs/petsc-3.7.6/lib/libzmumps.a
step-2: /home/dealii/libs/petsc-3.7.6/lib/libmumps_common.a
step-2: /home/dealii/libs/petsc-3.7.6/lib/libpord.a
step-2: /home/dealii/libs/parmetis-4.0.3/lib/libparmetis.so
step-2: /home/dealii/libs/parmetis-4.0.3/lib/libmetis.so
step-2: /home/dealii/libs/petsc-3.7.6/lib/libHYPRE.a
step-2: /home/dealii/libs/petsc-3.7.6/lib/libscalapack.a
step-2: /usr/lib/liblapack.so
step-2: /usr/lib/libblas.so
step-2: /usr/lib/x86_64-linux-gnu/libhwloc.so
step-2: /usr/lib/openmpi/lib/libmpi_usempif08.so
step-2: /usr/lib/openmpi/lib/libmpi_usempi_ignore_tkr.so
step-2: /usr/lib/openmpi/lib/libmpi_mpifh.so
step-2: /usr/lib/openmpi/lib/libmpi_cxx.so
step-2: /usr/lib/openmpi/lib/libmpi.so
step-2: /home/dealii/libs/sundials-3.1.0/lib/libsundials_idas.so
step-2: /home/dealii/libs/sundials-3.1.0/lib/libsundials_arkode.so
step-2: /home/dealii/libs/sundials-3.1.0/lib/libsundials_kinsol.so
step-2: /home/dealii/libs/sundials-3.1.0/lib/libsundials_nvecserial.so
step-2: /home/dealii/libs/sundials-3.1.0/lib/libsundials_nvecparallel.so
step-2: CMakeFiles/step-2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/dealii/P2.4_seed/exercises/lab02/step-2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable step-2"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/step-2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/step-2.dir/build: step-2

.PHONY : CMakeFiles/step-2.dir/build

CMakeFiles/step-2.dir/requires: CMakeFiles/step-2.dir/step-2.cc.o.requires

.PHONY : CMakeFiles/step-2.dir/requires

CMakeFiles/step-2.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/step-2.dir/cmake_clean.cmake
.PHONY : CMakeFiles/step-2.dir/clean

CMakeFiles/step-2.dir/depend:
	cd /home/dealii/P2.4_seed/exercises/lab02/step-2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dealii/P2.4_seed/exercises/lab02/step-2 /home/dealii/P2.4_seed/exercises/lab02/step-2 /home/dealii/P2.4_seed/exercises/lab02/step-2/build /home/dealii/P2.4_seed/exercises/lab02/step-2/build /home/dealii/P2.4_seed/exercises/lab02/step-2/build/CMakeFiles/step-2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/step-2.dir/depend

