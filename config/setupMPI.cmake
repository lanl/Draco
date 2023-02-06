# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   config/setupMPI.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2016 Sep 22
# brief  Setup MPI Vendors
# note   Copyright (C) 2014-2023 Triad National Security, LLC., All rights reserved.
#
# Try to find MPI in the default locations (look for mpic++ in PATH)
#
# See cmake --help-module FindMPI for details on variables set and published targets. Additionally,
# this module will set the following variables:
#
# * DRACO_C4   MPI|SCALAR
# * C4_SCALAR  BOOL
# * C4_MPI     BOOL
# * MPI_FLAVOR openmpi|mpih|intel|mvapich2|spectrum|msmpi
# * MPI_VERSION NN.NN.NN
#
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)
include(FeatureSummary)

# cmake-lint: disable=R0912,R0915,W0106
#
# * too many branches
# * function too long

# ------------------------------------------------------------------------------------------------ #
# Set MPI flavor and vendor version
#
# Returns (as cache variables)
#
# * MPI_VERSION
# * MPI_FLAVOR = {openmpi, mpich, cray, spectrum, mvapich2, intel}
#
# ------------------------------------------------------------------------------------------------ #
function(setMPIflavorVer)

  # First attempt to determine MPI flavor -- scape flavor from full path (this ususally works for
  # HPC or systems with modules)
  if(CMAKE_CXX_COMPILER_WRAPPER STREQUAL CrayPrgEnv)
    set(MPI_FLAVOR "cray")
  elseif(
    "${MPIEXEC_EXECUTABLE}" MATCHES "openmpi"
    OR "${MPIEXEC_EXECUTABLE}" MATCHES "smpi"
    OR ("${MPIEXEC_EXECUTABLE}" MATCHES "srun" AND "${MPI_C_COMPILER}" MATCHES "openmpi"))
    set(MPI_FLAVOR "openmpi")
  elseif("${MPIEXEC_EXECUTABLE}" MATCHES "mpich" OR "${MPI_C_HEADER_DIR}" MATCHES "mpich")
    set(MPI_FLAVOR "mpich")
  elseif("${MPIEXEC_EXECUTABLE}" MATCHES "impi" OR "${MPIEXEC_EXECUTABLE}" MATCHES "clusterstudio")
    set(MPI_FLAVOR "intel")
  elseif("${MPIEXEC_EXECUTABLE}" MATCHES "mvapich2")
    set(MPI_FLAVOR "mvapich2")
  elseif(
    "${MPIEXEC_EXECUTABLE}" MATCHES "spectrum-mpi"
    OR "${MPIEXEC_EXECUTABLE}" MATCHES "lrun"
    OR "${MPIEXEC_EXECUTABLE}" MATCHES "jsrun")
    set(MPI_FLAVOR "spectrum")
  endif()

  if(CMAKE_CXX_COMPILER_WRAPPER STREQUAL CrayPrgEnv)
    if(DEFINED ENV{CRAY_MPICH2_VER})
      set(MPI_VERSION $ENV{CRAY_MPICH2_VER})
    endif()
  elseif(DEFINED ENV{LMOD_MPI_VERSION})
    # Toss3 with srun
    string(REGEX REPLACE "-[a-z0-9]+" "" MPI_VERSION "$ENV{LMOD_MPI_VERSION}")
  elseif(DEFINED ENV{LMOD_FAMILY_MPI_VERSION})
    # ATS-2 with lrun
    string(REGEX REPLACE "-[a-z0-9]+" "" MPI_VERSION "$ENV{LMOD_FAMILY_MPI_VERSION}")
  else()
    execute_process(
      COMMAND ${MPIEXEC_EXECUTABLE} --version
      OUTPUT_VARIABLE DBS_MPI_VER_OUT
      ERROR_VARIABLE DBS_MPI_VER_ERR)
    set(DBS_MPI_VER "${DBS_MPI_VER_OUT}${DBS_MPI_VER_ERR}")
    string(REPLACE "\n" ";" TEMP ${DBS_MPI_VER})
    foreach(line ${TEMP})
      # extract the version...
      if(${line} MATCHES "Version"
         OR ${line} MATCHES "OpenRTE"
         OR ${line} MATCHES "Open MPI"
         OR ${line} MATCHES "Spectrum MPI")
        set(DBS_MPI_VER "${line}")
        if("${DBS_MPI_VER}" MATCHES "[0-9]+[.][0-9]+[.][0-9]+")
          string(REGEX REPLACE ".*([0-9]+)[.]([0-9]+)[.]([0-9]+).*" "\\1" DBS_MPI_VER_MAJOR
                               ${DBS_MPI_VER})
          string(REGEX REPLACE ".*([0-9]+)[.]([0-9]+)[.]([0-9]+).*" "\\2" DBS_MPI_VER_MINOR
                               ${DBS_MPI_VER})
          string(REGEX REPLACE ".*([0-9]+)[.]([0-9]+)[.]([0-9]+).*" "\\3" DBS_MPI_VER_PATCH
                               ${DBS_MPI_VER})
          set(MPI_VERSION "${DBS_MPI_VER_MAJOR}.${DBS_MPI_VER_MINOR}.${DBS_MPI_VER_PATCH}")
        elseif("${DBS_MPI_VER}" MATCHES "[0-9]+[.][0-9]+")
          string(REGEX REPLACE ".*([0-9]+)[.]([0-9]+).*" "\\1" DBS_MPI_VER_MAJOR ${DBS_MPI_VER})
          string(REGEX REPLACE ".*([0-9]+)[.]([0-9]+).*" "\\2" DBS_MPI_VER_MINOR ${DBS_MPI_VER})
          set(MPI_VERSION "${DBS_MPI_VER_MAJOR}.${DBS_MPI_VER_MINOR}")
        endif()
      endif()

      # if needed, make a 2nd pass at identifying the MPI flavor
      if(NOT DEFINED MPI_FLAVOR)
        if("${line}" MATCHES "HYDRA")
          set(MPI_FLAVOR "mpich")
        elseif("${line}" MATCHES "OpenRTE")
          set(MPI_FLAVOR "openmpi")
        elseif("${line}" MATCHES "intel-mpi" OR "${line}" MATCHES "Intel[(]R[)] MPI Library")
          set(MPI_FLAVOR "intel")
        endif()
      endif()

      # Once we have the needed information stop parsing...
      if(DEFINED MPI_FLAVOR AND DEFINED MPI_VERSION)
        break()
      endif()
    endforeach()

  endif()

  # if the FindMPI.cmake module didn't set the version, then try to do so here.
  if(NOT DEFINED MPI_VERSION AND DEFINED MPI_C_VERSION)
    set(MPI_VERSION ${MPI_C_VERSION})
  endif()

  # Return the discovered values to the calling scope
  if(DEFINED MPI_FLAVOR)
    set(MPI_FLAVOR
        "${MPI_FLAVOR}"
        CACHE STRING "Vendor brand of MPI" FORCE)
  endif()
  if(DEFINED MPI_VERSION)
    set(MPI_VERSION
        "${MPI_VERSION}"
        CACHE STRING "Vendor version of MPI" FORCE)
  endif()

endfunction()

# ------------------------------------------------------------------------------------------------ #
# Set Draco specific MPI variables
# ------------------------------------------------------------------------------------------------ #
macro(setupDracoMPIVars)

  # Set Draco build system variables based on what we know about MPI.
  if(MPI_CXX_FOUND OR MPI_C_FOUND)
    set(DRACO_C4 "MPI")
  else()
    set(DRACO_C4 "SCALAR")
  endif()

  # Save the result in the cache file.
  set(DRACO_C4
      "${DRACO_C4}"
      CACHE STRING "C4 communication mode (SCALAR or MPI)")
  # Provide a constrained pull down list in cmake-gui
  set_property(CACHE DRACO_C4 PROPERTY STRINGS SCALAR MPI)

  if("${DRACO_C4}" MATCHES "MPI")
    set(C4_SCALAR 0)
    set(C4_MPI 1)
  else()
    set(C4_SCALAR 1)
    set(C4_MPI 0)
  endif()
  set(C4_SCALAR
      ${C4_SCALAR}
      CACHE STRING "Are we building a scalar-only version (no mpi in c4/config.h)?" FORCE)
  set(C4_MPI
      ${C4_MPI}
      CACHE STRING "Are we building an MPI aware version? (c4/config.h)" FORCE)
  mark_as_advanced(C4_MPI C4_SCALAR)

endmacro()

# ------------------------------------------------------------------------------------------------ #
# Query CPU topology
#
# Returns:
#
# * MPI_CORES_PER_CPU
# * MPI_CPUS_PER_NODE
# * MPI_PHYSICAL_CORES
# * MPI_MAX_NUMPROCS_PHYSICAL
# * MPI_HYPERTHREADING
#
# See also:
#
# * Try running 'lstopo' for a graphical view of the local topology or 'lscpu' for a text version.
# * EAP's flags can be found in Test.rh/General/run_job.pl (look for $other_args).  In particular,
#   it may be useful to examine EAP's options for srun or aprun.
# ------------------------------------------------------------------------------------------------
# #
macro(query_topology)

  # These cmake commands, while useful, don't provide the topology detail that we are interested in
  # (i.e. number of sockets per node). We could use the results of these queries to know if
  # hyper-threading is enabled (if logical != physical cores)
  #
  # * cmake_host_system_information(RESULT MPI_PHYSICAL_CORES QUERY NUMBER_OF_PHYSICAL_CORES)
  # * cmake_host_system_information(RESULT MPI_LOGICAL_CORES  QUERY NUMBER_OF_LOGICAL_CORES)

  # start with default values
  set(MPI_CORES_PER_CPU 4)
  set(MPI_PHYSICAL_CORES 1)

  if("${SITENAME}" STREQUAL "Trinitite" OR "${SITENAME}" STREQUAL "Trinity")
    # Backend is different than build-node
    if($ENV{CRAY_CPU_TARGET} MATCHES "mic-knl")
      set(MPI_CORES_PER_CPU 17)
      set(MPI_PHYSICAL_CORES 4)
      set(MPIEXEC_MAX_NUMPROCS
          68
          CACHE STRING "Max procs on node." FORCE)
    else()
      set(MPI_CORES_PER_CPU 16)
      set(MPI_PHYSICAL_CORES 2)
      set(MPIEXEC_MAX_NUMPROCS
          32
          CACHE STRING "Max procs on node." FORCE)
    endif()
  elseif(SITENAME MATCHES "rznevada" OR SITENAME MATCHES "rzvernal")
    set(MPI_CORES_PER_CPU 2)
    set(MPI_PHYSICAL_CORES 64)
    set(MPIEXEC_MAX_NUMPROCS
        64
        CACHE STRING "Max procs on node." FORCE)
  elseif(EXISTS "/proc/cpuinfo")
    # read the system's cpuinfo...
    file(READ "/proc/cpuinfo" cpuinfo_data)
    string(REGEX REPLACE "\n" ";" cpuinfo_data "${cpuinfo_data}")
    foreach(line ${cpuinfo_data})
      if("${line}" MATCHES "cpu cores")
        string(REGEX REPLACE ".* ([0-9]+).*" "\\1" MPI_CORES_PER_CPU "${line}")
      elseif("${line}" MATCHES "physical id")
        string(REGEX REPLACE ".* ([0-9]+).*" "\\1" tmp "${line}")
        if(${tmp} GREATER ${MPI_PHYSICAL_CORES})
          set(MPI_PHYSICAL_CORES ${tmp})
        endif()
      endif()
    endforeach()
    # correct 0-based indexing
    math(EXPR MPI_PHYSICAL_CORES "${MPI_PHYSICAL_CORES} + 1")
  endif()

  math(EXPR MPI_CPUS_PER_NODE "${MPIEXEC_MAX_NUMPROCS} / ${MPI_CORES_PER_CPU}")
  set(MPI_CPUS_PER_NODE
      ${MPI_CPUS_PER_NODE}
      CACHE STRING "Number of multi-core CPUs per node" FORCE)
  set(MPI_CORES_PER_CPU
      ${MPI_CORES_PER_CPU}
      CACHE STRING "Number of cores per cpu" FORCE)

  #
  # Check for hyper-threading - This is important for reserving threads for OpenMP tests...
  #
  math(EXPR MPI_MAX_NUMPROCS_PHYSICAL "${MPI_PHYSICAL_CORES} * ${MPI_CORES_PER_CPU}")
  if("${MPI_MAX_NUMPROCS_PHYSICAL}" STREQUAL "${MPIEXEC_MAX_NUMPROCS}")
    set(MPI_HYPERTHREADING
        "OFF"
        CACHE BOOL "Are we using hyper-threading?" FORCE)
  else()
    set(MPI_HYPERTHREADING
        "ON"
        CACHE BOOL "Are we using hyper-threading?" FORCE)
  endif()
endmacro()

# ------------------------------------------------------------------------------------------------ #
# Setup OpenMPI
# ------------------------------------------------------------------------------------------------ #
macro(setupOpenMPI)

  # sanity check, these OpenMPI flags (below) require version >= 1.8
  if(${MPI_VERSION} VERSION_LESS 1.8)
    message(FATAL_ERROR "OpenMPI version < 1.8 found.")
  endif()

  # Find cores/cpu, cpu/node, hyper-threading
  query_topology()

  # Extra options provided from the environment or by cmake
  if(DEFINED ENV{MPIEXEC_PREFLAGS})
    set(MPIEXEC_PREFLAGS "$ENV{MPIEXEC_PREFLAGS}")
  endif()

  if("${MPIEXEC_EXECUTABLE}" MATCHES "srun")
    set(preflags " --overlap") # -N 1 --cpu_bind=verbose,cores
    set(MPIEXEC_PREFLAGS ${preflags})
    set(MPIEXEC_PREFLAGS_PERFBENCH ${preflags})
    set(MPIEXEC_OMP_PREFLAGS "${MPIEXEC_PREFLAGS} -c ${MPI_CORES_PER_CPU}")
  else()
    # Notes:
    #
    # * For PERFBENCH that use Quo, we need '--map-by socket:SPAN' instead of '-bind-to none'.  The
    #   'bind-to none' is required to pack a node.
    # * Adding '--debug-daemons' is often requested by the OpenMPI dev team in conjunction with
    #   'export OMPI_MCA_btl_base_verbose=100' to obtain debug traces from openmpi.
    set(MPIEXEC_PREFLAGS_PERFBENCH "${MPIEXEC_PREFLAGS} --map-by socket:SPAN")
    if(NOT MPIEXEC_PREFLAGS MATCHES " -bind-to none")
      string(APPEND MPIEXEC_PREFLAGS " -bind-to none")
    endif()
    # Setup for OMP plus MPI
    if(NOT APPLE AND NOT MPIEXEC_OMP_PREFLAGS MATCHES "--map-by ppr")
      # -bind-to fails on OSX, See #691
      set(MPIEXEC_OMP_PREFLAGS
          "${MPIEXEC_PREFLAGS} --map-by ppr:${MPI_CORES_PER_CPU}:socket --report-bindings")
    endif()

    # Spectrum-MPI on darwin
    #
    # Limit communication to on-node via '-intra sm' or 'intra vader'
    # https://www.ibm.com/support/knowledgecenter/SSZTET_EOS/eos/guide_101.pdf
    if("${MPIEXEC_EXECUTABLE}" MATCHES "smpi" AND NOT MPIEXEC_PREFLAGS MATCHES "-intra sm")
      string(REPLACE "-bind-to none" "-bind-to core" MPIEXEC_PREFLAGS ${MPIEXEC_PREFLAGS})
      # string(REPLACE "-bind-to none" "-bind-to core" MPIEXEC_OMP_PREFLAGS ${MPIEXEC_OMP_PREFLAGS})
      set(smpi-sm-only "-intra sm -aff off --report-bindings")
      string(APPEND MPIEXEC_PREFLAGS " ${smpi-sm-only}")
      string(APPEND MPIEXEC_OMP_PREFLAGS " ${smpi-sm-only}")
      unset(smpi-sm-only)
    endif()
  endif()

  # Cache the result
  set(MPIEXEC_PREFLAGS
      "${MPIEXEC_PREFLAGS}"
      CACHE STRING "extra mpirun flags (list)." FORCE)
  set(MPIEXEC_PREFLAGS_PERFBENCH
      "${MPIEXEC_PREFLAGS_PERFBENCH}"
      CACHE STRING "extra mpirun flags (list)." FORCE)
  set(MPIEXEC_OMP_PREFLAGS
      "${MPIEXEC_OMP_PREFLAGS}"
      CACHE STRING "extra mpirun flags (list)." FORCE)

  mark_as_advanced(MPI_CPUS_PER_NODE MPI_CORES_PER_CPU MPI_PHYSICAL_CORES MPI_MAX_NUMPROCS_PHYSICAL
                   MPI_HYPERTHREADING)

endmacro()

# ------------------------------------------------------------------------------------------------ #
# Setup Mpich
# ------------------------------------------------------------------------------------------------ #
macro(setupMpichMPI)

  # Find cores/cpu, cpu/node, hyper-threading
  query_topology()

endmacro()

# ------------------------------------------------------------------------------------------------ #
# Setup Intel MPI
# ------------------------------------------------------------------------------------------------ #
macro(setupIntelMPI)

  # Find cores/cpu, cpu/node, hyper-threading
  query_topology()

endmacro()

# ------------------------------------------------------------------------------------------------ #
# Setup Cray MPI and wrappers
# ------------------------------------------------------------------------------------------------ #
macro(setupCrayMPI)

  query_topology()

  # salloc/sbatch options:
  #
  # * -N      limit job to a single node.
  # * --gres=craynetwork:0 This option allows more than one srun to be running at the same time on
  #   the Cray. There are 4 gres “tokens” available. If unspecified, each srun invocation will
  #   consume all of them. Setting the value to 0 means consume none and allow the user to run as
  #   many concurrent jobs as there are cores available on the node. This should only be specified
  #   on the salloc/sbatch command.  Gabe doesn't recommend this option for regression testing.
  # * --vm-overcommit=disable|enable Do not allow overcommit of heap resources.
  # * -p knl    Limit allocation to KNL nodes.
  #
  # srun options:
  #
  # * --cpu_bind=verbose,cores Bind MPI ranks to cores and print a summary of binding when run
  # * --exclusive This option will keep concurrent jobs from running on the same cores. If you want
  #   to background tasks to have them run simultaneously, this option is required to be set or they
  #   will stomp on the same cores.
  # * --overlap Allow steps to overlap each other on the CPUs. By default steps do not share CPUs
  #   with other parallel steps.

  set(preflags " ") # -N 1 --cpu_bind=verbose,cores
  string(APPEND preflags " --gres=craynetwork:0 --overlap")
  set(MPIEXEC_PREFLAGS
      ${preflags}
      CACHE STRING "extra mpirun flags (list)." FORCE)
  # consider adding '-m=cyclic'
  set(MPIEXEC_PREFLAGS_PERFBENCH
      ${preflags}
      CACHE STRING "extra mpirun flags (list)." FORCE)
  set(MPIEXEC_OMP_PREFLAGS
      "${MPIEXEC_PREFLAGS} -c ${MPI_CORES_PER_CPU}"
      CACHE STRING "extra mpirun flags (list)." FORCE)
endmacro()

# ------------------------------------------------------------------------------------------------ #
# Setup Spectrum MPI wrappers (Sierra, Rzansel) *
# https://www.ibm.com/support/knowledgecenter/SSZTET_EOS/eos/guide_101.pdf
# ------------------------------------------------------------------------------------------------ #
macro(setupSpectrumMPI)

  # Find cores/cpu, cpu/node, hyper-threading
  query_topology()

  # jsrun options
  #
  # Ref: https://hpc.llnl.gov/training/tutorials/using-lcs-sierra-system#jsrun
  #
  # Basic: -n Total number of cores
  #
  # Using resource sets: -a Number of tasks per resource set -n Total number of resource sets for
  # the job -c Number of CPUs (cores) per resource set -g Number of GPUs per resource set -r Number
  # of resource sets per node.
  #
  # Examples: - jsrun -a4 -c1 -g1  => 4 tasks on 4 cores that share 1 gpu. - jsrun -r4 -n8      => 2
  # nodes, 4 resource sets per node, 8 resource sets total - jsrun -a4 -c16 -g2 => 4 tasks, 16
  # cores, 2 gpus

  if(MPIEXEC_EXECUTABLE MATCHES jsrun)
    # 1 resource set; 1 thread; no gpu; no binding --nrs 1
    set(MPIEXEC_PREFLAGS "-c 1 -g 0 --bind none")
  elseif(MPIEXEC_EXECUTABLE MATCHES lrun)
    set(MPIEXEC_PREFLAGS "--pack --threads=1 -v") # --bind=off
  else()
    message(FATAL_ERROR "Unexpected mpirun: ${MPIEXEC_EXECUTABLE}")
  endif()
  # --pack ==> -c 1 -g 0.  This is actually bad for us. Disable lrun -n 2 -c 10 --threads=10
  # --bind=off ==> jsrun --np 1 --nrs 1 -c ALL_CPUS -g ALL_GPUS -d plane:1 -b rs -X 1 consider:
  # jsrun --np 2 --nrs 1 -c 10 -g 0 -bind none

  #
  # Setup for OMP plus MPI
  #

  if(DEFINED ENV{OMP_NUM_THREADS})
    set(ont $ENV{OMP_NUM_THREADS})
  else()
    set(ont 4)
  endif()
  if(MPIEXEC_EXECUTABLE MATCHES jsrun)
    # 1 resource set; OMP_NUM_THREADS tasks; no gpu; no binding
    set(MPIEXEC_OMP_PREFLAGS "--nrs 1 -c ${ont} -g 0 --bind none")
  elseif(MPIEXEC_EXECUTABLE MATCHES lrun)
    set(MPIEXEC_OMP_PREFLAGS "--threads=${ont} -c${ont} --bind=none -v")
    # --pack ? --smt=4 ?
  else()
    message(FATAL_ERROR "Unexpected mpirun: ${MPIEXEC_EXECUTABLE}")
  endif()

  set(MPIEXEC_OMP_PREFLAGS
      ${MPIEXEC_OMP_PREFLAGS}
      CACHE STRING "extra mpirun flags (list)." FORCE)

  mark_as_advanced(MPI_CPUS_PER_NODE MPI_CORES_PER_CPU MPI_PHYSICAL_CORES MPI_MAX_NUMPROCS_PHYSICAL
                   MPI_HYPERTHREADING)

endmacro()

# ------------------------------------------------------------------------------------------------ #
# Setup MPI when on Linux
# ------------------------------------------------------------------------------------------------ #
macro(setupMPILibrariesUnix)

  # MPI ---------------------------------------------------------------------
  if(NOT "${DRACO_C4}" STREQUAL "SCALAR")

    message(STATUS "Looking for MPI...")

    # Preserve data that may already be set.
    if(DEFINED ENV{MPIRUN})
      set(MPIEXEC_EXECUTABLE
          $ENV{MPIRUN}
          CACHE STRING "Program to execute MPI parallel programs.")
    elseif(DEFINED ENV{MPIEXEC_EXECUTABLE})
      set(MPIEXEC_EXECUTABLE
          $ENV{MPIEXEC_EXECUTABLE}
          CACHE STRING "Program to execute MPI parallel programs.")
    elseif(DEFINED ENV{MPIEXEC})
      set(MPIEXEC_EXECUTABLE
          $ENV{MPIEXEC}
          CACHE STRING "Program to execute MPI parallel programs.")
    endif()

    # If this is a Cray system and the Cray MPI compile wrappers are used, or if this is CTS-1 with
    # Toss3, then do some special setup:

    if(CMAKE_CXX_COMPILER_WRAPPER MATCHES CrayPrgEnv OR IS_DIRECTORY "/usr/projects/hpcsoft/toss3/")
      if(NOT EXISTS ${MPIEXEC_EXECUTABLE})
        find_program(MPIEXEC_EXECUTABLE srun)
      endif()
      set(MPIEXEC_EXECUTABLE
          ${MPIEXEC_EXECUTABLE}
          CACHE STRING "Program to execute MPI parallel programs." FORCE)
      set(MPIEXEC_NUMPROC_FLAG
          "-n"
          CACHE STRING "mpirun flag used to specify the number of processors to use")
    elseif(DEFINED ENV{SYS_TYPE} AND "$ENV{SYS_TYPE}" MATCHES "ppc64le_ib_p9") # ATS-2
      if(NOT EXISTS ${MPIEXEC_EXECUTABLE})
        find_program(MPIEXEC_EXECUTABLE lrun)
      endif()
      set(MPIEXEC_EXECUTABLE
          ${MPIEXEC_EXECUTABLE}
          CACHE STRING "Program to execute MPI parallel programs." FORCE)
      set(MPIEXEC_NUMPROC_FLAG
          "--np"
          CACHE STRING "mpirun flag used to specify the number of processors to use" FORCE)
    endif()

    # Call the standard CMake FindMPI macro.
    find_package(MPI QUIET)

    # Try to discover the MPI flavor and the vendor version. Returns MPI_VERSION, MPI_FLAVOR as
    # cache variables
    setmpiflavorver()

    # Set additional flags, environments that are MPI vendor specific.
    if("${MPI_FLAVOR}" MATCHES "openmpi")
      setupopenmpi()
    elseif("${MPI_FLAVOR}" MATCHES "mpich")
      setupmpichmpi()
    elseif("${MPI_FLAVOR}" MATCHES "intel")
      setupintelmpi()
    elseif("${MPI_FLAVOR}" MATCHES "spectrum")
      setupspectrummpi()
    elseif("${MPI_FLAVOR}" MATCHES "cray")
      setupcraympi()
    else()
      message(
        FATAL_ERROR
          "
The Draco build system doesn't know how to configure the build for
  MPIEXEC_EXECUTABLE         = ${MPIEXEC_EXECUTABLE}
  DBS_MPI_VER                = ${DBS_MPI_VER}
  CMAKE_CXX_COMPILER_WRAPPER = ${CMAKE_CXX_COMPILER_WRAPPER}")
    endif()

    # Mark some of the variables created by the above logic as 'advanced' so that they do not show
    # up in the 'simple' ccmake view.
    mark_as_advanced(MPI_EXTRA_LIBRARY MPI_LIBRARY)

    message(STATUS "Looking for MPI.......found ${MPIEXEC_EXECUTABLE}")

    # Sanity Checks for DRACO_C4==MPI
    if("${MPI_CORES_PER_CPU}x" STREQUAL "x")
      message(FATAL_ERROR "setupMPILibrariesUnix:: MPI_CORES_PER_CPU " "is not set!")
    endif()

  endif()

  # Set DRACO_C4 and other variables. Returns DRACO_C4, C4_SCALAR, C4_MPI
  setupdracompivars()

  set_package_properties(
    MPI PROPERTIES
    URL "http://www.open-mpi.org/"
    DESCRIPTION "A High Performance Message Passing Library"
    TYPE RECOMMENDED
    PURPOSE "If not available, all Draco components will be built as scalar applications.")

  mark_as_advanced(MPIEXEC_OMP_PREFLAGS MPI_LIBRARIES)

endmacro()

# ------------------------------------------------------------------------------------------------ #
# setupMPILibrariesWindows
# ------------------------------------------------------------------------------------------------ #
macro(setupMPILibrariesWindows)

  set(verbose FALSE)

  # MPI ---------------------------------------------------------------------
  if(NOT "${DRACO_C4}" STREQUAL "SCALAR")

    message(STATUS "Looking for MPI...")
    set(MPI_ASSUME_NO_BUILTIN_MPI TRUE)
    set(MPI_GUESS_LIBRARY_NAME "MSMPI")
    find_package(MPI QUIET)

    if(EXISTS "$ENV{MSMPI_INC}")
      # if msmpi is installed via vcpkg, then use the vcpkg include path instead of the system one.
      file(TO_CMAKE_PATH $ENV{MSMPI_INC} MSMPI_INC)
      unset(tmp)
      foreach(ipath ${MPI_C_INCLUDE_DIRS})
        if("${ipath}" STREQUAL "${MSMPI_INC}")
          continue()
        endif()
        list(APPEND tmp ${ipath})
      endforeach()
      set(MPI_C_INCLUDE_DIRS ${tmp})
      unset(tmp)
    endif()

    if(verbose)
      message(
        "
        MPI_C_FOUND        = ${MPI_C_FOUND}
        MPIEXEC_EXECUTABLE = ${MPIEXEC_EXECUTABLE}
        MPI_C_LIBRARIES    = ${MPI_C_LIBRARIES}
        MPI_CXX_LIBRARIES  = ${MPI_CXX_LIBRARIES}
        MPI_Fortran_LIBRARIES = ${MPI_Fortran_LIBRARIES}
        MPI_C_LIB_NAMES    = ${MPI_C_LIB_NAMES}
        CMAKE_GENERATOR    = ${CMAKE_GENERATOR}
        MPI_C_INCLUDE_DIRS = ${MPI_C_INCLUDE_DIRS}
        MPI_C_VERSION_MAJOR = ${MPI_C_VERSION_MAJOR}")
    endif()

    # If this macro is called from a MinGW builds system (for a CAFS subdirectory) and is trying to
    # discover MS-MPI, the above check will fail (when CMake > 3.12). However, MS-MPI is known to be
    # good when linking with Visual Studio so override the 'failed' report.
    if("${MPI_C_LIBRARIES}" MATCHES "msmpi" AND "${CMAKE_GENERATOR}" STREQUAL "MinGW Makefiles")
      if(EXISTS "${MPI_C_LIBRARIES}" AND EXISTS "${MPI_C_INCLUDE_DIRS}")
        set(MPI_C_FOUND TRUE)
        set(MPI_Fortran_FOUND TRUE)
      endif()
    endif()

    if(verbose)
      message("
        MPI_C_FOUND       = ${MPI_C_FOUND}
        MPI_Fortran_FOUND = ${MPI_Fortran_FOUND}")
    endif()

    # For MS-MPI, mpifptr.h is architecture dependent. Figure out what arch this is and save this
    # path to MPI_Fortran_INCLUDE_PATH
    list(GET MPI_C_LIBRARIES 0 first_c_mpi_library)
    if(first_c_mpi_library AND NOT MPI_Fortran_INCLUDE_PATH)
      get_filename_component(MPI_Fortran_INCLUDE_PATH "${first_c_mpi_library}" DIRECTORY)
      if(EXISTS "${MPI_Fortran_INCLUDE_PATH}/../include/mpifptr.h")
        get_filename_component(MPI_Fortran_INCLUDE_PATH "${MPI_Fortran_INCLUDE_PATH}/../include"
                               REALPATH)
      elseif(EXISTS "${MPI_Fortran_INCLUDE_PATH}/../../include/mpifptr.h")
        get_filename_component(MPI_Fortran_INCLUDE_PATH "${MPI_Fortran_INCLUDE_PATH}/../../include"
                               REALPATH)
      else()
        string(REGEX REPLACE "[Ll]ib" "Include" MPI_Fortran_INCLUDE_PATH
                             ${MPI_Fortran_INCLUDE_PATH})
        set(MPI_Fortran_INCLUDE_PATH
            "${MPI_C_INCLUDE_PATH};${MPI_Fortran_INCLUDE_PATH}"
            CACHE STRING "Location for MPI include files for Fortran.")
      endif()
      if(verbose)
        message("    MPI_Fortran_INCLUDE_PATH = ${MPI_Fortran_INCLUDE_PATH}")
      endif()
    endif()

    setupdracompivars()

    # Find the version. [This is not working (hardwire it for now)]
    execute_process(
      COMMAND "${MPIEXEC_EXECUTABLE}" -help
      OUTPUT_VARIABLE DBS_MPI_VER_OUT
      ERROR_VARIABLE DBS_MPI_VER_ERR ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
                     ERROR_STRIP_TRAILING_WHITESPACE)
    if("${DBS_MPI_VER_OUT}" MATCHES "Microsoft MPI Startup Program")
      string(REGEX REPLACE ".*Version ([0-9.]+).*" "\\1" DBS_MPI_VER
                           "${DBS_MPI_VER_OUT}${DBS_MPI_VER_ERR}")
      string(REPLACE "." ";" DBS_MPI_VER ${DBS_MPI_VER})
      list(GET DBS_MPI_VER 0 DBS_MPI_VER_MAJOR)
      list(GET DBS_MPI_VER 1 DBS_MPI_VER_MINOR)
      set(DBS_MPI_VER "${DBS_MPI_VER_MAJOR}.${DBS_MPI_VER_MINOR}")
    else()
      set(DBS_MPI_VER "5.0")
    endif()
    if(verbose)
      message("    DBS_MPI_VER = ${DBS_MPI_VER}")
    endif()

    set_package_properties(
      MPI PROPERTIES
      URL "https://msdn.microsoft.com/en-us/library/bb524831%28v=vs.85%29.aspx"
      DESCRIPTION "Microsoft MPI"
      TYPE RECOMMENDED
      PURPOSE "If not available, all Draco components will be built as scalar applications.")

    # Check flavor and add optional flags
    if("${MPIEXEC_EXECUTABLE}" MATCHES "Microsoft MPI")
      set(MPI_FLAVOR
          "MicrosoftMPI"
          CACHE STRING "Flavor of MPI.")

      # Use wmic to learn about the current machine
      execute_process(COMMAND wmic cpu get NumberOfCores
                      OUTPUT_VARIABLE MPI_CORES_PER_CPU OUTPUT_STRIP_TRAILING_WHITESPACE)
      execute_process(COMMAND wmic computersystem get NumberOfLogicalProcessors
                      OUTPUT_VARIABLE MPIEXEC_MAX_NUMPROCS OUTPUT_STRIP_TRAILING_WHITESPACE)
      execute_process(COMMAND wmic computersystem get NumberOfProcessors
                      OUTPUT_VARIABLE MPI_CPUS_PER_NODE OUTPUT_STRIP_TRAILING_WHITESPACE)
      string(REGEX REPLACE ".*[\n]([0-9]+$)" "\\1" MPI_CORES_PER_CPU ${MPI_CORES_PER_CPU})
      string(REGEX REPLACE ".*[\n]([0-9]+$)" "\\1" MPIEXEC_MAX_NUMPROCS ${MPIEXEC_MAX_NUMPROCS})
      string(REGEX REPLACE ".*[\n]([0-9]+$)" "\\1" MPI_CPUS_PER_NODE ${MPI_CPUS_PER_NODE})

      set(MPI_CPUS_PER_NODE
          ${MPI_CPUS_PER_NODE}
          CACHE STRING "Number of multi-core CPUs per node" FORCE)
      set(MPI_CORES_PER_CPU
          ${MPI_CORES_PER_CPU}
          CACHE STRING "Number of cores per CPU" FORCE)
      set(MPIEXEC_MAX_NUMPROCS
          ${MPIEXEC_MAX_NUMPROCS}
          CACHE STRING "Total number of available MPI ranks" FORCE)

      # Check for hyper-threading - This is important for reserving threads for OpenMP tests...

      math(EXPR MPI_MAX_NUMPROCS_PHYSICAL "${MPI_CPUS_PER_NODE} * ${MPI_CORES_PER_CPU}")
      if("${MPI_MAX_NUMPROCS_PHYSICAL}" STREQUAL "${MPIEXEC_MAX_NUMPROCS}")
        set(MPI_HYPERTHREADING
            "OFF"
            CACHE BOOL "Are we using hyper-threading?" FORCE)
      else()
        set(MPI_HYPERTHREADING
            "ON"
            CACHE BOOL "Are we using hyper-threading?" FORCE)
      endif()

      set(MPIEXEC_OMP_PREFLAGS
          "-exitcodes"
          CACHE STRING "extra mpirun flags (list)." FORCE)
    endif()

  else()
    # Set DRACO_C4 and other variables
    setupdracompivars()
  endif() # NOT "${DRACO_C4}" STREQUAL "SCALAR"

  # Found MPI_C, but not MPI_CXX -- create a duplicate to satisfy link targets.
  if(TARGET MPI::MPI_C AND NOT TARGET MPI::MPI_CXX)

    if(verbose)
      message(
        "    Found target MPI::MPI_C but not MPI::MPI_CXX

      set_target_properties(MPI::MPI_CXX PROPERTIES
        IMPORTED_LOCATION_RELEASE         \"${MPI_C_LIBRARIES}\"
        IMPORTED_IMPLIB                   \"${MPI_C_LIBRARIES}\"
        INTERFACE_INCLUDE_DIRECTORIES     \"${MPI_C_INCLUDE_DIRS}\"
        IMPORTED_CONFIGURATIONS           Release
        IMPORTED_LINK_INTERFACE_LANGUAGES \"CXX\" )
      ")
    endif()

    # Windows systems with dll libraries.
    add_library(MPI::MPI_CXX SHARED IMPORTED)

    # Windows with dlls, but only Release libraries.
    set_target_properties(
      MPI::MPI_CXX
      PROPERTIES
      IMPORTED_LOCATION_RELEASE "${MPI_C_LIBRARIES}"
      IMPORTED_IMPLIB "${MPI_C_LIBRARIES}"
      INTERFACE_INCLUDE_DIRECTORIES "${MPI_C_INCLUDE_DIRS}"
      IMPORTED_CONFIGURATIONS Release
      IMPORTED_LINK_INTERFACE_LANGUAGES "CXX")

  endif()

  # Don't link to the C++ MS-MPI library when compiling with MinGW gfortran.  Instead, link to
  # libmsmpi.a that was created via gendef.exe and dlltool.exe from msmpi.dll.  Ref:
  # http://www.geuz.org/pipermail/getdp/2012/001520.html, or
  # https://github.com/KineticTheory/Linux-HPC-Env/wiki/Setup-Win32-development-environment

  # Preparing Microsoft's MPI to work with x86_64-w64-mingw32-gfortran by creating libmsmpi.a. (Last
  # tested: 2017-08-31)
  #
  # 1. You need MinGW/MSYS. Please make sure that the Devel tools are installed.
  # 2. Download and install Microsoft's MPI. You need the main program and the SDK.
  # 3. In the file %MSMPI_INC%\mpif.h, replace INT_PTR_KIND() by 8
  # 4. Create a MSYS version of the MPI library:
  # ~~~
  #     cd %TEMP%
  #     copy c:\Windows\System32\msmpi.dll msmpi.dll
  #     gendef msmpi.dll
  #     dlltool -d msmpi.def -l libmsmpi.a -D msmpi.dll
  #     del msmpi.def
  #     copy libmsmpi.a %MSMPI_LIB32%/libmsmpi.a
  # ~~~

  if(WIN32
     AND DEFINED CMAKE_Fortran_COMPILER
     AND TARGET MPI::MPI_Fortran)

    if(verbose)
      message("    Win32 AND CMAKE_Fortran_COMPILER AND TARGET MPI::MPI_Fortran
      ")
    endif()

    # only do this if we are in a CMakeAddFortranSubdirectory directive when the main Generator is
    # Visual Studio and the Fortran subdirectory uses gfortran with Makefiles.

    # MS-MPI has an architecture specific include directory that FindMPI.cmake doesn't seem to
    # pickup correctly.  Add it here.
    get_target_property(mpilibdir MPI::MPI_Fortran INTERFACE_LINK_LIBRARIES)
    get_target_property(mpiincdir MPI::MPI_Fortran INTERFACE_INCLUDE_DIRECTORIES)
    foreach(arch x86 x64)
      string(FIND "${mpilibdir}" "lib/${arch}" found_lib_arch)
      string(FIND "${mpiincdir}" "include/${arch}" found_inc_arch)
      if(${found_lib_arch} GREATER 0 AND ${found_inc_arch} LESS 0)
        if(IS_DIRECTORY "${mpiincdir}/${arch}")
          list(APPEND mpiincdir "${mpiincdir}/${arch}")
        endif()
      endif()
    endforeach()

    # Reset the include directories for MPI::MPI_Fortran to pull in the extra $arch locations (if
    # any)

    list(REMOVE_DUPLICATES mpiincdir)
    set_target_properties(MPI::MPI_Fortran PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${mpiincdir}")

    if(verbose)
      message("set_target_properties(MPI::MPI_Fortran
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES \"${mpiincdir}\")")
    endif()

  endif()

  if(${MPI_C_FOUND})
    message(STATUS "Looking for MPI.......found ${MPIEXEC_EXECUTABLE}")
  else()
    message(STATUS "Looking for MPI.......not found")
  endif()

  mark_as_advanced(MPI_FLAVOR MPIEXEC_OMP_PREFLAGS MPI_LIBRARIES)

endmacro(setupMPILibrariesWindows)

# ------------------------------------------------------------------------------------------------ #
# Helper
# ------------------------------------------------------------------------------------------------ #
macro(setupMPILibraries)
  if(UNIX)
    setupmpilibrariesunix()
  elseif(WIN32)
    setupmpilibrarieswindows()
  endif()
endmacro()

# ------------------------------------------------------------------------------------------------ #
# End of setupMPI.cmake
# ------------------------------------------------------------------------------------------------ #
