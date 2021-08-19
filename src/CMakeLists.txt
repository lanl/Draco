#--------------------------------------------*-cmake-*---------------------------------------------#
# file   src/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 April 28
# brief  Instructions for building src level Makefile.
# note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.17.0)
project(draco_src_dir)

# Provide helper functions used by component CMakeLists.txt files
include(component_macros)

# Extra 'draco-only' flags
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 6.0)
  toggle_compiler_flag(TRUE "-Wdouble-promotion" "CXX" "DEBUG")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  string(REPLACE "/W2" "/W4" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  string(REPLACE "/W2" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif()

# <experimental/mdspan> doesn't work with some tool chains...
option(USE_MDSPAN "Use experimental/mdspan" ON)
if((MSVC AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
   OR (CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.0))
  set(USE_MDSPAN OFF)
endif()

# Level 1
message(" ")
message(STATUS "Configuring Level 1 packages --")
add_dir_if_exists(ds++)

# Only test selected components if CI_CLANG_TIDY is running

# Level 2
message(" ")
message(STATUS "Configuring Level 2 packages --")
add_dir_if_exists(c4) #< needs ds++
add_dir_if_exists(cdi) #< needs ds++
add_dir_if_exists(device) #< needs ds++
add_dir_if_exists(FortranChecks) #< needs ds++
add_dir_if_exists(lapack_wrap) #< needs ds++
add_dir_if_exists(linear) #< needs ds++
add_dir_if_exists(memory) #< needs ds++
add_dir_if_exists(mesh_element) #< needs ds++
add_dir_if_exists(ode) #< needs ds++
add_dir_if_exists(units) #< needs ds++

# Level 3
message(" ")
message(STATUS "Configuring Level 3 packages --")
if(USE_MDSPAN)
  add_dir_if_exists(experimental) #< only the tests need ds++
  add_dir_if_exists(cdi_CPEloss) #< needs cdi, units, experimental
endif()
add_dir_if_exists(cdi_eospac) #< needs parser, ode, cdi
add_dir_if_exists(cdi_ipcress) #< needs cdi
add_dir_if_exists(compton_interface) #< needs c4, CSK
add_dir_if_exists(diagnostics) #< needs c4
add_dir_if_exists(fit) #< needs linear
add_dir_if_exists(kde) #< needs c4
add_dir_if_exists(meshReaders) #< needs c4
add_dir_if_exists(min) #< needs linear
add_dir_if_exists(norms) #< needs c4
add_dir_if_exists(parser) #< needs units, c4, mesh_element
add_dir_if_exists(rng) #< needs ds++, device
add_dir_if_exists(roots) #< needs linear
add_dir_if_exists(VendorChecks) #< needs c4
add_dir_if_exists(viz) #< needs c4

# Level 4
message(" ")
message(STATUS "Configuring Level 4 packages --")
add_dir_if_exists(cdi_analytic) #< needs parser, roots, ode, cdi
add_dir_if_exists(cdi_ndi) #< needs ds++, rng, cdi
add_dir_if_exists(RTT_Format_Reader) #< needs meshReaders
add_dir_if_exists(special_functions) #< needs roots, units, GSL

# Level 5
message(" ")
message(STATUS "Configuring Level 5 packages --")
add_dir_if_exists(compton_tools) #< needs ds++, units, cdi_analtyic
add_dir_if_exists(mesh) #< needs c4, ds++, mesh_element
add_dir_if_exists(quadrature) #< needs mesh_element, parser, special_functions

# Summary

if(VERBOSE OR ENV{VERBOSE})
  message(" ")
  feature_summary(WHAT ALL INCLUDE_QUIET_PACKAGES FATAL_ON_MISSING_REQUIRED_PACKAGES QUIET_ON_EMPTY)
endif()

message(
  "
Draco build summary:

Version     : Draco ${Draco_VERSION_FULL}
Build type  : ${Draco_BUILD_TYPE}
            : DRACO_DBC_LEVEL   = ${DRACO_DBC_LEVEL}
            : DRACO_DIAGNOSTICS = ${DRACO_DIAGNOSTICS}
            : DRACO_TIMING      = ${DRACO_TIMING}
C4 Model    : ${DRACO_C4}
Prefix dir  : ${CMAKE_INSTALL_PREFIX}
Source dir  : ${Draco_SOURCE_DIR}
Build dir   : ${Draco_BINARY_DIR}

CXX Compiler: ${CMAKE_CXX_COMPILER}")
include(string_manip)
if(CMAKE_CONFIGURATION_TYPES)
  hanging_indent(100 5 "CXX FLAGS: ${CMAKE_CXX_FLAGS} ")
  hanging_indent(100 5 "CXX Debug FLAGS: ${CMAKE_CXX_FLAGS_DEBUG}")
  hanging_indent(100 5 "CXX Release FLAGS: ${CMAKE_CXX_FLAGS_RELEASE}")
  if(_LANGUAGES_ MATCHES Fortran)
    hanging_indent(100 5 "Fortran: ${CMAKE_Fortran_COMPILER}")
    hanging_indent(100 5 "FC FLAGS: ${CMAKE_Fortran_FLAGS}")
    hanging_indent(100 5 "FC Debug FLAGS: ${CMAKE_Fortran_FLAGS_DEBUG}")
    hanging_indent(100 5 "FC Release FLAGS: ${CMAKE_Fortran_FLAGS_RELEASE}")
  endif()
  if(_LANGUAGES_ MATCHES CUDA)
    message("CUDA        : ${CMAKE_CUDA_COMPILER}")
    message("Cuda FLAGS  : ${CMAKE_CUDA_FLAGS}")
    message("Cuda Debug F: ${CMAKE_CUDA_FLAGS_DEBUG}")
    message("Cuda Release: ${CMAKE_CUDA_FLAGS_RELEASE}")
  endif()
else()
  hanging_indent(100 5 "C FLAGS: ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${Draco_BUILD_TYPE}}")
  hanging_indent(100 5 "CXX FLAGS: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${Draco_BUILD_TYPE}}")
  if(_LANGUAGES_ MATCHES Fortran)
    message("Fortran: ${CMAKE_Fortran_COMPILER}")
    hanging_indent(100 5
                   "FC FLAGS: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${Draco_BUILD_TYPE}}")
  endif()
  if(_LANGUAGES_ MATCHES CUDA)
    message("CUDA        : ${CMAKE_CUDA_COMPILER} (${CMAKE_CUDA_COMPILER_ID})")
    message("Cuda FLAGS  : ${CMAKE_CUDA_FLAGS} ${CMAKE_CUDA_FLAGS_${Draco_BUILD_TYPE}}")
  endif()
endif()
if(CAFS_Fortran_COMPILER)
  message("CAFS Fortran: ${CAFS_Fortran_COMPILER}")
endif()
if("${DRACO_C4}" STREQUAL "MPI")
  hanging_indent(
    100 5 "mpirun command: ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} N ${MPIEXEC_PREFLAGS}")
endif()
message("Library Type: ${DRACO_LIBRARY_TYPE}
")

if(CMAKE_CXX_COMPILER_WRAPPER STREQUAL CrayPrgEnv)
  message("Cray system detected: CC -craype-verbose -V|--version:
")
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Cray" AND NOT DEFINED ENV{CRAY_PE_USE_CLANG})
    set(ver_opt "-V")
  else()
    set(ver_opt "--version")
  endif()
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -craype-verbose ${ver_opt} OUTPUT_VARIABLE tmp)
  message("${tmp}")
endif()

# Push some variables up one level
set(Draco_EXPORT_TARGET_PROPERTIES
    "${Draco_EXPORT_TARGET_PROPERTIES}"
    PARENT_SCOPE)

# ------------------------------------------------------------------------------------------------ #
# Debug build system state:
# ------------------------------------------------------------------------------------------------ #
option(DBS_PRINT_DEBUG_INFO "Print extra build system data" OFF)
if(DBS_PRINT_DEBUG_INFO)
  include(debug_macros)
  # Debug variables
  message("CMake variables and values:")
  echo_all_cmake_variable_values()

  # Debug target properties (run 'make help' to see a list of targets)
  message("CMake targets and properties:")
  set(dbs_targets # Lib_dsxx
      # Lib_dsxx_test
      Exe_Ipcress_Interpreter)
  echo_targets(${dbs_targets})
endif()

# ------------------------------------------------------------------------------------------------ #
# End of src/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
