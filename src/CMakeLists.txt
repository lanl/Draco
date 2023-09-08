#--------------------------------------------*-cmake-*---------------------------------------------#
# file   src/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 April 28
# brief  Instructions for building src level Makefile.
# note   Copyright (C) 2010-2023 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
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
force_compiler_flags_to_cache("C;CXX")

# <experimental/mdspan> doesn't work with some tool chains...
option(USE_MDSPAN "Use experimental/mdspan" ON)
if((MSVC AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
   OR (CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.0))
  set(USE_MDSPAN OFF)
endif()

# Level 1
message(" ")
message(STATUS "Configuring Level 1 packages --")
add_subdirectory(ds++)

# Only test selected components if CI_CLANG_TIDY is running

# Level 2
message(" ")
message(STATUS "Configuring Level 2 packages --")
add_subdirectory(c4) #< needs ds++
add_subdirectory(cdi) #< needs ds++
add_subdirectory(device) #< needs ds++
add_subdirectory(FortranChecks) #< needs ds++
add_subdirectory(linear) #< needs ds++
add_subdirectory(memory) #< needs ds++
add_subdirectory(mesh_element) #< needs ds++
add_subdirectory(ode) #< needs ds++
add_subdirectory(units) #< needs ds++

# Level 3
message(" ")
message(STATUS "Configuring Level 3 packages --")
if(USE_MDSPAN)
  add_subdirectory(experimental) #< only the tests need ds++
  add_subdirectory(cdi_CPEloss) #< needs cdi, units, experimental
endif()
add_subdirectory(cdi_eospac) #< needs parser, ode, cdi
add_subdirectory(cdi_ipcress) #< needs cdi
add_subdirectory(diagnostics) #< needs c4
add_subdirectory(fit) #< needs linear
add_subdirectory(kde) #< needs c4
add_subdirectory(meshReaders) #< needs c4
add_subdirectory(min) #< needs linear
add_subdirectory(norms) #< needs c4
add_subdirectory(parser) #< needs units, c4, mesh_element
if(RANDOM123_FOUND)
  add_subdirectory(rng) #< needs ds++, device
endif()
add_subdirectory(roots) #< needs linear
add_subdirectory(VendorChecks) #< needs c4
add_subdirectory(viz) #< needs c4
add_subdirectory(predict)

# Level 4
message(" ")
message(STATUS "Configuring Level 4 packages --")
add_subdirectory(cdi_analytic) #< needs parser, roots, ode, cdi
add_subdirectory(cdi_ndi) #< needs ds++, rng, cdi
add_subdirectory(RTT_Format_Reader) #< needs meshReaders
add_subdirectory(special_functions) #< needs roots, units, GSL

# Level 5
message(" ")
message(STATUS "Configuring Level 5 packages --")
add_subdirectory(compton_tools) #< needs ds++, units, cdi_analtyic
add_subdirectory(mesh) #< needs c4, ds++, mesh_element
add_subdirectory(quadrature) #< needs mesh_element, parser, special_functions

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
  if(_LANGUAGES_ MATCHES HIP)
    message("HIP         : ${CMAKE_HIP_COMPILER}")
    message("HIP FLAGS   : ${CMAKE_HIP_FLAGS}")
    message("HIP Debug FL: ${CMAKE_HIP_FLAGS_DEBUG}")
    message("HIP Release : ${CMAKE_HIP_FLAGS_RELEASE}")
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
  if(_LANGUAGES_ MATCHES HIP)
    message("HIP         : ${CMAKE_HIP_COMPILER} (${CMAKE_HIP_COMPILER_ID})")
    message("HIP FLAGS   : ${CMAKE_HIP_FLAGS} ${CMAKE_HIP_FLAGS_${Draco_BUILD_TYPE}}")
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
