#--------------------------------------------*-cmake-*---------------------------------------------#
# file   dsxx/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 Sep 24
# brief  Generate build project files for dsxx.
# note   Copyright (C) 2024-2026 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(dsxx CXX ${GPU_DBS_STRING})

# ------------------------------------------------------------------------------------------------ #
# Generate config.h and device_config.h (only occurs when cmake is run)

# ------------------------------------------------------------------------------------------------ #
# GPU-specific configuration file
# ------------------------------------------------------------------------------------------------ #
if(DEFINED GPU_DBS_STRING)
  if(GPU_DBS_STRING STREQUAL "HIP")
    set(USE_HIP ON)
  elseif(GPU_DBS_STRING STREQUAL "CUDA")
    set(USE_CUDA ON)
  endif()
endif()
set(TEST_KERNEL_BINDIR
    ${PROJECT_BINARY_DIR}/test
    CACHE PATH "GPU kernel binary install location")
configure_file(device_config.h.in ${PROJECT_BINARY_DIR}/dsxx/device_config.h)

# ------------------------------------------------------------------------------------------------ #
# Build GPU_Device
# ------------------------------------------------------------------------------------------------ #
if(HAVE_GPU AND DRACO_USE_GPU)
  set(include_dirs PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>)
  if(EXISTS "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
    string(APPEND include_dirs ";$<BUILD_INTERFACE:${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}>")
    string(APPEND include_dirs ";$<INSTALL_INTERFACE:${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}>")
  elseif(EXISTS "$ENV{ROCM_PATH}")
    string(APPEND include_dirs ";$<BUILD_INTERFACE:$ENV{ROCM_PATH}/include>")
    string(APPEND include_dirs ";$<INSTALL_INTERFACE:$ENV{ROCM_PATH}/include>")
  endif()
endif()
# ------------------------------------------------------------------------------------------------ #

# Review config/platform_checks.cmake for many of these settings.

# Choose platforms for fpe_trap
if(${CMAKE_SYSTEM} MATCHES "OSF")
  set(FPETRAP_OSF_ALPHA 1)
elseif(${CMAKE_SYSTEM} MATCHES "Linux")
  # If the OS is Linux, fpe_trap should work as long as feenableexcept is available.  That is, it
  # should work for both GNU and Intel compiler suites.
  set(FPETRAP_LINUX_X86 1)

  # Turn fpe_trap off for PGI compilers.
  if(${CMAKE_C_COMPILER} MATCHES pgcc)
    set(FPETRAP_LINUX_X86 0)
  endif()
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  if(${_CMAKE_OSX_MACHINE} MATCHES "ppc")
    set(FPETRAP_DARWIN_PPC 1)
  else()
    set(FPETRAP_DARWIN_INTEL 1)
  endif()

elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set(FPETRAP_WINDOWS_X86 1)
endif()

# Choose settings for DRACO_DIAGNOSTICS. Uses a bit mask:
#
# * 001 - basic print messages (negligible performance hit)
# * 010 - detailed print messages (some performance hit or intrusive output)
# * 100 - additional computation required for diagnostic messages (including fpe_trap)
set(DRACO_DIAGNOSTICS
    "0"
    CACHE STRING "Diagnostics verbosity level (0-7). Default is 0.")
# provide a constrained drop down menu for cmake-gui
# cmake-format: off
set_property(CACHE DRACO_DIAGNOSTICS PROPERTY STRINGS 0 1 2 3 4 5 6 7)
# cmake-format: on

# Provide an enum from CMAKE_CXX_COMPILER_VERSION
if(NOT "${CMAKE_CXX_COMPILER_VERSION}x" STREQUAL "x")
  string(REPLACE "." "" CMAKE_CXX_COMPILER_VERSION_ENUM ${CMAKE_CXX_COMPILER_VERSION})
endif()

if(CMAKE_CXX_COMPILER_WRAPPER STREQUAL CrayPrgEnv)
  # Is this CCE
  if(DEFINED ENV{CRAY_PE_CCE_VARIANT} OR MPI_CXX_LIB_NAMES MATCHES "mpi_cray")
    set(CRAYPE_CCE TRUE)
  endif()
endif()

# Save the configure date
string(TIMESTAMP Draco_BUILD_DATE "%Y %b %d")

# Information about the current machine
cmake_host_system_information(RESULT COMPILED_HOSTNAME QUERY HOSTNAME)
cmake_host_system_information(RESULT COMPILED_FQDN QUERY FQDN)
cmake_host_system_information(RESULT COMPILED_OS_NAME QUERY OS_NAME)
cmake_host_system_information(RESULT COMPILED_OS_RELEASE QUERY OS_RELEASE)
cmake_host_system_information(RESULT COMPILED_OS_VERSION QUERY OS_VERSION)
cmake_host_system_information(RESULT COMPILED_OS_PLATFORM QUERY OS_PLATFORM)
cmake_host_system_information(RESULT draco_is64bit QUERY IS_64BIT)
cmake_host_system_information(RESULT HAVE_HARDWARE_FPU QUERY HAS_FPU)
cmake_host_system_information(RESULT HAVE_HARDWARE_MMX QUERY HAS_MMX)
cmake_host_system_information(RESULT HAVE_HARDWARE_MMX_PLUS QUERY HAS_MMX_PLUS)
cmake_host_system_information(RESULT HAVE_HARDWARE_SSE QUERY HAS_SSE)
cmake_host_system_information(RESULT HAVE_HARDWARE_SSE2 QUERY HAS_SSE2)
cmake_host_system_information(RESULT HAVE_HARDWARE_SSE_FP QUERY HAS_SSE_FP)
cmake_host_system_information(RESULT HAVE_HARDWARE_SSE_MMX QUERY HAS_SSE_MMX)

# SYS_TYPE
#
# * RHEL7: blueos_3_ppc64le_ib_p9
# * REHL8: toss_4_ppc64le_ib
if(${SITENAME_FAMILY} MATCHES "ATS-2" AND "$ENV{SYS_TYPE}" MATCHES "blueos_3")
  set(draco_isATS2_rhel7 1)
endif()

# Create dsxx/config.h
configure_file(config.h.in ${PROJECT_BINARY_DIR}/dsxx/config.h)

# Generate the dll_declspec.h file that defines DLL_PUBLIC_<pkg> CPP macros.
include(generate_dll_declspec)
set(components
    # cmake-format: sortable
    c4
    cdi
    cdi_analytic
    cdi_CPEless_test
    cdi_CPEloss
    cdi_eospac
    cdi_ipcress
    cdi_ipcress_test
    cdi_ndi
    cdi_ndi_test
    cdi_test
    device
    diagnostics
    dsxx
    FC_Derived_Type
    fit
    linear
    mesh
    mesh_element
    meshReaders
    min
    norms
    ode
    parser
    quadrature
    quadrature_test
    rng
    roots
    RTT_Format_Reader
    special_functions
    units
    viz)
generate_dll_declspec("dsxx" "${components}")

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
file(GLOB template_implementations *.t.hh *.i.hh)
file(GLOB sources *.cc)
file(GLOB explicit_instantiations *_pt.cc)
file(GLOB dsxx_headers *.hh)
list(APPEND dsxx_headers ${PROJECT_BINARY_DIR}/dsxx/config.h
     ${PROJECT_BINARY_DIR}/dsxx/device_config.h ${PROJECT_BINARY_DIR}/dsxx/dll_declspec.h)
set(dsxx_terminal_headers ${PROJECT_SOURCE_DIR}/terminal/terminal.h
                          ${PROJECT_SOURCE_DIR}/terminal/terminal_base.h)
list(REMOVE_ITEM dsxx_headers ${template_implementations})

if(${CMAKE_GENERATOR} MATCHES "Xcode")
  # Xcode does not support per file, per configuration compile definitions, so we add this compile
  # definition for all files
  add_compile_definitions(BUILD_TYPE="$<CONFIG>")
else()
  set_source_files_properties(Release.cc PROPERTIES COMPILE_DEFINITIONS BUILD_TYPE="$<CONFIG>")
endif()

# ------------------------------------------------------------------------------------------------ #
# Build package library
# ------------------------------------------------------------------------------------------------ #
list(APPEND include_dirs PUBLIC $<BUILD_INTERFACE:${Draco_SOURCE_DIR}/src>
     $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> $<INSTALL_INTERFACE:${DBSCFG_IMPORT_PREFIX}>)

set(target_deps)

if(TARGET coverage_config)
  list(APPEND target_deps coverage_config)
endif()

add_component_library(
  TARGET Lib_dsxx
  LIBRARY_NAME dsxx
  SOURCES "${sources}"
  HEADERS "${dsxx_headers};${dsxx_terminal_headers}"
  INCLUDE_DIRS "${include_dirs}"
  TARGET_DEPS "${target_deps}")

# ------------------------------------------------------------------------------------------------ #
# Installation instructions
# ------------------------------------------------------------------------------------------------ #
install(TARGETS coverage_config EXPORT draco-targets)
install(FILES ${dsxx_headers} ${template_implementations} DESTINATION ${DBSCFGDIR}include/dsxx)
install(FILES ${dsxx_terminal_headers} ${PROJECT_SOURCE_DIR}/terminal/LICENSE
              ${PROJECT_SOURCE_DIR}/terminal/README.md DESTINATION ${DBSCFGDIR}include/terminal)

# ------------------------------------------------------------------------------------------------ #
# Unit tests
# ------------------------------------------------------------------------------------------------ #
if(BUILD_TESTING)
  add_subdirectory(test)
endif()
if(WIN32 AND (NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
  add_subdirectory(bin)
endif()

# ------------------------------------------------------------------------------------------------ #
# Autodoc
# ------------------------------------------------------------------------------------------------ #
process_autodoc_pages()

# ------------------------------------------------------------------------------------------------ #
# End dsxx/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
