#--------------------------------------------*-cmake-*---------------------------------------------#
# file   quadrature/ftest/CMakeLists.txt
# author Allan Wollaber
# date   Mon May 23 16:15:44 MDT 2016
# brief  Generate a Fortran testing library for quadrature
# note   Copyright (C) 2016-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(quadrature_ftest C Fortran)

# ------------------------------------------------------------------------------------------------ #
# If this directory is built via CAFS (ie. MSVC or XCode), then we need to do some extra setup. When
# build via CAFS, target Lib_dsxx won't be defined yet.
# ------------------------------------------------------------------------------------------------ #
if(NOT TARGET Lib_dsxx)
  include(CMakeAddFortranSubdirectory)
  cafs_fortran_dir_boilerplate_setup()
  cafs_fortran_dir_mpi_setup() # returns CAFS_MPI_DEPS

  # CAFS setup unique to this directory. 1. create import targets that we will link against.
  cafs_create_imported_targets(Lib_quadrature "rtt_quadrature" "${draco_BINARY_DIR}/src/quadrature"
                               CXX)
  cafs_create_imported_targets(Lib_quadrature_test "rtt_quadrature_test"
                               "${draco_BINARY_DIR}/src/quadrature/test" CXX)
  cafs_create_imported_targets(Lib_quadrature_fquadrature "librtt_fquadrature"
                               "${draco_BINARY_DIR}/src/quadrature/fquadrature" CXX)

  # 1. Vendor discovery is not needed for this package setupVendorLibraries()
else()
  set(draco_BINARY_DIR ${Draco_BINARY_DIR})
endif()

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
set(f90sources ${PROJECT_SOURCE_DIR}/tstquadrature_interfaces.f90)

# ------------------------------------------------------------------------------------------------ #
# Build library for test directory
# ------------------------------------------------------------------------------------------------ #
set(target_deps Lib_quadrature_test Lib_quadrature_fquadrature Lib_quadrature Lib_dsxx
                ${CAFS_MPI_DEPS})
set(include_dirs PUBLIC "$<BUILD_INTERFACE:${quadrature_ftest_BINARY_DIR}>"
                 "$<BUILD_INTERFACE:${draco_BINARY_DIR}/src/quadrature/fquadrature>") # .mod files
add_component_library(
  TARGET Lib_quadrature_ftest
  LIBRARY_NAME quadrature_ftest
  SOURCES "${f90sources}"
  LINK_LANGUAGE Fortran
  TARGET_DEPS "${target_deps}"
  INCLUDE_DIRS "${include_dirs}"
  NOEXPORT)

# ------------------------------------------------------------------------------------------------ #
# end quadrature/ftest/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
