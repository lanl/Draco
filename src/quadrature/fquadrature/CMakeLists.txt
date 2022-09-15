#--------------------------------------------*-cmake-*---------------------------------------------#
# file   quadrature/fquadrature/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2020 October 13
# brief  Generate build project files for Fortran.
# note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(fquadrature C Fortran)

# ------------------------------------------------------------------------------------------------ #
# If this directory is built via CAFS (ie. MSVC or XCode), then we need to do some extra setup. When
# build via CAFS, target Lib_dsxx won't be defined yet.
# ------------------------------------------------------------------------------------------------ #
if(NOT TARGET Lib_dsxx)
  include(CMakeAddFortranSubdirectory)
  cafs_fortran_dir_boilerplate_setup()
  cafs_fortran_dir_mpi_setup() # returns CAFS_MPI_DEPS
endif()

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
file(GLOB f90sources *.f90 *.F90)

# ------------------------------------------------------------------------------------------------ #
# Build library
# ------------------------------------------------------------------------------------------------ #
add_component_library(
  TARGET Lib_quadrature_fquadrature
  LIBRARY_NAME fquadrature
  SOURCES "${f90sources}"
  LINK_LANGUAGE Fortran
  INCLUDE_DIRS "PUBLIC;$<INSTALL_INTERFACE:${DBSCFG_IMPORT_PREFIX}/fquadrature>"
  NOEXPORT)
if(WIN32)
  # The library must be installed in the master project's RUNTIME_OUTPUT_DIR
  set_target_properties(
    Lib_quadrature_fquadrature
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${draco_BINARY_DIR}/${CMAKE_BUILD_TYPE})
endif()

# ------------------------------------------------------------------------------------------------ #
# Install
# ------------------------------------------------------------------------------------------------ #
if(TARGET Objlib_quadrature_fquadrature)
  install(
    TARGETS Objlib_quadrature_fquadrature
    EXPORT draco-targets
    LIBRARY DESTINATION "${DBSCFGDIR_LIBRARY}"
    ARCHIVE DESTINATION "${DBSCFGDIR}lib"
    RUNTIME DESTINATION "${DBSCFGDIR}bin"
    PUBLIC_HEADER DESTINATION "${DBSCFGDIR}include/fquadrature")
endif()
install(
  TARGETS Lib_quadrature_fquadrature
  EXPORT draco-targets
  LIBRARY DESTINATION "${DBSCFGDIR_LIBRARY}"
  ARCHIVE DESTINATION "${DBSCFGDIR}lib"
  RUNTIME DESTINATION "${DBSCFGDIR}bin"
  PUBLIC_HEADER DESTINATION "${DBSCFGDIR}include/fquadrature")
install(
  DIRECTORY ${fquadrature_BINARY_DIR}
  DESTINATION "${DBSCFGDIR}include"
  FILES_MATCHING
  PATTERN "*.mod")

# ------------------------------------------------------------------------------------------------ #
# End quadrature/fquadrature/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
