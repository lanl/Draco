#--------------------------------------------*-cmake-*---------------------------------------------#
# file   c4/fc4/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2020 April 20
# brief  Generate build project files for Fortran.
# note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(fc4 C Fortran)

# ------------------------------------------------------------------------------------------------ #
# If this directory is built via CAFS (ie. MSVC or XCode), then we need to do some extra setup. When
# build via CAFS, target Lib_dsxx won't be defined yet.
# ------------------------------------------------------------------------------------------------ #
set(target_deps Lib_c4 Lib_dsxx)
if(NOT TARGET Lib_dsxx)
  include(CMakeAddFortranSubdirectory)
  cafs_fortran_dir_boilerplate_setup()
  cafs_fortran_dir_mpi_setup() # returns CAFS_MPI_DEPS
elseif("${DRACO_C4}" STREQUAL "MPI")
  list(APPEND target_deps MPI::MPI_Fortran)
  add_definitions(-DC4_MPI)
endif()

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
file(GLOB f90sources *.f90 *.F90)

# gfortran-10 needs to be told to allow 'invalid' initialization value.
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL
                                              "10.0")
  set_source_files_properties(Draco_MPI.F90 PROPERTIES COMPILE_FLAGS "-fallow-invalid-boz")
endif()

# ------------------------------------------------------------------------------------------------ #
# Build library
# ------------------------------------------------------------------------------------------------ #
add_component_library(
  TARGET Lib_c4_fc4
  LIBRARY_NAME fc4
  SOURCES "${f90sources}"
  LINK_LANGUAGE Fortran
  TARGET_DEPS "${target_deps}"
  INCLUDE_DIRS "PUBLIC;$<INSTALL_INTERFACE:${DBSCFG_IMPORT_PREFIX}/fc4>"
  NOEXPORT)
if(WIN32)
  # The library must be installed in the master project's RUNTIME_OUTPUT_DIR
  set_target_properties(
    Lib_c4_fc4
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${draco_BINARY_DIR}/${CMAKE_BUILD_TYPE})
endif()

# ------------------------------------------------------------------------------------------------ #
# Install
# ------------------------------------------------------------------------------------------------ #
if(TARGET Objlib_c4_fc4)
  install(
    TARGETS Objlib_c4_fc4
    EXPORT draco-targets
    LIBRARY DESTINATION "${DBSCFGDIR_LIBRARY}"
    ARCHIVE DESTINATION "${DBSCFGDIR}lib"
    RUNTIME DESTINATION "${DBSCFGDIR}bin"
    PUBLIC_HEADER DESTINATION "${DBSCFGDIR}include/fc4")
endif()
install(
  TARGETS Lib_c4_fc4
  EXPORT draco-targets
  LIBRARY DESTINATION "${DBSCFGDIR_LIBRARY}"
  ARCHIVE DESTINATION "${DBSCFGDIR}lib"
  RUNTIME DESTINATION "${DBSCFGDIR}bin"
  PUBLIC_HEADER DESTINATION "${DBSCFGDIR}include/fc4")
install(
  DIRECTORY ${fc4_BINARY_DIR}
  DESTINATION "${DBSCFGDIR}include"
  FILES_MATCHING
  PATTERN "*.mod")

# ------------------------------------------------------------------------------------------------ #
# End c4/fc4/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
