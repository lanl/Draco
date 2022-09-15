#--------------------------------------------*-cmake-*---------------------------------------------#
# file   FortranChecks/f90sub/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2013 June 24
# brief  Generate build project files for Fortran.
# note   Copyright (C) 2014-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(FortranChecks_f90sub C Fortran)

# ------------------------------------------------------------------------------------------------ #
# If this directory is built via CAFS (ie. MSVC or XCode), then we need to do some extra setup. When
# build via CAFS, target Lib_dsxx won't be defined yet.
# ------------------------------------------------------------------------------------------------ #
if(NOT TARGET Lib_dsxx)
  include(CMakeAddFortranSubdirectory)
  cafs_fortran_dir_boilerplate_setup()

  # CAFS setup unique to this directory.
  cafs_create_imported_targets(Lib_FC_Derived_Type "rtt_FC_Derived_Type"
                               "${draco_BINARY_DIR}/src/FortranChecks" CXX)
endif()

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
file(GLOB f90sources *.f90 *.F90)

# ------------------------------------------------------------------------------------------------ #
# Build library for test directory
# ------------------------------------------------------------------------------------------------ #

# Xcode: Since cafs_create_imported_targets does not set all of the dependencies between the
# imported libraries, we need to list them explicitly.
set(target_deps Lib_FC_Derived_Type;Lib_dsxx)
add_component_library(
  TARGET Lib_FC_f90sub
  LIBRARY_NAME FC_f90sub
  SOURCES "${f90sources}"
  LINK_LANGUAGE "Fortran"
  TARGET_DEPS "${target_deps}"
  NOEXPORT)

# ------------------------------------------------------------------------------------------------ #
# End f90sub/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
