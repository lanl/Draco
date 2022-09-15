#--------------------------------------------*-cmake-*---------------------------------------------#
# file   FortranChecks/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for FortranChecks.
# note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(FortranChecks CXX)

# When using CXX/Fortran cross language linking and function calling the CMake module
# FortranCInterface.cmake can be used to deal with name mangling issues.  However, Draco has chosen
# to use ISO_C_BINDING instead and this makes manual name demangling unnecessary.

if(NOT BUILD_TESTING OR NOT HAVE_Fortran)
  return()
endif()

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
file(GLOB Csources *.cc)

# ------------------------------------------------------------------------------------------------ #
# Build package libraries
# ------------------------------------------------------------------------------------------------ #
add_component_library(
  TARGET Lib_FC_Derived_Type
  TARGET_DEPS "Lib_dsxx"
  LIBRARY_NAME FC_Derived_Type
  SOURCES "${Csources}"
  NOEXPORT)
set_target_properties(
  Lib_FC_Derived_Type
  PROPERTIES
  FOLDER FortranChecks)

# ------------------------------------------------------------------------------------------------ #
# Build a library from Fortran sources
# ------------------------------------------------------------------------------------------------ #

# For Unix, this will simply run add_subdirectory(f90sub). For Visual Studio or XCode, an
# ExternalProject will be created that uses an alternate gfortran via Makefiles to generate the
# library external to the main project.
include(CMakeAddFortranSubdirectory)
init_build_system_state() # returns ${build_system_state} used below.

# CMake will put some run scripts in $build_dir/f90sub Another directory for the build is generated
# at $build_dir/FortranChecks_f90sub_build-prefix
cmake_add_fortran_subdirectory(
  f90sub # directory name
  PROJECT FortranChecks_f90sub # project name in top level CMakeLists.txt Creates target named
                               # 'FortranChecks_f90sub_build'
  ARCHIVE_DIR f90sub # .lib location relative to root binary tree
  RUNTIME_DIR f90sub # .dll location relative to root binary tree
  LIBRARIES rtt_FC_f90sub # create a library with this name. --> librtt_FC_f90.dll.
  TARGET_NAMES Lib_FC_f90sub
  DEPENDS "Lib_FC_Derived_Type;Lib_c4;Lib_dsxx"
  NO_EXTERNAL_INSTALL
  CMAKE_COMMAND_LINE "${build_system_state}"
  # VERBOSE
)

add_subdirectory(test)

# ------------------------------------------------------------------------------------------------ #
# Autodoc
# ------------------------------------------------------------------------------------------------ #
process_autodoc_pages()

# ------------------------------------------------------------------------------------------------ #
# end FortranChecks/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
