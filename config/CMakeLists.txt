#--------------------------------------------*-cmake-*---------------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 Nov 2
# brief  Install the config directory in the target directory
# note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(config)

include(CMakePackageConfigHelpers)

file(GLOB CMake_src *.cmake)
file(GLOB config_helper_sh *.sh)
file(GLOB Python_src *.py)
set(CMake_in cmake_uninstall.cmake.in doxygen_config.in)
set(CMake_AFSD CMakeAddFortranSubdirectory/config_cafs_proj.cmake.in)

# include(debug_macros) echo_all_cmake_variable_values()

# ------------------------------------------------------------------------------------------------ #
# Prepare data that will be saved in draco-config.cmake
# ------------------------------------------------------------------------------------------------ #

string(TOUPPER "${CMAKE_BUILD_TYPE}" upper_build_type)
set(Draco_CXX_COMPILER_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${upper_build_type}}")
set(Draco_C_COMPILER_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${upper_build_type}}")
set(Draco_Fortran_COMPILER_FLAGS
    "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${upper_build_type}}")

# What MPI import targets should be found by clients (captured in draco-config.cmake).
if("${DRACO_C4}" STREQUAL "MPI")
  # c4 requires MPI::MPI_CXX if DRACO_C4==MPI.
  set(Draco_MPI_LANG_LIST "CXX")
  # If Draco found other MPI imported targets, then ensure they are defined for client software
  # packages (e.g.: MPI::MPI_Fortran)
  foreach(lang C Fortran)
    if(TARGET MPI::MPI_${lang})
      list(APPEND Draco_MPI_LANG_LIST ${lang})
    endif()
  endforeach()
endif()

set(Draco_TPL_LIST GSL Random123) # except MPI (which is treated separately).
if(TARGET EOSPAC::eospac)
  list(APPEND Draco_TPL_LIST EOSPAC)
endif()

# Process the config file
configure_file(draco-config-install.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/draco-config.cmake" @ONLY)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/draco-config-version.cmake"
  VERSION ${Draco_VERSION}
  COMPATIBILITY AnyNewerVersion)

# Install scripts and macros to make them available by other projects.
set(file_list
    ${CMake_src}
    ${CMake_in}
    ${Python_src}
    ${CMAKE_CURRENT_BINARY_DIR}/draco-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/draco-config-version.cmake
    ${Draco_SOURCE_DIR}/autodoc/html/doxygen.css
    ${Draco_SOURCE_DIR}/.gitlab/ci/valgrind_suppress.txt)

install(FILES ${file_list} DESTINATION cmake)
install(FILES ${CMake_AFSD} DESTINATION cmake/CMakeAddFortranSubdirectory)
install(PROGRAMS ${config_helper_sh} DESTINATION cmake)

# ---------------------------------------------------------------------------##
# End of config/CMakeLists.txt
# ---------------------------------------------------------------------------##
