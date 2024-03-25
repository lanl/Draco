#--------------------------------------------*-cmake-*---------------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 Nov 2
# brief  Install the config directory in the target directory
# note   Copyright (C) 2010-2024 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(config)

include(CMakePackageConfigHelpers)

file(GLOB CMake_src *.cmake)
file(GLOB config_helper_sh *.sh)
file(GLOB Python_src *.py)
file(GLOB Json_src *.json)
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
    ${Json_src}
    ${CMAKE_CURRENT_BINARY_DIR}/draco-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/draco-config-version.cmake
    ${Draco_SOURCE_DIR}/autodoc/html/doxygen.css
    ${Draco_SOURCE_DIR}/.gitlab/ci/valgrind_suppress.txt)

install(FILES ${file_list} DESTINATION cmake)
install(FILES ${CMake_AFSD} DESTINATION cmake/CMakeAddFortranSubdirectory)
install(PROGRAMS ${config_helper_sh} DESTINATION cmake)

# ---------------------------------------------------------------------------##
# Rules to submit build to CDash
# ---------------------------------------------------------------------------##

add_custom_command(
  OUTPUT ${PROJECT_BINARY_DIR}/always_out_of_date.log
  COMMAND "${CMAKE_PROGRAM_PATH}" -E make_directory "${Draco_BINARY_DIR}/cdashrep"
  COMMAND
    ${Draco_SOURCE_DIR}/.gitlab/ci/gitlab-nightly-regress.sh -e -g "${CMAKE_GENERATOR}" -i
    "${Draco_BINARY_DIR}/cdashrep" -p "Draco" -t "${CMAKE_BUILD_TYPE}" -s "${Draco_SOURCE_DIR}" -n
    "${MPIEXEC_MAX_NUMPROCS}"
    # -x 24 -v "-VV" -m "Configure,Build,Test,Submit" (defaulted)
  COMMENT "starting...")

add_custom_target(
  cdashrep
  DEPENDS ${PROJECT_BINARY_DIR}/always_out_of_date.log
  COMMENT
    "Report build to CDash: https://rtt.lanl.gov/cdash3/index.php?project=Draco

This target runs the script:

  ${Draco_SOURCE_DIR}/.gitlab/ci/gitlab-nightly-regress.sh -e -g \"${CMAKE_GENERATOR}\"
  -i ${Draco_BINARY_DIR}/cdashrep -p Draco -t ${CMAKE_BUILD_TYPE} -s ${Draco_SOURCE_DIR}
  -n ${MPIEXEC_MAX_NUMPROCS}

Required environment variables:
  cdash_upload_auth_token")

# ---------------------------------------------------------------------------##
# End of config/CMakeLists.txt
# ---------------------------------------------------------------------------##
