#--------------------------------------------*-cmake-*---------------------------------------------#
# file   c4/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for c4
# note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(c4 CXX)

# ------------------------------------------------------------------------------------------------ #
# Quo related setup that will be saved to c4/config.h
# ------------------------------------------------------------------------------------------------ #
if(TARGET LIBQUO::libquo AND DRACO_C4 STREQUAL MPI)
  set(HAVE_LIBQUO 1)
else()
  set(HAVE_LIBQUO 0)
endif()
set(HAVE_LIBQUO
    ${HAVE_LIBQUO}
    CACHE BOOL "Does Draco require the libquo library?" FORCE)

# ------------------------------------------------------------------------------------------------ #
# Generate config.h (only occurs when cmake is run) many c4 and MPI values are set in
# config/setupMPI.cmake.
# ------------------------------------------------------------------------------------------------ #
configure_file(config.h.in ${PROJECT_BINARY_DIR}/c4/config.h)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
file(GLOB sources *.cc)
file(GLOB headers *.hh *.h)
list(APPEND headers ${PROJECT_BINARY_DIR}/c4/config.h)

# ------------------------------------------------------------------------------------------------ #
# Build package library
# ------------------------------------------------------------------------------------------------ #
if(PAPI_LIBRARY)
  set(c4_extra_VENDORS VENDOR_LIST;PAPI;VENDOR_LIBS;${PAPI_LIBRARY})
endif()

set(target_deps Lib_dsxx)
if("${DRACO_C4}" STREQUAL "MPI")
  if(HAVE_LIBQUO)
    list(APPEND target_deps LIBQUO::libquo)
  endif()
  list(APPEND target_deps MPI::MPI_CXX MPI::MPI_C)
endif()

# Manually link to the OpenMP runtime for cl-clang under MSVC.  This is work-around for a cmake bug
# and hopefully can be removed once the bug is fixed.
if(WIN32
   AND (CMAKE_C_COMPILER_ID STREQUAL "Clang")
   AND OpenMP_FOUND)
  list(APPEND target_deps OpenMP::OpenMP_C)
endif()

add_component_library(
  ${c4_extra_VENDORS}
  TARGET Lib_c4
  TARGET_DEPS ${target_deps}
  INCLUDE_DIRS "PUBLIC;$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>"
  LIBRARY_NAME c4
  SOURCES "${sources}"
  HEADERS "${headers}")

# xthi, ythi binaries
add_subdirectory(bin)

# ------------------------------------------------------------------------------------------------ #
# Fortran-based MPI wrapper
# ------------------------------------------------------------------------------------------------ #
if(HAVE_Fortran)

  set(deplist Lib_c4 Lib_dsxx)

  # For Unix, this will simply run add_subdirectory(fc4). For Visual Studio or XCode, an
  # ExternalProject will be created that uses an alternate gfortran via Makefiles to generate the
  # library external to the main project.
  include(CMakeAddFortranSubdirectory)
  init_build_system_state() # returns ${build_system_state} used below.

  cmake_add_fortran_subdirectory(
    fc4 #< directory name
    PROJECT c4_fc4 #< project name in top level CMakeLists.txt. Creates tgt 'c4_fc4_build'
    ARCHIVE_DIR fc4 #< .lib location relative to root binary tree
    RUNTIME_DIR fc4 #< .dll location relative to root binary tree
    LIBRARIES rtt_fc4 #< new library will have this name. --> librtt_fc4.dll.
    TARGET_NAMES Lib_c4_fc4
    DEPENDS "${deplist}"
    CMAKE_COMMAND_LINE "${build_system_state}"
    # VERBOSE
  )
  string(APPEND Draco_EXPORT_TARGET_PROPERTIES "${CAFS_EXPORT_DEFINE_IMPORT_PREFIX}"
         "${CAFS_EXPORT_TARGET_PROPERTIES}")
  set(Draco_EXPORT_TARGET_PROPERTIES
      "${Draco_EXPORT_TARGET_PROPERTIES}"
      PARENT_SCOPE)

endif()

# ------------------------------------------------------------------------------------------------ #
# Installation instructions
# ------------------------------------------------------------------------------------------------ #
install(FILES ${headers} DESTINATION ${DBSCFGDIR}include/c4)

# ------------------------------------------------------------------------------------------------ #
# Unit tests
# ------------------------------------------------------------------------------------------------ #
if(BUILD_TESTING)
  add_subdirectory(test)

  if(HAVE_Fortran)
    cmake_add_fortran_subdirectory(
      ftest #< directory name
      NO_EXTERNAL_INSTALL
      # VERBOSE
      PROJECT c4_ftest #< project name in top level CMakeLists.txt. Creates target named
                       #< 'c4_ftest_build'
      ARCHIVE_DIR ftest #< .lib location relative to root binary tree
      RUNTIME_DIR ftest #< .dll location relative to root binary tree
      LIBRARIES rtt_c4_ftest #< new library will have this name. --> librtt_c4_ftest.dll.
      TARGET_NAMES Lib_c4_ftest
      DEPENDS "Lib_c4_fc4;${deplist}"
      CMAKE_COMMAND_LINE "${build_system_state}")
  endif()

endif()

# ------------------------------------------------------------------------------------------------ #
# Autodoc
# ------------------------------------------------------------------------------------------------ #
process_autodoc_pages()

# ------------------------------------------------------------------------------------------------ #
# End c4/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
