#--------------------------------------------*-cmake-*---------------------------------------------#
# file   quadrature/CMakeLists.txt
# brief  Instructions for building quadrature Makefiles.
# note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(quadrature CXX)

# Local variables ---------------------------------------------------------------------------------#
string(CONCAT header_prefix_source "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>"
              "$<INSTALL_INTERFACE:${DBSCFG_IMPORT_PREFIX}/quadrature>")

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
file(
  GLOB sources
  RELATIVE "${PROJECT_SOURCE_DIR}"
  *.cc)
file(
  GLOB headers
  RELATIVE "${PROJECT_SOURCE_DIR}"
  *.hh)

# prefix header files with appropriate paths (build vs install).
foreach(file ${headers})
  list(APPEND headers_public ${header_prefix_source}/${file})
endforeach()

# ------------------------------------------------------------------------------------------------ #
# Build package library
# ------------------------------------------------------------------------------------------------ #
add_component_library(
  TARGET Lib_quadrature
  TARGET_DEPS "Lib_special_functions;Lib_parser;Lib_mesh_element"
  LIBRARY_NAME ${PROJECT_NAME}
  SOURCES "${sources}"
  HEADERS "${headers_public}")
# NOEXPORT) # handle installation instructions below.

# target_sources( Lib_quadrature PRIVATE "${headers_public}" "${sources}" ) set_target_properties(
# Lib_quadrature PROPERTIES PUBLIC_HEADER "${headers_public}") if(TARGET Objlib_quadrature)
# target_sources( Objlib_quadrature PRIVATE "${headers_public}" "${sources}" )
# set_target_properties( Objlib_quadrature PROPERTIES PUBLIC_HEADER "${headers_public}") endif()

# ------------------------------------------------------------------------------------------------ #
# Fortran side of the API
# ------------------------------------------------------------------------------------------------ #
if(HAVE_Fortran)

  # For Unix, this will simply run add_subdirectory(fc4). For Visual Studio or XCode, an
  # ExternalProject will be created that uses an alternate gfortran via Makefiles to generate the
  # library external to the main project.
  include(CMakeAddFortranSubdirectory)
  init_build_system_state() # returns ${build_system_state} used below.

  cmake_add_fortran_subdirectory(
    fquadrature # directory name
    PROJECT quadrature_fquadrature # project name in top level CMakeLists.txt. Creates target named
            # 'quadrature_fquadrature_build'
    ARCHIVE_DIR fquadrature # .lib location relative to root binary tree
    RUNTIME_DIR fquadrature # .dll location relative to root binary tree
    LIBRARIES rtt_fquadrature # new library will have this name. --> librtt_fquadrature.dll.
    TARGET_NAMES Lib_quadrature_fquadrature
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
# install( TARGETS Lib_quadrature EXPORT draco-targets LIBRARY DESTINATION "${DBSCFGDIR_LIBRARY}"
# ARCHIVE DESTINATION "${DBSCFGDIR}lib" RUNTIME DESTINATION "${DBSCFGDIR}bin" PUBLIC_HEADER
# DESTINATION "${DBSCFGDIR}include/quadrature" ) if(TARGET Objlib_quadrature) install( TARGETS
# Objlib_quadrature EXPORT draco-targets LIBRARY DESTINATION "${DBSCFGDIR_LIBRARY}" ARCHIVE
# DESTINATION "${DBSCFGDIR}lib" RUNTIME DESTINATION "${DBSCFGDIR}bin" PUBLIC_HEADER DESTINATION
# "${DBSCFGDIR}include/quadrature" ) endif()
install(FILES ${headers_public} DESTINATION ${DBSCFGDIR}include/quadrature)

# ------------------------------------------------------------------------------------------------ #
# Unit tests
# ------------------------------------------------------------------------------------------------ #
if(BUILD_TESTING)
  add_subdirectory(test)

  if(HAVE_Fortran)
    cmake_add_fortran_subdirectory(
      ftest # directory name
      PROJECT quadrature_ftest # project name in top level CMakeLists.txt. Creates target named
                               # 'quadrature_ftest_build'
      ARCHIVE_DIR ftest # .lib location relative to root binary tree
      RUNTIME_DIR ftest # .dll location relative to root binary tree
      LIBRARIES rtt_quadrature_ftest # new library will have this name.
                # --> librtt_quadrature_ftest.dll.
      TARGET_NAMES Lib_quadrature_ftest
      DEPENDS "Lib_quadrature_test;Lib_quadrature;Lib_quadrature_fquadrature" NO_EXTERNAL_INSTALL
      CMAKE_COMMAND_LINE "${build_system_state}"
      # VERBOSE
    )
  endif()

endif()

# ------------------------------------------------------------------------------------------------ #
# Autodoc
# ------------------------------------------------------------------------------------------------ #
process_autodoc_pages()

# ------------------------------------------------------------------------------------------------ #
# End quadrature/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
