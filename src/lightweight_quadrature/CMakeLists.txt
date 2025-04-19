#--------------------------------------------*-cmake-*---------------------------------------------#
# file   lightweight_quadrature/CMakeLists.txt
# author Jim Warsa <warsa@lanl.gov>
# date   2025 Apr 1
# brief  Generate build project files for lightweight_quadrature
# note   Copyright (C) 2025 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.21.0)
project(lightweight_quadrature CXX)

if(TARGET GSL::gsl)

  # ---------------------------------------------------------------------------------------------- #
  # Source files
  # ---------------------------------------------------------------------------------------------- #
  file(GLOB template_implementations *.t.hh *.i.hh)
  file(GLOB sources *.cc)
  file(GLOB explicit_instantiations *_pt.cc)
  file(GLOB headers *.hh)
  list(REMOVE_ITEM headers ${template_implementations})

  # ---------------------------------------------------------------------------------------------- #
  # Build package library
  # ---------------------------------------------------------------------------------------------- #
  add_component_library(
    TARGET Lib_lightweight_quadrature
    TARGET_DEPS "Lib_dsxx;Lib_c4;Lib_units;GSL::gsl"
    LIBRARY_NAME ${PROJECT_NAME}
    HEADERS "${headers}"
    SOURCES "${sources}")

  # ---------------------------------------------------------------------------------------------- #
  # Installation instructions
  # ---------------------------------------------------------------------------------------------- #
  install(FILES ${headers} ${template_implementations}
          DESTINATION ${DBSCFGDIR}include/lightweight_quadrature)

  # ---------------------------------------------------------------------------------------------- #
  # Unit tests
  # ---------------------------------------------------------------------------------------------- #
  if(BUILD_TESTING)
    add_subdirectory(test)
  endif()

  # ---------------------------------------------------------------------------------------------- #
  # Autodoc
  # ---------------------------------------------------------------------------------------------- #
  process_autodoc_pages()

endif(TARGET GSL::gsl)

# ------------------------------------------------------------------------------------------------ #
# End lightweight_quadrature/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
