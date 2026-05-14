#--------------------------------------------*-cmake-*---------------------------------------------#
# file   special_functions/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for special_functions.
# note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(special_functions CXX)

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
    TARGET Lib_special_functions
    TARGET_DEPS "Lib_units;Lib_roots;GSL::gsl"
    LIBRARY_NAME ${PROJECT_NAME}
    HEADERS "${headers}"
    SOURCES "${sources}")

  # ---------------------------------------------------------------------------------------------- #
  # Installation instructions
  # ---------------------------------------------------------------------------------------------- #
  install(FILES ${headers} ${template_implementations}
          DESTINATION ${DBSCFGDIR}include/special_functions)

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
# End special_functions/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
