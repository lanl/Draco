#--------------------------------------------*-cmake-*---------------------------------------------#
#
# file   predict/CMakeLists.txt
# author Mathew cleveland
# date   Jan. 11 2021
# brief  Provided an interface for generating predictions from pre-trained neural networks.
# note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved.
#
#--------------------------------------------------------------------------------------------------#

project(predict CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #

file(GLOB sources *.cc)
file(GLOB headers *.hh)

# Make the header files available in the IDE.
if(MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode)
  list(APPEND sources ${headers})
endif()

# ------------------------------------------------------------------------------------------------ #
# Build package library
# ------------------------------------------------------------------------------------------------ #
set(predict_target_deps Lib_dsxx Lib_c4)
if(TORCH_FOUND)
  list(APPEND predict_target_deps torch)
endif()

add_component_library(
  TARGET Lib_predict
  TARGET_DEPS ${predict_target_deps}
  INCLUDE_DIRS ${include_dirs}
  LIBRARY_NAME ${PROJECT_NAME}
  HEADERS "${headers}"
  SOURCES "${sources}")

# ------------------------------------------------------------------------------------------------ #
# Installation instructions
# ------------------------------------------------------------------------------------------------ #

install(FILES ${headers} DESTINATION ${DBSCFGDIR}include/predict)

# ------------------------------------------------------------------------------------------------ #
# Unit tests
# ------------------------------------------------------------------------------------------------ #

if(BUILD_TESTING)
  add_subdirectory(test)
endif()

# ------------------------------------------------------------------------------------------------ #
# Autodoc
# ------------------------------------------------------------------------------------------------ #

process_autodoc_pages()
# ------------------------------------------------------------------------------------------------ #
# End predict/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
