#--------------------------------------------*-cmake-*---------------------------------------------#
#
# file   predict/CMakeLists.txt
# author Mathew cleveland
# date   Jan. 11 2021
# brief  Provided an interface for generating predictions from pre-trained neural networks.
# note   Copyright (C) 2018-2020 Triad National Security, LLC.
#        All rights reserved.
#--------------------------------------------------------------------------------------------------#
project( predict CXX )

#--------------------------------------------------------------------------------------------------#
# Source files
#--------------------------------------------------------------------------------------------------#

file( GLOB sources *.cc )
file( GLOB headers *.hh )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers})
endif()

#--------------------------------------------------------------------------------------------------#
# Build package library
#--------------------------------------------------------------------------------------------------#

if(NOT Torch_FOUND)
    message(STATUS "libtorch - not found")
    add_component_library(
        TARGET       Lib_predict
        TARGET_DEPS  Lib_dsxx Lib_c4
        INCLUDE_DIRS ${include_dirs}
        LIBRARY_NAME ${PROJECT_NAME}
        HEADERS      "${headers}"
        SOURCES      "${sources}"
    )
else()
    add_definitions(-DLIBTORCH)
    add_component_library(
        TARGET       Lib_predict
        TARGET_DEPS  torch Lib_dsxx Lib_c4
        INCLUDE_DIRS ${include_dirs}
        LIBRARY_NAME ${PROJECT_NAME}
        HEADERS      "${headers}"
        SOURCES      "${sources}"
    )
endif()


#--------------------------------------------------------------------------------------------------#
# Installation instructions
#--------------------------------------------------------------------------------------------------#

install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/predict )

#--------------------------------------------------------------------------------------------------#
# Unit tests
#--------------------------------------------------------------------------------------------------#

if( BUILD_TESTING )
   add_subdirectory( test )
endif()

#--------------------------------------------------------------------------------------------------#
# Autodoc
#--------------------------------------------------------------------------------------------------#

process_autodoc_pages()
#--------------------------------------------------------------------------------------------------#
# End predict/CMakeLists.txt
#--------------------------------------------------------------------------------------------------#

