#--------------------------------------------*-cmake-*---------------------------------------------#
# file   device/CMakeLists.txt
# author Gabriel Rockefeller
# date   2011 June 13
# brief  Instructions for building device Makefile.
# note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(device CXX ${CUDA_DBS_STRING})

# ------------------------------------------------------------------------------------------------ #
# Project configuration file
# ------------------------------------------------------------------------------------------------ #
set(TEST_KERNEL_BINDIR
    ${PROJECT_BINARY_DIR}/test
    CACHE PATH "GPU kernel binary install location")
configure_file(config.h.in ${PROJECT_BINARY_DIR}/device/config.h)
set(headers ${PROJECT_BINARY_DIR}/device/config.h)

# ------------------------------------------------------------------------------------------------ #
# Build one of DACS_Device or GPU_Device
# ------------------------------------------------------------------------------------------------ #
if(HAVE_CUDA AND USE_CUDA)

  set(sources GPU_Device.cc)
  list(APPEND headers GPU_Device.hh device_cuda.h)
  set(include_dirs PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
                   $<BUILD_INTERFACE:${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}>)
  add_component_library(
    TARGET Lib_device
    TARGET_DEPS Lib_dsxx
    INCLUDE_DIRS "${include_dirs}"
    LIBRARY_NAME device
    LIBRARY_TYPE STATIC
    SOURCES "${sources}"
    HEADERS "${headers}")

endif()

install(FILES ${headers} DESTINATION ${DBSCFGDIR}include/device)

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

# Push some variables up one level
set(Draco_EXPORT_TARGET_PROPERTIES
    "${Draco_EXPORT_TARGET_PROPERTIES}"
    PARENT_SCOPE)

# ------------------------------------------------------------------------------------------------ #
# End devices/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
