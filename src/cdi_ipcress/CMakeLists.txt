#--------------------------------------------*-cmake-*---------------------------------------------#
# file   cdi_ipcress/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for cdi_ipcress.
# note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(cdi_ipcress CXX)
set(CMAKE_INSTALL_DO_STRIP NO)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
set(sources IpcressFile.cc IpcressDataTable.cc IpcressGrayOpacity.cc IpcressMultigroupOpacity.cc)
set(headers IpcressFile.hh IpcressMaterial.hh IpcressFile.t.hh IpcressDataTable.hh
            IpcressGrayOpacity.hh IpcressMultigroupOpacity.hh)

# ------------------------------------------------------------------------------------------------ #
# Build package library
# ------------------------------------------------------------------------------------------------ #
add_component_library(
  TARGET Lib_cdi_ipcress
  TARGET_DEPS "Lib_cdi"
  LIBRARY_NAME ${PROJECT_NAME}
  HEADERS "${headers}"
  SOURCES "${sources}")

add_component_executable(
  TARGET Exe_Ipcress_Interpreter
  TARGET_DEPS "Lib_cdi_ipcress"
  SOURCES ${PROJECT_SOURCE_DIR}/IpcressInterpreter.cc)

# ------------------------------------------------------------------------------------------------ #
# Installation instructions
# ------------------------------------------------------------------------------------------------ #
install(FILES ${headers} DESTINATION ${DBSCFGDIR}include/cdi_ipcress)
install(
  TARGETS Exe_Ipcress_Interpreter
  EXPORT draco-targets
  DESTINATION ${DBSCFGDIR}bin)
install(PROGRAMS python/ipcress_reader.py python/plot_ipcress_opacity.py
        DESTINATION ${DBSCFGDIR}bin)

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
# End cdi_ipcress/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
