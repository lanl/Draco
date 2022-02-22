#--------------------------------------------*-cmake-*---------------------------------------------#
# file   diagnostics/qt/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2015 June 26
# brief  Generate build project files for diagnostics/qt
# note   Copyright (C) 2014-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.17.0)
project(diagnostics_qt CXX)

# ------------------------------------------------------------------------------------------------ #
# Special options for Qt applications
# ------------------------------------------------------------------------------------------------ #

# Instruct CMake to run MOC automatically when needed (only for subdirectories that need Qt)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #

set(sources main.cc mainwindow.cc diWidget.cc)
set(headers mainwindow.hh diWidget.hh)
set(ui_files mainwindow.ui)
# file( GLOB resources *.qrc )

# ------------------------------------------------------------------------------------------------ #
# Build package library
# ------------------------------------------------------------------------------------------------ #

# List Qt5::Widgets as a dependencies but not a VENDOR_LIB to prevent DBS from listing this library
# as a DRACO_TPL_LIBRARY.  It is not intended for export into other projects.
add_component_executable(
  TARGET Exe_draco_info_gui
  TARGET_DEPS "Lib_diagnostics;Qt5::Widgets"
  SOURCES "${sources};${ui_files}"
  FOLDER diagnostics NOCOMMANDWINDOW)

# Disable IPO for applications that link to Qt SDK, as this causes compile warnings for the Intel
# compiler (and it refuses to run IPO for this application anyway).
if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set_target_properties(Exe_draco_info_gui PROPERTIES INTERPROCEDURAL_OPTIMIZATION_RELEASE OFF)
endif()

# Copy necessary dll files to build directory
if(WIN32)
  add_custom_command(
    TARGET Exe_draco_info_gui
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:Lib_dsxx>
            $<TARGET_FILE_DIR:Exe_draco_info_gui>
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:Lib_c4>
            $<TARGET_FILE_DIR:Exe_draco_info_gui>
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:Lib_diagnostics>
            $<TARGET_FILE_DIR:Exe_draco_info_gui>
    COMMENT "")
endif()

# ------------------------------------------------------------------------------------------------ #
# end diagnostics/qt/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
