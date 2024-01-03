#--------------------------------------------*-cmake-*---------------------------------------------#
# file   cdi_ipcress/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for cdi_ipcress/test.
# note   Copyright (C) 2011-2024 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(cdi_ipcress_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #

set(test_sources tIpcressFile.cc tIpcressOpacity.cc tIpcressWithCDI.cc)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #

add_scalar_tests(SOURCES "${test_sources}" DEPS "Lib_cdi_ipcress")

include(ApplicationUnitTest)
# Test multi-group opacities found in an ipcressfile
add_app_unit_test(
  DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/tIpcress_Interpreter.py
  APP $<TARGET_FILE_DIR:Exe_Ipcress_Interpreter>/$<TARGET_FILE_NAME:Exe_Ipcress_Interpreter>
  TEST_ARGS "-v;-h;${PROJECT_SOURCE_DIR}/two-mats.ipcress;missingfile.ipcress"
  STDINFILE "${PROJECT_SOURCE_DIR}/IpcressInterpreter.stdin"
  GOLDFILE "${PROJECT_SOURCE_DIR}/Ipcress_Interpreter-two-mats.ipcress.gold"
  LABELS "nomemcheck")
# Test gray opacities found in an ipcressfile
add_app_unit_test(
  DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/tIpcress_Interpreter.py
  APP $<TARGET_FILE_DIR:Exe_Ipcress_Interpreter>/$<TARGET_FILE_NAME:Exe_Ipcress_Interpreter>
  TEST_ARGS "${PROJECT_SOURCE_DIR}/odfregression10.ipcress"
  STDINFILE "${PROJECT_SOURCE_DIR}/IpcressInterpreter.odf.stdin"
  GOLDFILE "${PROJECT_SOURCE_DIR}/Ipcress_Interpreter-odfregression10.ipcress.gold"
  LABELS "nomemcheck")

# ------------------------------------------------------------------------------------------------ #
# End cdi_ipcress/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
