#--------------------------------------------*-cmake-*---------------------------------------------#
#
# file   predict/test/CMakeLists.txt
# author Mathew Cleveland
# date   Jan. 12 2020
# brief  CMake configuration for nnetwork_interface tests
# note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
#
#--------------------------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.21.0)
project(predict_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #

file(GLOB test_sources *.cc)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #

add_parallel_tests(
  SOURCES "${test_sources}"
  DEPS Lib_predict
  PE_LIST "1;3")

# Provide libtorch network
file(GLOB nn_files ${PROJECT_SOURCE_DIR}/*.pt ${PROJECT_SOURCE_DIR}/*.nnb)
provide_aux_files(FILES "${nn_files}" FOLDER predict_test)

# ------------------------------------------------------------------------------------------------ #
# end kde/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
