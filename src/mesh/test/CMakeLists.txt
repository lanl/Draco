#--------------------------------------------*-cmake-*---------------------------------------------#
# file   mesh/test/CMakeLists.txt
# author Ryan Wollaeger <wollaeger@lanl.gov>
# date   Thursday, Jun 07, 2018, 15:44 pm
# brief  Generate build project files for mesh/test.
# note   Copyright (C) 2018-2023 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(mesh_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #

# create unit test lists
set(one_pe_tests
    ${PROJECT_SOURCE_DIR}/tstDraco_Mesh.cc ${PROJECT_SOURCE_DIR}/tstDraco_Mesh_Builder.cc
    ${PROJECT_SOURCE_DIR}/tstX3D_Draco_Mesh_Reader.cc ${PROJECT_SOURCE_DIR}/tstCAMR_Mesh.cc)
set(two_pe_tests ${PROJECT_SOURCE_DIR}/tstDraco_Mesh_DD.cc ${PROJECT_SOURCE_DIR}/tstTopology.cc)
set(two_pe_only_tests ${PROJECT_SOURCE_DIR}/tstCAMR_Mesh_DD.cc)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #

set(libs_for_tests Lib_c4 Lib_mesh Lib_quadrature)

add_parallel_tests(
  SOURCES "${one_pe_tests}"
  PE_LIST "1"
  DEPS "${libs_for_tests}")

if(NOT ${DRACO_C4} STREQUAL "SCALAR")
  add_parallel_tests(
    SOURCES "${two_pe_tests}"
    PE_LIST "2;4"
    DEPS "${libs_for_tests}")
endif()

if(NOT ${DRACO_C4} STREQUAL "SCALAR")
  add_parallel_tests(
    SOURCES "${two_pe_only_tests}"
    PE_LIST "2"
    DEPS "${libs_for_tests}")
endif()

# ------------------------------------------------------------------------------------------------ #
# End mesh/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
