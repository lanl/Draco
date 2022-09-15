#--------------------------------------------*-cmake-*---------------------------------------------#
# file   src/VendorsChecks/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2016 May 11
# brief  Instructions for building VendorsChecks/test level Makefile.
# note   Copyright (C) 2016-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(VendorsChecks_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #

if(TARGET METIS::metis)
  add_scalar_tests(SOURCES "tstMetis.cc" DEPS "Lib_dsxx;METIS::metis")
endif()

if(TARGET ParMETIS::parmetis AND ${DRACO_C4} STREQUAL "MPI")
  add_parallel_tests(
    SOURCES "tstParmetis.cc"
    DEPS "Lib_c4;ParMETIS::parmetis"
    PE_LIST "3")
endif()

# ------------------------------------------------------------------------------------------------ #
# end src/VendorsChecks/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
