#--------------------------------------------*-cmake-*---------------------------------------------#
# file   cdi_ndi/test/CMakeLists.txt
# author Ben. R. Ryan <brryan@lanl.gov>
# date   2020 Feb 4
# brief  Generate build project files for cdi_ndi/test.
# note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(cdi_ndi_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
set(test_sources tstNDI_AtomicMass.cc tstNDI_TNReaction.cc tstNDI_CP_Eloss.cc)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #

# Both tests read a local test data file, 'ndi_data'. To avoid resource issues, tell ctest to avoid
# running these two test at the same time with the keyword RESOURCE_LOCK.
if(TARGET NDI::ndi)
  add_scalar_tests(SOURCES "${test_sources}" DEPS "Lib_cdi_ndi" RESOURCE_LOCK "ndi_data")
endif()

# ------------------------------------------------------------------------------------------------ #
# End cdi_ndi/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
