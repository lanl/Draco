#--------------------------------------------*-cmake-*---------------------------------------------#
# file   cdi_cpt/test/CMakeLists.txt
# author Ben. R. Ryan <brryan@lanl.gov>
# date   2024 Oct 9
# brief  Generate build project files for cdi_cpt/test.
# note   Copyright (C) 2025-2026 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(cdi_cpt_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
set(test_sources tstCPT_Collision_Source.cc)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #

if(TARGET shacl::physics)
  add_scalar_tests(SOURCES "${test_sources}" DEPS "Lib_cdi_cpt;avalanche::avalanche_test_resources")
endif()

# ------------------------------------------------------------------------------------------------ #
# End cdi_cpt/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
