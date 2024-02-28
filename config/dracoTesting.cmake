# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   config/compilerEnv.cmake
# brief  Default CMake build parameters
# note   Copyright (C) 2011-2024 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)

include(FeatureSummary)
option(BUILD_TESTING "Should we compile the tests?" ON)
add_feature_info(BUILD_TESTING BUILD_TESTING
                 "Turn off to prevent the compilation of unit tests (ctest).")

# enable ctest funcitons and run ctest in parallel if we have multiple cores.
if(BUILD_TESTING)
  include(CTest)
  enable_testing()
endif()

# ------------------------------------------------------------------------------------------------ #
# End dracoTesting.cmake
# ------------------------------------------------------------------------------------------------ #
