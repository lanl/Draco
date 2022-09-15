#--------------------------------------------*-cmake-*---------------------------------------------#
# file   experimental/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   Saturday, Dec 28, 2019, 15:35 pm
# brief  Provide future C++ features that can be made to work when limiting C++
#        standard to 14.
# note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(experimental CXX)

# These are header only implementations...

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
set(experimental_headers ${PROJECT_SOURCE_DIR}/mdspan)
file(GLOB supporting_headers_p0009_bits ${PROJECT_SOURCE_DIR}/__p0009_bits/*.hpp)

# ------------------------------------------------------------------------------------------------ #
# Installation instructions
# ------------------------------------------------------------------------------------------------ #
install(FILES ${experimental_headers} DESTINATION ${DBSCFGDIR}include/experimental)
install(FILES ${supporting_headers_p0009_bits}
        DESTINATION ${DBSCFGDIR}include/experimental/__p0009_bits)

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
# End experimental/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
