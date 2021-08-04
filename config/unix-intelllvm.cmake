# -------------------------------------------*-cmake-*---------------------------------------------#
# file   config/unix-intelllvm.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   Wednesday, May 26, 2021, 13:44 pm
# brief  Establish flags for Linux64 - Intel C++ (LLVM)
# note   Copyright (C) 2021 Triad National Security, LLC., All rights reserved.
# -------------------------------------------------------------------------------------------------#

#
# Compiler Flags
#

if(NOT CXX_FLAGS_INITIALIZED)
  set(CXX_FLAGS_INITIALIZED
      "yes"
      CACHE INTERNAL "using draco settings.")

  string(APPEND CMAKE_C_FLAGS " -g")
  if(DBS_GENERATE_OBJECT_LIBRARIES)
    string(APPEND CMAKE_C_FLAGS " -ipo")
  endif()
  string(CONCAT CMAKE_C_FLAGS_DEBUG "-O0 -DDEBUG -Wno-potentially-evaluated-expression")
  # Useful for profilers
  #
  # * -gline-tables-only
  # * -fdebug-info-for-profiling
  string(CONCAT CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
  set(CMAKE_C_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_RELEASE}")
  string(CONCAT CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 ")

  string(APPEND CMAKE_CXX_FLAGS " ${CMAKE_C_FLAGS}")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}")

  # Use C99 standard.
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

endif()

# -------------------------------------------------------------------------------------------------#
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_C_FLAGS)
deduplicate_flags(CMAKE_CXX_FLAGS)
force_compiler_flags_to_cache("C;CXX")

# Exceptions for -xHost
#
# * If this is a Cray, the compile wrappers take care of any xHost flags that are needed.
# * On ccs-net we mix & match cpu types (vendor libs are built for sandybridge) and -xHost causes an
#   ICE when linking libraries.
if(NOT CMAKE_CXX_COMPILER_WRAPPER STREQUAL CrayPrgEnv AND NOT IS_DIRECTORY "/ccs/opt/bin")
  set(HAS_XHOST TRUE)
  toggle_compiler_flag(HAS_XHOST "-xHost" "C;CXX" "")
endif()
toggle_compiler_flag(OPENMP_FOUND ${OpenMP_C_FLAGS} "C;CXX" "")

# -------------------------------------------------------------------------------------------------#
# End config/unix-intelllvm.cmake
# -------------------------------------------------------------------------------------------------#
