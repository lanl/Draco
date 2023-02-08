# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   config/unix-intelllvm.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   Wednesday, May 26, 2021, 13:44 pm
# brief  Establish flags for Linux64 - Intel C++ (LLVM)
# note   Copyright (C) 2021-2023 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

# cmake-lint: disable=C0301

#
# Compiler Flags
#
# * https://www.intel.com/content/www/us/en/develop/documentation/oneapi-dpcpp-cpp-compiler-dev-guide-and-reference/top/compiler-reference/macros/use-predefined-macros-for-intel-compilers.html
# * https://community.intel.com/t5/Intel-C-Compiler/Cause-for-remark-11074-Inlining-inhibited-by-limit-max-size/m-p/1073635

if(NOT CXX_FLAGS_INITIALIZED)
  set(CXX_FLAGS_INITIALIZED
      "yes"
      CACHE INTERNAL "using draco settings.")

  string(APPEND CMAKE_C_FLAGS " -g -fp-model=precise -Wno-tautological-constant-compare")
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
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wno-undefined-var-template")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}")

  # Use C99 standard.
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}") # -std=c99

endif()

# ------------------------------------------------------------------------------------------------ #
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_C_FLAGS)
deduplicate_flags(CMAKE_CXX_FLAGS)

# Exceptions for -xHost
#
# * If this is a Cray, the compile wrappers take care of any xHost flags that are needed.
# * On ccs-net we mix & match cpu types (vendor libs are built for sandybridge) and -xHost causes an
#   ICE when linking libraries.
if(NOT CMAKE_CXX_COMPILER_WRAPPER STREQUAL CrayPrgEnv AND NOT IS_DIRECTORY "/ccs/opt/bin")
  set(HAS_XHOST TRUE)
  toggle_compiler_flag(HAS_XHOST "-xHost" "C;CXX" "")
endif()
toggle_compiler_flag(OpenMP_FOUND "${OpenMP_C_FLAGS}" "C;CXX" "")
force_compiler_flags_to_cache("C;CXX")

# ------------------------------------------------------------------------------------------------ #
# End config/unix-intelllvm.cmake
# ------------------------------------------------------------------------------------------------ #
