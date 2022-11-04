# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   config/Linux-NVHPC.cmake
# brief  Establish flags for Linux NVHPC
# note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)

# Note: In config/compilerEnv.cmake, the build system sets flags for
#
# 1. the language standard (C++14, C99, etc)
# 2. interprocedural optimization.

# Command line option information:
#
# * https://docs.nvidia.com/hpc-sdk/compilers/hpc-compilers-user-guide/index.html#cmdln-options-use
# * man nvc++

#
# Compiler Flags
#
if(NOT CXX_FLAGS_INITIALIZED)
  set(CXX_FLAGS_INITIALIZED
      "yes"
      CACHE INTERNAL "using draco settings.")

  string(APPEND CMAKE_C_FLAGS " -march=native -Mint128")
  # Use '--display_error_number' to find the warning numbers used by pragmas.
  set(CMAKE_C_FLAGS_DEBUG "-g -O0 -DDEBUG")
  string(CONCAT CMAKE_C_FLAGS_RELEASE "-gopt -O3 -Munroll=c:1 -Mnoframe -Mlre -Mautoinline"
                " -Mcache_align -Mflushz -DNDEBUG")
  set(CMAKE_C_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_RELEASE}")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g -O3")

  string(APPEND CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS}")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}")

endif()

# ------------------------------------------------------------------------------------------------ #
# Toggle for OpenMP support
if(OpenMP_C_FLAGS)
  toggle_compiler_flag(OpenMP_FOUND "${OpenMP_C_FLAGS}" "C" "")
endif()
if(OpenMP_CXX_FLAGS)
  toggle_compiler_flag(OpenMP_FOUND "${OpenMP_CXX_FLAGS}" "CXX" "")
endif()

# ------------------------------------------------------------------------------------------------ #
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_C_FLAGS)
deduplicate_flags(CMAKE_CXX_FLAGS)
force_compiler_flags_to_cache("C;CXX;EXE_LINKER")

# ------------------------------------------------------------------------------------------------ #
# End config/unix-clang.cmake
# ------------------------------------------------------------------------------------------------ #
