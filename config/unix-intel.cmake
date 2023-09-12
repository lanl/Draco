# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   config/unix-intel.cmake
# author Kelly Thompson
# date   2010 Nov 1
# brief  Establish flags for Linux64 - Intel C++
# note   Copyright (C) 2010-2023 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

# Compiler Flags

if(NOT CXX_FLAGS_INITIALIZED)
  set(CXX_FLAGS_INITIALIZED
      "yes"
      CACHE INTERNAL "using draco settings.")

  # * [KT 2015-07-10] I would like to turn on -w2, but this generates many warnings from Trilinos
  #   headers that I can't suppress easily (warning 191: type qualifier is meaningless on cast type)
  # * [KT 2015-07-10] -diag-disable 11060 -- disable warning that is issued when '-ip' is turned on
  #   and a library has no symbols (this occurs when capsaicin links some trilinos libraries.)
  string(APPEND CMAKE_C_FLAGS " -g -w1 -vec-report0 -diag-disable=remark -shared-intel -no-ftz"
         " -diag-disable=11060")
  if(DBS_GENERATE_OBJECT_LIBRARIES)
    string(APPEND CMAKE_C_FLAGS " -ipo")
  endif()
  string(CONCAT CMAKE_C_FLAGS_DEBUG "-O0 -inline-level=0 -ftrapuv -check=uninit -fp-model=precise"
                " -fp-speculation=safe -debug inline-debug-info -fno-omit-frame-pointer -DDEBUG")
  # -qno-opt-dynamic-align fixed LAP numerical sensitivity for Hypre (Gaber, 12/01/2020).
  string(CONCAT CMAKE_C_FLAGS_RELEASE "-O3 -fp-speculation=fast -fp-model=precise -pthread"
                " -qno-opt-dynamic-align -DNDEBUG")

  set(CMAKE_C_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_RELEASE}")
  string(CONCAT CMAKE_C_FLAGS_RELWITHDEBINFO "-debug inline-debug-info -O3 -pthread"
                " -fp-model=precise -fp-speculation=safe -fno-omit-frame-pointer")

  string(APPEND CMAKE_CXX_FLAGS " ${CMAKE_C_FLAGS}")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -early-template-check")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}")

  # Use C99 standard.
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

  # Disable FMA at the compile level if desired
  if(DEFINED FMA_NEVER_HARDWARE)
    string(APPEND CMAKE_C_FLAGS " -no-fma")
    string(APPEND CMAKE_CXX_FLAGS " -no-fma")
  else()
    string(APPEND CMAKE_C_FLAGS " -fma")
    string(APPEND CMAKE_CXX_FLAGS " -fma")
  endif()

endif()

find_library(INTEL_LIBM m)
mark_as_advanced(INTEL_LIBM)

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
force_compiler_flags_to_cache("C;CXX")

# ------------------------------------------------------------------------------------------------ #
# End config/unix-intel.cmake
# ------------------------------------------------------------------------------------------------ #
