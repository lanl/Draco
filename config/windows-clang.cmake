# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   windows-clang.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2020 March 27
# brief  Establish flags for Windows - MSVC+Clang
# note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)

if(NOT CXX_FLAGS_INITIALIZED)
  # Settings that are common for Visual Studio are in a shared file.  These settings are valid for
  # both cl.exe and clang-cl.exe.
  include(windows-msvc)

  # Notes on options:
  #
  # * /W[0-4] Warning levels. Draco is currently using /W4 (see src/CMakeLists.txt), but clients
  #   default to /W2 for cl.exe and /W2 for clang-cl.exe. /W0 suppresses all warnings.
  string(APPEND CMAKE_C_FLAGS " /W4")
  string(APPEND CMAKE_CXX_FLAGS " /W4")

  set(CXX_FLAGS_INITIALIZED
      "yes"
      CACHE INTERNAL "using Draco settings.")
endif()

# ------------------------------------------------------------------------------------------------ #
# Ensure cache values always match current selection

deduplicate_flags(CMAKE_C_FLAGS)
deduplicate_flags(CMAKE_CXX_FLAGS)
set(DRACO_LINK_OPTIONS
    ${DRACO_LINK_OPTIONS}
    CACHE STRING "link flags" FORCE)
force_compiler_flags_to_cache("C;CXX")

# ------------------------------------------------------------------------------------------------ #
# End windows-clang.cmake
# ------------------------------------------------------------------------------------------------ #
