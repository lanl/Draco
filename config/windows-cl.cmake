# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   windows-cl.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 June 5
# brief  Establish flags for Windows - MSVC
# note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)

if(NOT CXX_FLAGS_INITIALIZED)
  # Settings that are common for Visual Studio are in a shared file.  These settings are valid for
  # both cl.exe and clang-cl.exe.
  include(windows-msvc)

  # Notes on options:
  #
  # * /W[0-4] Warning levels. Draco is currently using /W4 (see src/CMakeLists.txt), but clients
  #   default to /W2. Use /W0 to disable all warnings.
  string(APPEND CMAKE_C_FLAGS " /W2")
  string(APPEND CMAKE_CXX_FLAGS " /W2")

  set(CXX_FLAGS_INITIALIZED
      "yes"
      CACHE INTERNAL "using Draco settings.")
endif()

# ------------------------------------------------------------------------------------------------ #
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_C_FLAGS)
deduplicate_flags(CMAKE_CXX_FLAGS)
set(DRACO_LINK_OPTIONS
    "${DRACO_LINK_OPTIONS}"
    CACHE STRING "link flags" FORCE)
force_compiler_flags_to_cache("C;CXX")

# ------------------------------------------------------------------------------------------------ #
# End windows-cl.cmake
# ------------------------------------------------------------------------------------------------ #
