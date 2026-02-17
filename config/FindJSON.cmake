#--------------------------------------------*-cmake-*---------------------------------------------#
# file   draco/config/FindJSON.cmake
# author Kyle S. Beling
# brief  Helper for finding nolohmann::json.
# note   Copyright (C) 2025-2026 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

# * Find JSON Find the JSON includes
#
# JSON_INCLUDE_DIRS   - where to find json.hpp, etc. JSON_FOUND          - True if JSON found.

find_path(
  JSON_INCLUDE_DIR
  NAMES "nlohmann/json.hpp"
  PATHS "${NLOHMANN_JSON_INC_DIR}" "$ENV{NLOHMANN_JSON_INC_DIR}"
  PATH_SUFFIXES include)
mark_as_advanced(JSON_INCLUDE_DIR)

# handle the QUIETLY and REQUIRED arguments and set JSON_FOUND to TRUE if all listed variables are
# TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(JSON DEFAULT_MSG JSON_INCLUDE_DIR)

if(JSON_FOUND)
  set(JSON_FOUND
      ${JSON_FOUND}
      CACHE BOOL "Did we find the JSON include directory?")
  set(JSON_INCLUDE_DIRS ${JSON_INCLUDE_DIR})
endif()

if(VERBOSE)
  message(
    "
    JSON_FOUND        = ${JSON_FOUND}
    JSON_INCLUDE_DIRS = ${JSON_INCLUDE_DIR}
  ")
endif()
