# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   .gitlab/ci/vs2022-toolchain.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   Tuesday, Jun 07, 2022, 13:35 pm
# brief  Toolchian settings for MSVC 2022 used by CI/regression
# note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

set(VSCMD_DEBUG 0)

# CSK support
list(APPEND CMAKE_PREFIX_PATH c:/work/vendors64/csk-0.6.0)

# Data
if(EXISTS "C:/work/vendors-ec/physical_data/eos")
  set(SESPATHU "C:/work/vendors-ec/physical_data/eos")
  set(SESAMEPATH "C:/work/vendors-ec/physical_data/eos")
endif()

# TPLs via vcpkg Use KT's fork of vcpkg to get R123 and eospac.
#
# * https://github.com/KineticTheory/vcpkg.
# * https://vcpkg.readthedocs.io/en/latest/about/faq/#can-i-use-my-own-cmake-toolchain-file-with
#   -vcpkgs-toolchain-file
#
if(EXISTS "C:/work/vcpkg22/scripts/buildsystems/vcpkg.cmake")
  # set(VCPKG_VERBOSE TRUE CACHE BOOL "opt")
  include("C:/work/vcpkg22/scripts/buildsystems/vcpkg.cmake")
  list(REMOVE_ITEM CMAKE_PREFIX_PATH "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/debug")
  set(ParMETIS_ROOT_DIR "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}")
  set(EOSPAC_ROOT_DIR "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}")
endif()

# Add sanity checks here?
