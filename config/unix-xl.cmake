# --------------------------------------------*-cmake-*------------------------------------------- #
# file   config/unix-xl.cmake
# author Gabriel Rockefeller, Kelly Thompson <kgt@lanl.gov>
# date   2012 Nov 1
# brief  Establish flags for Linux64 - IBM XL C++
# note   Copyright (C) 2012-2023 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)

# cmake-lint: disable=C0301,W0106
# ~~~
# Ref:https://www.ibm.com/support/knowledgecenter/en/SSXVZZ_16.1.1/com.ibm.xlcpp1611.lelinux.doc/compiler_ref/rucmpopt.html
# ~~~

#
# Compiler flag checks
#
include(platform_checks)

#
# Compiler Flags
#

if(NOT CXX_FLAGS_INITIALIZED)
  set(CXX_FLAGS_INITIALIZED
      "yes"
      CACHE INTERNAL "using draco settings.")

  string(APPEND CMAKE_C_FLAGS " -g -qflttrap -qmaxmem=-1")

  if(EXISTS /usr/gapps)
    string(APPEND CMAKE_C_FLAGS " --gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1") # ATS-2
  elseif(DEFINED ENV{CMAKE_CXX_COMPILER_CONFIG_FILE} AND EXISTS
                                                         "$ENV{CMAKE_CXX_COMPILER_CONFIG_FILE}")
    # Darwin power9
    set(CMAKE_CXX_COMPILER_CONFIG_FILE
        $ENV{CMAKE_CXX_COMPILER_CONFIG_FILE}
        CACHE FILEPATH "XL config file" FORCE)
    string(APPEND CMAKE_C_FLAGS " -F${CMAKE_CXX_COMPILER_CONFIG_FILE}")
  endif()

  set(CMAKE_C_FLAGS_DEBUG "-O0 -qfullpath -DDEBUG")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -qstrict=nans:operationprecision")
  set(CMAKE_C_FLAGS_RELEASE "-O2 -qstrict=nans:operationprecision -DNDEBUG")
  set(CMAKE_C_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_RELEASE}")

  if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "ppc64le")
    string(APPEND CMAKE_C_FLAGS " -qarch=pwr9 -qtune=pwr9")
  endif()

  # Email from Roy Musselman <roymuss@us.ibm.com, 2019-03-21: For C++14, add
  # -qxflag=disable__cplusplusOverride
  string(APPEND CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -qxflag=disable__cplusplusOverride"
         " -Wno-undefined-var-template")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}")

endif()

# ------------------------------------------------------------------------------------------------ #
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_C_FLAGS)
deduplicate_flags(CMAKE_CXX_FLAGS)
force_compiler_flags_to_cache("C;CXX")

# ------------------------------------------------------------------------------------------------ #
# End config/unix-xl.cmake
# ------------------------------------------------------------------------------------------------ #
