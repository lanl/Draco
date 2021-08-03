# --------------------------------------------*-cmake-*------------------------------------------- #
# file   config/unix-xl.cmake
# author Gabriel Rockefeller, Kelly Thompson <kgt@lanl.gov>
# date   2012 Nov 1
# brief  Establish flags for Linux64 - IBM XL C++
# note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)

# Ref:
# https://www.ibm.com/support/knowledgecenter/en/SSXVZZ_16.1.1/com.ibm.xlcpp1611.lelinux.doc/compiler_ref/rucmpopt.html
#
# Compiler flag checks
#
include(platform_checks)
query_openmp_availability()

#
# Compiler Flags
#

if(NOT CXX_FLAGS_INITIALIZED)
  set(CXX_FLAGS_INITIALIZED
      "yes"
      CACHE INTERNAL "using draco settings.")

  string(APPEND CMAKE_C_FLAGS " -g -qflttrap -qmaxmem=-1")
  if(EXISTS /usr/gapps)
    # ATS-2
    string(APPEND CMAKE_C_FLAGS " --gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1")
  elseif(EXISTS /projects/opt/ppc64le/ibm)
    # Darwin power9 - extract version from module environment.
    string(REPLACE ":" ";" modules $ENV{LOADEDMODULES})
    foreach(module ${modules})
      if(${module} MATCHES "^gcc")
        string(REGEX REPLACE "[^0-9]*([0-9]+).([0-9]+).([0-9]+)" "\\1.\\2.\\3" gcc_version
                             ${module})
      elseif(NOT DEFINED xlc_version AND ${module} MATCHES "^ibm/xlc")
        string(REGEX REPLACE "[^0-9]*([0-9]+).([0-9]+).([0-9]+).([0-9]+).*" "\\1.\\2.\\3.\\4"
                             xlc_version ${module})
        string(REGEX REPLACE "[^0-9]*([0-9]+).([0-9]+).([0-9]+).*" "\\1.\\2.\\3" xlc_version_3
                             ${xlc_version})
      elseif(${module} MATCHES "^cuda")
        if(NOT DEFINED cuda_version)
          string(REGEX REPLACE "[^0-9]*([0-9]+).([0-9]+)" "\\1.\\2" cuda_version ${module})
        endif()
      endif()
    endforeach()
    # Only redhat 7.7 is currently supported
    file(READ /etc/redhat-release rhr)
    string(REGEX REPLACE "[^0-9]*([0-9]+).([0-9]+).*" "\\1.\\2" redhat_version "${rhr}")
    if(NOT DEFINED cuda_version)
      # if no cuda module is loaded, we still need to point to the config file that ends in
      # *.cuda.11.0.
      set(cuda_version "11.0")
    endif()
    string(
      CONCAT CMAKE_CXX_COMPILER_CONFIG_FILE
             "/projects/opt/ppc64le/ibm/xlc-${xlc_version}/xlC/${xlc_version_3}/etc/xlc.cfg.rhel."
             "${redhat_version}.gcc.${gcc_version}.cuda.${cuda_version}")
    set(CMAKE_CXX_COMPILER_CONFIG_FILE
        ${CMAKE_CXX_COMPILER_CONFIG_FILE}
        CACHE FILEPATH "XL config file" FORCE)
    if(EXISTS ${CMAKE_CXX_COMPILER_CONFIG_FILE})
      string(APPEND CMAKE_C_FLAGS " -F${CMAKE_CXX_COMPILER_CONFIG_FILE}")
    else()
      message(FATAL_ERROR "IBM XLC selected (Darwin), but requested config file was not found."
                          "\nCMAKE_CXX_COMPILER_CONFIG_FILE = ${CMAKE_CXX_COMPILER_CONFIG_FILE}")
    endif()
    unset(gcc_version)
    unset(xl_version)
    unset(config_file)
  endif()

  set(CMAKE_C_FLAGS_DEBUG "-O0 -qsmp=omp:noopt -qfullpath -DDEBUG")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -qsmp=omp -qstrict=nans:operationprecision")
  set(CMAKE_C_FLAGS_RELEASE "-O2 -qsmp=omp -qstrict=nans:operationprecision -DNDEBUG")
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

# Toggle for OpenMP support
if(OpenMP_C_FLAGS)
  toggle_compiler_flag(OPENMP_FOUND "${OpenMP_C_FLAGS}" "C" "")
endif()
if(OpenMP_CXX_FLAGS)
  toggle_compiler_flag(OPENMP_FOUND "${OpenMP_CXX_FLAGS}" "CXX" "")
endif()

# CMake will set OpenMP_C_FLAGS to '-qsmp.'  This option turns on OpenMP but also activates the
# auto-parallelizer.  We don't want to enable the 2nd feature so we need to specify the OpenMP flag
# to be '-qsmp=omp.'
if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13.0)
  toggle_compiler_flag(OPENMP_FOUND "-qsmp=omp" "C;CXX;EXE_LINKER" "")
endif()

# ------------------------------------------------------------------------------------------------ #
# End config/unix-xl.cmake
# ------------------------------------------------------------------------------------------------ #
