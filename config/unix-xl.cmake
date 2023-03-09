# --------------------------------------------*-cmake-*------------------------------------------- #
# file   config/unix-xl.cmake
# author Gabriel Rockefeller, Kelly Thompson <kgt@lanl.gov>
# date   2012 Nov 1
# brief  Establish flags for Linux64 - IBM XL C++
# note   Copyright (C) 2012-2023 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)

# cmake-lint: disable=C0301
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
    # ATS-2
    string(APPEND CMAKE_C_FLAGS " --gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1")

  elseif(EXISTS /projects/opt/ppc64le/ibm AND NOT ("${CMAKE_CXX_FLAGS}" MATCHES "-F")) # Darwin
    # ~~~
    # Darwin power9 requires specialized config files for xlC + cuda builds...
    # ~~~
    # Only RHEL7 needs this special flag
    file(READ /etc/redhat-release rhr)
    string(REGEX REPLACE "[^0-9]*([0-9]+).([0-9]+).*" "\\1.\\2" redhat_version "${rhr}")
    if(redhat_version MATCHES "^7.([0-9]+)")
      string(REPLACE ":" ";" modules $ENV{LOADEDMODULES})
      # extract version from module environment.
      foreach(module ${modules})
        if(${module} MATCHES "^gcc")
          string(REGEX REPLACE "[^0-9]*([0-9]+).([0-9]+).([0-9]+)" "\\1.\\2.\\3" gcc_version
                               ${module})
        elseif(NOT DEFINED xlc_version AND (${module} MATCHES "^ibm/xlc" OR ${module} MATCHES "^xl"
                                           ))
          string(REGEX REPLACE "[^0-9]*([0-9]+).([0-9]+).([0-9]+).([0-9]+).*" "\\1.\\2.\\3.\\4"
                               xlc_version ${module})
          string(REGEX REPLACE "[^0-9]*([0-9]+).([0-9]+).([0-9]+).*" "\\1.\\2.\\3" xlc_version_3
                               ${xlc_version})
        elseif(${module} MATCHES "^cuda")
          if(NOT DEFINED cuda_version)
            string(REGEX REPLACE "[^0-9]*([0-9]+).([0-9]+).*" "\\1" cuda_version_major ${module})
            string(REGEX REPLACE "[^0-9]*([0-9]+).([0-9]+).*" "\\2" cuda_version_minor ${module})
            string(REGEX REPLACE "[^0-9]*([0-9]+).([0-9]+).([0-9]+)*.*" "\\3" cuda_version_patch
                                 ${module})
            if("${cuda_version_patch}" STREQUAL "${module}")
              set(cuda_version "${cuda_version_major}.${cuda_version_minor}")
            else()
              set(cuda_version "${cuda_version_major}.${cuda_version_minor}.${cuda_version_patch}")
            endif()
          endif()
        endif()
      endforeach()

      if(NOT DEFINED cuda_version)
        # if no cuda module is loaded, we still need to point to the config file that ends in
        # *.cuda.11.XX.
        set(cuda_version "11.0") # config files only exist for cuda-11.0
      endif()
      set(CMAKE_CXX_COMPILER_CONFIG_FILE "/projects/opt/ppc64le/ibm/")
      # Look for filename variations
      string(
        CONCAT CMAKE_CXX_COMPILER_CONFIG_FILE3 "${CMAKE_CXX_COMPILER_CONFIG_FILE}"
               "xlc-${xlc_version}/xlC/${xlc_version_3}/etc/xlc.cfg.rhel.${redhat_version}.gcc."
               "${gcc_version}.cuda.${cuda_version}")
      string(
        CONCAT CMAKE_CXX_COMPILER_CONFIG_FILE2 "${CMAKE_CXX_COMPILER_CONFIG_FILE}"
               "xlc-${xlc_version}/xlC/${xlc_version_3}/etc/xlc.cfg.rhel.${redhat_version}.gcc."
               "${gcc_version}.cuda.${cuda_version_major}.${cuda_version_minor}")
      if(EXISTS ${CMAKE_CXX_COMPILER_CONFIG_FILE3})
        set(CMAKE_CXX_COMPILER_CONFIG_FILE ${CMAKE_CXX_COMPILER_CONFIG_FILE3})
      else()
        set(CMAKE_CXX_COMPILER_CONFIG_FILE ${CMAKE_CXX_COMPILER_CONFIG_FILE2})
      endif()
      unset(CMAKE_CXX_COMPILER_CONFIG_FILE3)
      unset(CMAKE_CXX_COMPILER_CONFIG_FILE2)
      set(CMAKE_CXX_COMPILER_CONFIG_FILE
          ${CMAKE_CXX_COMPILER_CONFIG_FILE}
          CACHE FILEPATH "XL config file" FORCE)
      if(EXISTS ${CMAKE_CXX_COMPILER_CONFIG_FILE})
        string(APPEND CMAKE_C_FLAGS " -F${CMAKE_CXX_COMPILER_CONFIG_FILE}")
      else()
        message(
          FATAL_ERROR
            "IBM XLC selected (Darwin), but requested config file was not found."
            "\nCMAKE_CXX_COMPILER_CONFIG_FILE = ${CMAKE_CXX_COMPILER_CONFIG_FILE}"
            "\nTry setting CXX_FLAGS=-F<filepath> in your enironment.")
      endif()
      unset(gcc_version)
      unset(xl_version)
      unset(config_file)
    endif()
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
