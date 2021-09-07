#--------------------------------------------*-cmake-*---------------------------------------------#
# file   config/compilerEnv.cmake
# brief  Default CMake build parameters
# note   Copyright (C) 2010-2021 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#

include_guard(GLOBAL)
include(FeatureSummary)
include(string_manip)

if(NOT DEFINED PLATFORM_CHECK_OPENMP_DONE OR NOT DEFINED CCACHE_CHECK_AVAIL_DONE)
  message("
Compiler Setup...
")
endif()

# --------------------------------------------------------------------------------------------------#
# PAPI
# --------------------------------------------------------------------------------------------------#
if(DEFINED ENV{PAPI_HOME})
  set(HAVE_PAPI
      1
      CACHE BOOL "Is PAPI available on this machine?")
  set(PAPI_INCLUDE
      $ENV{PAPI_INCLUDE}
      CACHE PATH "PAPI headers at this location")
  set(PAPI_LIBRARY
      $ENV{PAPI_LIBDIR}/libpapi.so
      CACHE FILEPATH "PAPI library.")
endif()

if(HAVE_PAPI)
  set(PAPI_INCLUDE
      ${PAPI_INCLUDE}
      CACHE PATH "PAPI headers at this location")
  set(PAPI_LIBRARY
      ${PAPI_LIBDIR}/libpapi.so
      CACHE FILEPATH "PAPI library.")
  if(NOT EXISTS ${PAPI_LIBRARY})
    message(FATAL_ERROR "PAPI requested, but library not found. Set PAPI_LIBDIR to correct path.")
  endif()
  mark_as_advanced(PAPI_INCLUDE PAPI_LIBRARY)
  add_feature_info(HAVE_PAPI HAVE_PAPI "Provide PAPI hardware counters if available.")
endif()

# --------------------------------------------------------------------------------------------------#
# Query OpenMP availability
#
# This feature is usually compiler specific and a compile flag must be added. For this to work the
# <platform>-<compiler>.cmake files (e.g.:  unix-g++.cmake) call this macro.
# --------------------------------------------------------------------------------------------------#
macro(query_openmp_availability)
  if(NOT PLATFORM_CHECK_OPENMP_DONE)
    set(PLATFORM_CHECK_OPENMP_DONE
        TRUE
        CACHE BOOL "Is check for OpenMP done?")
    mark_as_advanced(PLATFORM_CHECK_OPENMP_DONE)
    message(STATUS "Looking for OpenMP...")
    find_package(OpenMP QUIET)
    if(OPENMP_FOUND)
      message(STATUS "Looking for OpenMP... ${OpenMP_C_FLAGS} (supporting the "
                     "${OpenMP_C_VERSION} standard)")
      if(OpenMP_C_VERSION VERSION_LESS 3.0)
        message(STATUS "OpenMP standard support is too old (< 3.0). Disabling OpenMP build "
                       "features.")
        set(OPENMP_FOUND FALSE)
        set(OpenMP_C_FLAGS
            ""
            CACHE BOOL "OpenMP disabled (too old)." FORCE)
      endif()
      set(OPENMP_FOUND
          ${OPENMP_FOUND}
          CACHE BOOL "Is OpenMP available?" FORCE)
    else()
      message(STATUS "Looking for OpenMP... not found")
    endif()
  endif()
endmacro()

# --------------------------------------------------------------------------------------------------#
# Force save compiler flags to CMakeCache.txt
# --------------------------------------------------------------------------------------------------#
function(force_compiler_flags_to_cache lang_list)
  if("${lang_list}x" STREQUAL "x")
    message(FATAL_ERROR "argument lang_list is required. eg. \"C;CXX\"")
  endif()
  foreach(lang ${lang_list})
    foreach(flag FLAGS FLAGS_DEBUG FLAGS_RELEASE FLAGS_MINSIZEREL FLAGS_RELWITHDEBINFO)
      set(CMAKE_${lang}_${flag}
          "${CMAKE_${lang}_${flag}}"
          CACHE STRING "compiler flags" FORCE)
    endforeach()
  endforeach()
  set(DRACO_LINK_OPTIONS
      "${DRACO_LINK_OPTIONS}"
      CACHE STRING "link flags" FORCE)
endfunction()

# --------------------------------------------------------------------------------------------------#
# De-duplicate compiler flags
#
# example: deduplicate_flags(CMAKE_C_FLAGS)
#
# * ${FLAGS} evaluates to a string like "CMAKE_C_FLAGS"
# * ${${FLAGS}} evaluates to a list of compiler options like "-Werror -O2"
# --------------------------------------------------------------------------------------------------#
function(deduplicate_flags FLAGS)
  set(flag_list ${${FLAGS}}) # ${FLAGS} is CMAKE_C_FLAGS, double ${${FLAGS}} is the string of flags.
  separate_arguments(flag_list)
  list(REMOVE_DUPLICATES flag_list)
  string(REGEX REPLACE "([^\\]|^);" "\\1 " _TMP_STR "${flag_list}")
  string(REGEX REPLACE "[\\](.)" "\\1" _TMP_STR "${_TMP_STR}") # fixes escaping
  set(${FLAGS}
      "${_TMP_STR}"
      PARENT_SCOPE)
endfunction()

# --------------------------------------------------------------------------------------------------#
# Setup compilers
# --------------------------------------------------------------------------------------------------#
macro(dbsSetupCompilers)

  if(NOT dbsSetupCompilers_done)
    set(dbsSetupCompilers_done "ON")

    # Bad platform
    if(NOT WIN32 AND NOT UNIX)
      message(FATAL_ERROR "Unsupported platform (not WIN32 and not UNIX ).")
    endif()

    # -------------------------------------------------------------------------------------------- #
    # Add user provided options:
    #
    # * Users may set environment variables
    #
    #   * C_FLAGS
    #   * CXX_FLAGS
    #   * Fortran_FLAGS
    #   * EXE_LINKER_FLAGS
    #
    # * Provide these as arguments to cmake as -DC_FLAGS="whatever".
    #
    # -------------------------------------------------------------------------------------------- #
    foreach(lang C CXX Fortran EXE_LINKER SHARED_LINKER CUDA)
      if(DEFINED ENV{${lang}_FLAGS})
        string(REPLACE "\"" "" tmp "$ENV{${lang}_FLAGS}")
        string(APPEND ${lang}_FLAGS " ${tmp}")
        unset(tmp)
      endif()
      if(DEFINED ENV{CMAKE_${lang}_FLAGS})
        string(REPLACE "\"" "" tmp "$ENV{CMAKE_${lang}_FLAGS}")
        string(APPEND CMAKE_${lang}_FLAGS " ${tmp}")
        unset(tmp)
      endif()
      if(${lang}_FLAGS)
        toggle_compiler_flag(TRUE "${${lang}_FLAGS}" ${lang} "")
      endif()
      if(CMAKE_${lang}_FLAGS)
        toggle_compiler_flag(TRUE "${CMAKE_${lang}_FLAGS}" ${lang} "")
      endif()
    endforeach()

    # Defaults for 1st pass:

    # shared or static libraries?
    if(${DRACO_LIBRARY_TYPE} MATCHES "STATIC")
      set(DRACO_SHARED_LIBS 0)
    elseif(${DRACO_LIBRARY_TYPE} MATCHES "SHARED")
      # This CPP symbol is used by config.h to signal if we are need to add declspec(dllimport) or
      # declspec(dllexport) for MSVC.
      set(DRACO_SHARED_LIBS 1)
    else()
      message(FATAL_ERROR "DRACO_LIBRARY_TYPE must be set to either STATIC or SHARED.")
    endif()
    set(DRACO_SHARED_LIBS
        "${DRACO_SHARED_LIBS}"
        CACHE STRING "Do we activate declspec(dllimport) or declspec(dllexport) for MSVC." FORCE)
    mark_as_advanced(DRACO_SHARED_LIBS)

    # -------------------------------------------------------------------------------------------- #
    # Setup common options for targets
    # -------------------------------------------------------------------------------------------- #

    # Control the use of interprocedural optimization. This used to be set by editing compiler flags
    # directly, but now that CMake has a universal toggle, we use it. This value is used in
    # component_macros.cmake when properties are assigned to individual targets.

    # See https://cmake.org/cmake/help/git-stage/policy/CMP0069.html
    if(WIN32)
      set(USE_IPO
          OFF
          CACHE BOOL "Enable Interprocedural Optimization for Release builds." FORCE)
    else()
      include(CheckIPOSupported)
      check_ipo_supported(RESULT USE_IPO)
    endif()

    # -------------------------------------------------------------------------------------------- #
    # Special build mode for Coverage (gcov+lcov+genthml)
    # https://github.com/codecov/example-cpp11-cmake
    # -------------------------------------------------------------------------------------------- #
    if(NOT TARGET coverage_config)
      add_library(coverage_config INTERFACE)
    endif()
    if(UNIX)
      option(CODE_COVERAGE "Enable coverage reporting" OFF)
      if(NOT CODE_COVERAGE)
        message(STATUS "Code coverage build ... disabled (CODE_COVERAGE=OFF)")
      endif()
      if(CODE_COVERAGE AND NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        message(STATUS "Code coverage build ... disabled (Compiler not GNU|Clang)")
      endif()
      if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")

        list(
          APPEND
          CODE_COVERAGE_IGNORE_REGEX
          /usr/*
          *test/*
          */opt/spack/*
          *terminal/*
          *FortranChecks/*)
        list(REMOVE_DUPLICATES CODE_COVERAGE_IGNORE_REGEX)
        set(CODE_COVERAGE_IGNORE_REGEX
            ${CODE_COVERAGE_IGNORE_REGEX}
            CACHE STRING "List of regex that lcov will ignore" FORCE)

        if(CMAKE_BUILD_TYPE STREQUAL Debug)
          # Add required flags (GCC & LLVM/Clang)
          target_compile_options(coverage_config INTERFACE --coverage)
          target_link_options(coverage_config INTERFACE --coverage)

          # If env variable is set use it, otherwise search for default name.
          find_program(LCOV NAMES "$ENV{LCOV}" lcov)
          find_program(GCOV NAMES "$ENV{GCOV}" gcov)
          if(EXISTS "${LCOV}" AND EXISTS "${GCOV}")
            # Add a custom target that prints the coverage report
            set(lcovopts1 --gcov-tool ${GCOV} --quiet)
            set(lcovopts2 ${lcovopts1} --output-file coverage.info)
            unset(lcov_ignore)
            foreach(myregex ${CODE_COVERAGE_IGNORE_REGEX})
              list(APPEND lcov_ignore '${myregex}')
            endforeach()
            if(EXISTS "${PROJECT_SOURCE_DIR}/config/capture_lcov.sh")
              set(captureLcov "${PROJECT_SOURCE_DIR}/config/capture_lcov.sh")
            elseif(EXISTS "${DRACO_DIR}/cmake/capture_lcov.sh")
              set(captureLcov "${DRACO_DIR}/cmake/capture_lcov.sh")
            else()
              message(
                FATAL_ERROR
                  "CODE_COVERAGE=ON, but required helper script capture_lcov.sh not found.  Looked "
                  "at ${PROJECT_SOURCE_DIR}/config/capture_lcov.sh and "
                  "${DRACO_DIR}/cmake/capture_lcov.sh")
            endif()
            add_custom_command(
              OUTPUT "${PROJECT_BINARY_DIR}/covrep_target_aways_out_of_date.txt"
              BYPRODUCTS "${PROJECT_BINARY_DIR}/coverage.info" "${PROJECT_BINARY_DIR}/coverage.txt"
              COMMAND ${LCOV} ${lcovopts2} --capture --directory .
              COMMAND ${LCOV} ${lcovopts2} --remove coverage.info ${lcov_ignore}
              COMMAND genhtml coverage.info --demangle-cpp --output-directory cov-html
              COMMAND "${captureLcov}" -g "${GCOV}" -l "${LCOV}")
            unset(captureLcov)
            add_custom_target(
              covrep
              DEPENDS "${PROJECT_BINARY_DIR}/covrep_target_aways_out_of_date.txt"
              COMMENT
                "
  ==> View text coverage report in file coverage.txt.
  ==> View HTML coverage report with command: firefox cov-html/index.html
  ==> Repeat text coverage report with command: lcov --list coverage.info
  ")
            message(STATUS "Code coverage build ... enabled ('make covrep' to see a text and/or a"
                           " html report)")
            block_indent(90 27 "CODE_COVERAGE_IGNORE_REGEX = ${CODE_COVERAGE_IGNORE_REGEX}")
          else() # EXISTS "${LCOV}" AND EXISTS "${GCOV}"
            message(STATUS "Code coverage build ... disabled (lcov and/or gcov not found)")
          endif() # EXISTS "${LCOV}" AND EXISTS "${GCOV}"

        else() # CMAKE_BUILD_TYPE STREQUAL Debug
          message(STATUS "Code coverage build ... disabled (CMAKE_BUILD_TYPE != Debug")
        endif() # CMAKE_BUILD_TYPE STREQUAL Debug
      endif() # CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang"
    endif(UNIX)

  endif(NOT dbsSetupCompilers_done)

endmacro()

# ------------------------------------------------------------------------------------------------ #
# Setup C++ Compiler
# ------------------------------------------------------------------------------------------------ #
macro(dbsSetupCxx)

  # Static or shared libraries? Set IPO options.
  dbssetupcompilers()

  # Do we have access to openMP?
  query_openmp_availability()

  # Deal with compiler wrappers
  if(${CMAKE_CXX_COMPILER} MATCHES "tau_cxx.sh")
    # When using the TAU profiling tool, the actual compiler vendor is hidden under the tau_cxx.sh
    # script.  Use the following command to determine the actual compiler flavor before setting
    # compiler flags (end of this macro).
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -tau:showcompiler OUTPUT_VARIABLE my_cxx_compiler)
  else()
    set(my_cxx_compiler ${CMAKE_CXX_COMPILER})
  endif()

  # Setup compiler flags
  get_filename_component(my_cxx_compiler "${my_cxx_compiler}" NAME)

  # If the CMake_<LANG>_COMPILER is a MPI wrapper...
  if("${my_cxx_compiler}" MATCHES "mpicxx")
    # MPI wrapper
    execute_process(COMMAND "${my_cxx_compiler}" --version
                    OUTPUT_VARIABLE mpicxx_version_output OUTPUT_STRIP_TRAILING_WHITESPACE)
    # make output safe for regex matching
    string(REPLACE "+" "x" mpicxx_version_output ${mpicxx_version_output})
    if(${mpicxx_version_output} MATCHES icpc)
      set(my_cxx_compiler icpc)
    endif()
  endif()

  # setup flags based on COMPILER_ID...
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI" OR "${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
    include(unix-pgi)
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "IntelLLVM" OR "${CMAKE_C_COMPILER_ID}" STREQUAL
                                                            "IntelLLVM")
    include(unix-intelllvm)
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" OR "${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    include(unix-intel)
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Cray" OR "${CMAKE_C_COMPILER_ID}" STREQUAL "Cray")
    include(unix-crayCC)
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
    if(APPLE)
      include(apple-clang)
    elseif(MSVC)
      include(windows-clang)
    else()
      include(unix-clang)
    endif()
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    include(unix-g++)
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC" OR "${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
    include(windows-cl)
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "XLClang" OR "${CMAKE_C_COMPILER_ID}" STREQUAL "XLCLang")
    include(unix-xl)
  else()
    # missing CMAKE_CXX_COMPILER_ID? - try to match the compiler path+name to a string.
    if("${my_cxx_compiler}" MATCHES "pgCC" OR "${my_cxx_compiler}" MATCHES "pgc[+][+]")
      include(unix-pgi)
    elseif("${my_cxx_compiler}" MATCHES "CC")
      message(
        FATAL_ERROR "I think the C++ compiler is a Cray compiler wrapper, but I don't know what "
                    "compiler is wrapped.  CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
    elseif("${my_cxx_compiler}" MATCHES "cl" AND WIN32)
      include(windows-cl)
    elseif("${my_cxx_compiler}" MATCHES "icpc")
      include(unix-intel)
    elseif("${my_cxx_compiler}" MATCHES "xl")
      include(unix-xl)
    elseif("${my_cxx_compiler}" MATCHES "clang" OR "${my_cxx_compiler}" MATCHES "llvm")
      if(APPLE)
        include(apple-clang)
      else()
        include(unix-clang)
      endif()
    elseif("${my_cxx_compiler}" MATCHES "[cg][+x]+")
      include(unix-g++)
    else()
      message(FATAL_ERROR "Build system does not support " "CXX=${my_cxx_compiler}, "
                          "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
    endif()
  endif()

  # To the greatest extent possible, installed versions of packages should record the configuration
  # options that were used when they were built.  For preprocessor macros, this is usually
  # accomplished via #define directives in config.h files.  A package's installed config.h file
  # serves as both a record of configuration options and a central location for macro definitions
  # that control features in the package.  Defining macros via the -D command-line option to the
  # preprocessor leaves no record of configuration choices (except in a build log, which may not be
  # preserved with the installation).
  #
  # Unfortunately, there are cases where a particular macro must be defined before some particular
  # system header file is included, or before any system header files are included.  In these
  # situations, using the config.h mechanism introduces sensitivity to the order of header files,
  # which can lead to brittleness; defining project-wide language- or system-feature macros via -D,
  # using CMake's add_definitions command, is an acceptable alternative.  Such definitions appear
  # below.

  if(NOT DEFINED CMAKE_REQUIRED_DEFINITIONS)
    set(CMAKE_REQUIRED_DEFINITIONS "")
  endif()

  # Enable the definition of UINT64_C in stdint.h (required by Random123).
  add_definitions(-D__STDC_CONSTANT_MACROS)
  set(CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS} -D__STDC_CONSTANT_MACROS")

  # Define _POSIX_C_SOURCE=200112 and _XOPEN_SOURCE=600, to enable definitions conforming to
  # POSIX.1-2001, POSIX.2, XPG4, SUSv2, SUSv3, and C99.  See the feature_test_macros(7) man page for
  # more information.
  add_definitions(-D_POSIX_C_SOURCE=200112 -D_XOPEN_SOURCE=600)
  set(CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS} -D_POSIX_C_SOURCE=200112")
  set(CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS} -D_XOPEN_SOURCE=600")
  if(APPLE)
    # Defining the above requires adding POSIX extensions, otherwise, include ordering still goes
    # wrong on Darwin, (i.e., putting fstream before iostream causes problems) see
    # https://code.google.com/p/wmii/issues/detail?id=89
    add_definitions(-D_DARWIN_C_SOURCE)
    set(CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS} -D_DARWIN_C_SOURCE ")
  endif()

  if(NOT CCACHE_CHECK_AVAIL_DONE)
    set(CCACHE_CHECK_AVAIL_DONE
        TRUE
        CACHE BOOL "Have we looked for ccache/f90cache?")
    mark_as_advanced(CCACHE_CHECK_AVAIL_DONE)
    # From https://crascit.com/2016/04/09/using-ccache-with-cmake/
    message(STATUS "Looking for ccache...")
    find_program(CCACHE_PROGRAM ccache)
    if(CCACHE_PROGRAM)
      message(STATUS "Looking for ccache... ${CCACHE_PROGRAM}")
      # Set up wrapper scripts
      set(CMAKE_C_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
      set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
      add_feature_info(CCache CCACHE_PROGRAM "Using ccache to speed up builds.")
    else()
      message(STATUS "Looking for ccache... not found.")
    endif()

    # From https://crascit.com/2016/04/09/using-ccache-with-cmake/
    message(STATUS "Looking for f90cache...")
    find_program(F90CACHE_PROGRAM f90cache)
    if(F90CACHE_PROGRAM)
      message(STATUS "Looking for f90cache... ${F90CACHE_PROGRAM}")
      set(CMAKE_Fortran_COMPILER_LAUNCHER "${F90CACHE_PROGRAM}")
      add_feature_info(F90Cache F90CACHE_PROGRAM "Using f90cache to speed up builds.")
    else()
      message(STATUS "Looking for f90cache... not found.")
    endif()

  endif()

endmacro()

# ------------------------------------------------------------------------------------------------ #
# Setup Static Analyzer (if any)
#
# Enable with: -DDRACO_STATIC_ANALYZER=[none|clang-tidy|iwyu|cppcheck|cpplint|iwyl]
#
# Default is 'none'
#
# Variables set by this macro
#
# * DRACO_STATIC_ANALYZER
# * CMAKE_CXX_CLANG_TIDY
# * CMAKE_CXX_INCLUDE_WHAT_YOU_USE
# * CMAKE_CXX_CPPCHECK
# * CMAKE_CXX_CPPLINT
# * CMAKE_CXX_LINK_WHAT_YOU_USE
#
# References:
#
# * https://blog.kitware.com/static-checks-with-cmake-cdash-iwyu-clang-tidy-lwyu-cpplint-and-cppcheck/
# * https://github.com/KratosMultiphysics/Kratos/wiki/How-to-use-Clang-Tidy-to-automatically-correct-code
# * https://www.kdab.com/clang-tidy-part-1-modernize-source-code-using-c11c14/
#
# ------------------------------------------------------------------------------------------------ #
macro(dbsSetupStaticAnalyzers)

  set(DRACO_STATIC_ANALYZER
      "none"
      CACHE STRING "Enable a static analysis tool")
  set_property(
    CACHE DRACO_STATIC_ANALYZER
    PROPERTY STRINGS
             "none"
             "clang-tidy"
             "iwyu"
             "cppcheck"
             "cpplint"
             "iwyl")

  # Sanity Checks
  if("${DRACO_STATIC_ANALYZER}" STREQUAL "clang-tidy" AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL
                                                          "Clang")
    message(FATAL_ERROR "When DRACO_STATIC_ANALYZER=clang-tidy, the CXX compiler must be clang.")
  endif()

  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

    # clang-tidy (https://clang.llvm.org/extra/clang-tidy/)
    #
    # Example:
    #
    # * cmake -DDRACO_STATIC_ANALYZER=clang-tidy -DCLANG_TIDY_CHECKS="-checks=generate-*"  ...
    if("${DRACO_STATIC_ANALYZER}" MATCHES "clang-tidy")
      if(NOT CMAKE_CXX_CLANG_TIDY)
        find_program(CMAKE_CXX_CLANG_TIDY clang-tidy)
        get_filename_component(CT_BPATH ${CMAKE_CXX_CLANG_TIDY} DIRECTORY)
        get_filename_component(CT_BPATH "${CT_BPATH}" DIRECTORY)
        string(CONCAT CLANG_TIDY_IPATH ${CT_BPATH} "/include/c++/v1")
        unset(CT_BPATH)
      endif()

      if(CMAKE_CXX_CLANG_TIDY)
        if(NOT CLANG_TIDY_OPTIONS)
          set(CLANG_TIDY_OPTIONS "-header-filter=.*[.]hh")
          # Only run clang tidy on src, test, examples and skip 3rd party libraries
          #
          # set(CLANG_TIDY_OPTIONS "\"-header-filter=.*\\b(src|test|examples)\\b\\/(?!lib).*\"")
        endif()
        set(CLANG_TIDY_OPTIONS
            "${CLANG_TIDY_OPTIONS}"
            CACHE STRING "clang-tidy extra options (eg: -header-filter=.*[.]hh;-fix)" FORCE)

        if(NOT CLANG_TIDY_CHECKS)
          # * -checks=mpi-*,bugprone-*,performance-*,modernize-*
          # * See full list: `clang-tidy -check=* -list-checks'
          # * Default: all modernize checks; except use-trailing-return-type.
          set(CLANG_TIDY_CHECKS "-checks=modernize-*,-modernize-use-trailing-return-type")
          if(DEFINED ENV{CI})
            string(REPLACE "-checks=" "--warnings-as-errors=" tmp ${CLANG_TIDY_CHECKS})
            string(APPEND CLANG_TIDY_CHECKS ";${tmp}")
          endif()
        endif()
        set(CLANG_TIDY_CHECKS
            "${CLANG_TIDY_CHECKS}"
            CACHE STRING "clang-tidy check options (eg: -checks=bugprone-*,mpi-*)" FORCE)

        set(CLANG_TIDY_IPATH
            "${CLANG_TIDY_IPATH}"
            CACHE STRING "clang-tidy extra include directories" FORCE)
        if(NOT "${CLANG_TIDY_CHECKS}" MATCHES "[-]checks[=]")
          message(FATAL_ERROR "Option CLANG_TIDY_CHECKS string must start with the string "
                              "'-check='")
        endif()
        # re-create clang-tidy command
        if("${CMAKE_CXX_CLANG_TIDY}" MATCHES "[-]checks[=]")
          list(GET CMAKE_CXX_CLANG_TIDY 0 CMAKE_CXX_CLANG_TIDY)
        endif()
        set(CMAKE_CXX_CLANG_TIDY
            "${CMAKE_CXX_CLANG_TIDY};${CLANG_TIDY_CHECKS};${CLANG_TIDY_OPTIONS}"
            CACHE STRING "Run clang-tidy on each source file before compile." FORCE)
      else()
        unset(CMAKE_CXX_CLANG_TIDY)
        unset(CMAKE_CXX_CLANG_TIDY CACHE)
      endif()
      # Sanity check
      if(NOT CLANG_TIDY_IPATH
         OR NOT CLANG_TIDY_CHECKS
         OR NOT CMAKE_CXX_CLANG_TIDY)
        message(
          FATAL_ERROR
            "clang-tidy mode requested but some required"
            " variables were not found:
           - CLANG_TIDY_IPATH     = ${CLANG_TIDY_IPATH}
           - CLANG_TIDY_CHECKS    = ${CLANG_TIDY_CHECKS}
           - CMAKE_CXX_CLANG_TIDY = ${CMAKE_CXX_CLANG_TIDY}")
      endif()
    endif()

    # include-what-you-use
    # https://github.com/include-what-you-use/include-what-you-use/blob/master/README.md
    if(${DRACO_STATIC_ANALYZER} STREQUAL "iwyu")
      find_program(CMAKE_CXX_INCLUDE_WHAT_YOU_USE iwyu)
      if(CMAKE_CXX_INCLUDE_WHAT_YOU_USE)
        if(NOT "${CMAKE_CXX_INCLUDE_WHAT_YOU_USE}" MATCHES "Xiwyu")
          set(CMAKE_CXX_INCLUDE_WHAT_YOU_USE
              "${CMAKE_CXX_INCLUDE_WHAT_YOU_USE};-Xiwyu;--transitive_includes_only"
              CACHE STRING "Run iwyu on each source file before compile." FORCE)
        endif()
      else()
        unset(CMAKE_CXX_INCLUDE_WHAT_YOU_USE)
        unset(CMAKE_CXX_INCLUDE_WHAT_YOU_USE CACHE)
      endif()
    endif()
  endif()

  # cppcheck http://cppcheck.sourceforge.net/ http://cppcheck.sourceforge.net/demo/
  if(${DRACO_STATIC_ANALYZER} STREQUAL "cppcheck")
    find_program(CMAKE_CXX_CPPCHECK cppcheck)
    if(CMAKE_CXX_CPPCHECK)
      if(NOT "${CMAKE_CXX_CPPCHECK}" MATCHES "-std=")
        set(CMAKE_CXX_CPPCHECK
            "${CMAKE_CXX_CPPCHECK};--std=c++14"
            CACHE STRING "Run cppcheck on each source file before compile." FORCE)
      endif()
    else()
      unset(CMAKE_CXX_CPPCHECK)
      unset(CMAKE_CXX_CPPCHECK CACHE)
    endif()
  endif()

  # cpplint https://github.com/cpplint/cpplint
  if(${DRACO_STATIC_ANALYZER} STREQUAL "cpplint")
    find_program(CMAKE_CXX_CPPLINT cpplint)
    if(CMAKE_CXX_CPPLINT)
      if(NOT "${CMAKE_CXX_CPPLINT}" MATCHES "linelength")
        set(CMAKE_CXX_CPPLINT
            "${CMAKE_CXX_CPPLINT};--linelength=81"
            CACHE STRING "Run cpplint on each source file before compile." FORCE)
      endif()
    else()
      unset(CMAKE_CXX_CPPLINT)
      unset(CMAKE_CXX_CPPLINT CACHE)
    endif()
  endif()

  # include-what-you-link
  # https://blog.kitware.com/static-checks-with-cmake-cdash-iwyu-clang-tidy-lwyu-cpplint-and-cppcheck
  if(${DRACO_STATIC_ANALYZER} MATCHES "iwyl" AND UNIX)
    option(CMAKE_LINK_WHAT_YOU_USE "Report if extra libraries are linked." TRUE)
  else()
    option(CMAKE_LINK_WHAT_YOU_USE "Report if extra libraries are linked." FALSE)
  endif()

  # Report

  if(NOT ${DRACO_STATIC_ANALYZER} STREQUAL "none")
    message("\nStatic Analyzer Setup...\n")

    if(NOT "${CMAKE_CXX_CLANG_TIDY}x" STREQUAL "x")
      message(STATUS "Enabling static analysis option: ${CMAKE_CXX_CLANG_TIDY}")
    endif()
    if(NOT "${CMAKE_CXX_INCLUDE_WHAT_YOU_USE}x" STREQUAL "x")
      message(STATUS "Enabling static analysis option: ${CMAKE_CXX_INCLUDE_WHAT_YOU_USE}")
    endif()
    if(NOT "${CMAKE_CXX_CPPCHECK}x" STREQUAL "x")
      message(STATUS "Enabling static analysis option: ${CMAKE_CXX_CPPCHECK}")
    endif()
    if(NOT "${CMAKE_CXX_CPPLINT}x" STREQUAL "x")
      message(STATUS "Enabling static analysis option: ${CMAKE_CXX_CPPLINT}")
    endif()
    if(CMAKE_LINK_WHAT_YOU_USE)
      message(STATUS "Enabling static analysis option: CMAKE_LINK_WHAT_YOU_USE")
    endif()
  endif()

endmacro()

# ------------------------------------------------------------------------------------------------ #
# Setup Fortran Compiler
#
# Use:
#
# * include( compilerEnv )
# * dbsSetupFortran( [QUIET] )
#
# Returns:
#
# * BUILD_SHARED_LIBS - bool
# * CMAKE_Fortran_COMPILER - fullpath
# * CMAKE_Fortran_FLAGS
# * CMAKE_Fortran_FLAGS_DEBUG
# * CMAKE_Fortran_FLAGS_RELEASE
# * CMAKE_Fortran_FLAGS_RELWITHDEBINFO
# * CMAKE_Fortran_FLAGS_MINSIZEREL
# * ENABLE_SINGLE_PRECISION - bool
# * DBS_FLOAT_PRECISION     - string (config.h)
# * PRECISION_DOUBLE | PRECISION_SINGLE - bool
#
# ------------------------------------------------------------------------------------------------ #
macro(dbsSetupFortran)

  dbssetupcompilers()

  # Toggle if we should try to build Fortran parts of the project.  This will be set to true if
  # $ENV{FC} points to a working compiler (e.g.: GNU or Intel compilers with Unix Makefiles) or if
  # the current project doesn't support Fortran but CMakeAddFortranSubdirectory can be used.
  option(HAVE_Fortran "Should we build Fortran parts of the project?" OFF)

  # Is Fortran enabled (it is considered 'optional' for draco)?
  get_property(_LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES)
  if(_LANGUAGES_ MATCHES Fortran)

    set(HAVE_Fortran ON)
    set(my_fc_compiler ${CMAKE_Fortran_COMPILER})

    # MPI wrapper
    if(${my_fc_compiler} MATCHES "mpif90")
      execute_process(COMMAND ${my_fc_compiler} --version
                      OUTPUT_VARIABLE mpifc_version_output OUTPUT_STRIP_TRAILING_WHITESPACE)
      if("${mpifc_version_output}" MATCHES "ifort")
        set(my_fc_compiler ifort)
      elseif("${mpifc_version_output}" MATCHES "GNU")
        set(my_fc_compiler gfortran)
      endif()
    endif()

    # setup flags
    if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "PGI")
      include(unix-pgf90)
    elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
      if("${CMAKE_Fortran_COMPILER}" MATCHES "ifx")
        include(unix-ifx)
      else()
        include(unix-ifort)
      endif()
    elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Cray")
      include(unix-crayftn)
    elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Clang")
      include(unix-flang)
    elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
      include(unix-gfortran)
    else()
      # missing CMAKE_Fortran_COMPILER_ID? - try to match the compiler path+name to a string.
      if(${my_fc_compiler} MATCHES "pgf9[05]" OR ${my_fc_compiler} MATCHES "pgfortran")
        include(unix-pgf90)
      elseif(${my_fc_compiler} MATCHES "ftn")
        message(
          FATAL_ERROR
            "I think the C++ compiler is a Cray compiler wrapper, but I don't know"
            " what compiler is wrapped. CMAKE_Fortran_COMPILER_ID = ${CMAKE_Fortran_COMPILER_ID}")
      elseif(${my_fc_compiler} MATCHES "ifort")
        include(unix-ifort)
      elseif(${my_fc_compiler} MATCHES "xl")
        include(unix-xlf)
      elseif(${my_fc_compiler} MATCHES "flang")
        include(unix-flang)
      elseif(${my_fc_compiler} MATCHES "gfortran")
        include(unix-gfortran)
      else()
        message(FATAL_ERROR "Build system does not support FC=${my_fc_compiler}")
      endif()
    endif()

    if(_LANGUAGES_ MATCHES "^C$" OR _LANGUAGES_ MATCHES CXX)
      include(FortranCInterface)
    endif()

  else()
    # If CMake doesn't know about a Fortran compiler, $ENV{FC}, then also look for a compiler to use
    # with CMakeAddFortranSubdirectory.
    message(STATUS "Looking for CMakeAddFortranSubdirectory Fortran " "compiler...")
    set(CAFS_Fortran_COMPILER "NOTFOUND")

    # Try to find a Fortran compiler (use MinGW gfortran for MSVC).
    find_program(
      CAFS_Fortran_COMPILER
      NAMES ${CAFS_Fortran_COMPILER} $ENV{CAFS_Fortran_COMPILER} gfortran
      PATHS c:/msys64/mingw64/bin c:/MinGW/bin c:/msys64/usr/bin)

    if(EXISTS ${CAFS_Fortran_COMPILER})
      set(HAVE_Fortran ON)
      message(STATUS "Looking for CMakeAddFortranSubdirectory Fortran "
                     "compiler... found ${CAFS_Fortran_COMPILER}")
    else()
      message(STATUS "Looking for CMakeAddFortranSubdirectory Fortran " "compiler... not found")
    endif()

  endif()

  set(HAVE_Fortran
      ${HAVE_Fortran}
      CACHE BOOL "Should we build Fortran portions of this project?" FORCE)

endmacro()

# ------------------------------------------------------------------------------------------------ #
# Setup Cuda Compiler
#
# Use:
#
# * include( compilerEnv )
# * dbsSetupCuda( [QUIET] )
#
# Helpers - these environment variables help cmake find/set CUDA envs.
#
# * ENV{CUDACXX}
# * ENV{CUDAFLAGS}
# * ENV{CUDAHOSTCXX}
#
# Returns:
#
# * BUILD_SHARED_LIBS - bool
# * CMAKE_CUDA_FLAGS
# * CMAKE_CUDA_FLAGS_DEBUG
# * CMAKE_CUDA_FLAGS_RELEASE
# * CMAKE_CUDA_FLAGS_RELWITHDEBINFO
# * CMAKE_CUDA_FLAGS_MINSIZEREL
#
# Notes:
#
# * https://devblogs.nvidia.com/tag/cuda/
# * https://devblogs.nvidia.com/building-cuda-applications-cmake/
#
# ------------------------------------------------------------------------------------------------ #
macro(dbsSetupCuda)

  # Toggle if we should try to build Cuda parts of the project. Will be enabled if ENV{CUDACXX} is
  # set.
  option(HAVE_CUDA "Should we build Cuda parts of the project?" OFF)

  # Is Cuda enabled (it is considered 'optional' for draco)?
  get_property(_LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES)
  if(_LANGUAGES_ MATCHES CUDA)
    # We found Cuda, keep track of this information.
    set(HAVE_CUDA ON)
    # User option to disable Cuda, even when it is available.
    option(USE_CUDA "Use Cuda?" ON)
  endif()

  # Save the results
  set(HAVE_CUDA
      ${HAVE_CUDA}
      CACHE BOOL "Should we build CUDA portions of this project?" FORCE)
  if(HAVE_CUDA AND USE_CUDA)
    # Use this string in 'project(foo ${CUDA_DBS_STRING})' commands to enable cuda per project.
    set(CUDA_DBS_STRING
        "CUDA"
        CACHE STRING "If CUDA is available, this variable is 'CUDA'")
    # Use this string as a toggle when calling add_component_library or add_scalar_tests to force
    # compiling with nvcc.
    set(COMPILE_WITH_CUDA LINK_LANGUAGE CUDA)

    # setup flags
    if("${CMAKE_CUDA_COMPILER_ID}" MATCHES "NVIDIA")
      include(unix-cuda)
    else()
      message(FATAL_ERROR "Build system does not support CUDACXX=${CMAKE_CUDA_COMPILER}")
    endif()
  endif()

  # Sanity Check
  if(USE_CUDA AND NOT HAVE_CUDA)
    message(FATAL_ERROR "==> Cuda requested but nvcc not found. Try loading a cuda module.")
  endif()

endmacro()

# ------------------------------------------------------------------------------------------------ #
# Setup profile tools: valgrind
#
# * cmake -DENABLE_MEMORYCHECK=ON \
#   -DCTEST_MEMORYCHECK_SUPPRESSIONS_FILE=/scratch/regress/ccsradregress/valgrind_suppress.txt ...
# * ctest -L memcheck
#
# ------------------------------------------------------------------------------------------------ #
function(dbsSetupProfilerTools)

  option(ENABLE_MEMORYCHECK "provide memorycheck tests" OFF)
  if(ENABLE_MEMORYCHECK)
    find_program(CMAKE_MEMORYCHECK_COMMAND valgrind)
    if(NOT CMAKE_MEMORYCHECK_COMMAND)
      message(WARNING "ENABLE_MEMPORYCHECK=ON, but valgrind not found. Disabling memorycheck "
                      "features.")
    endif()
    if(EXISTS "$ENV{CTEST_MEMORYCHECK_SUPPRESSIONS_FILE}")
      set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE "$ENV{CTEST_MEMORYCHECK_SUPPRESSIONS_FILE}")
    endif()
    if(EXISTS "${CTEST_MEMORYCHECK_SUPPRESSIONS_FILE}")
      set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE
          "${CTEST_MEMORYCHECK_SUPPRESSIONS_FILE}"
          CACHE STRING "valgrind suppressions file")
    endif()
    if(DEFINED ENV{CMAKE_MEMORYCHECK_COMMAND_OPTIONS})
      set(CMAKE_MEMORYCHECK_COMMAND_OPTIONS $ENV{CMAKE_MEMORYCHECK_COMMAND_OPTIONS})
    endif()
    if("${CMAKE_MEMORYCHECK_COMMAND_OPTIONS}notset" STREQUAL "notset")
      string(CONCAT CMAKE_MEMORYCHECK_COMMAND_OPTIONS
                    "-q --tool=memcheck --trace-children=yes --leak-check=full --num-callers=20 "
                    "--gen-suppressions=all ") # --show-reachable=yes
      if(EXISTS "${CTEST_MEMORYCHECK_SUPPRESSIONS_FILE}")
        string(APPEND CMAKE_MEMORYCHECK_COMMAND_OPTIONS
               " --suppressions=${CTEST_MEMORYCHECK_SUPPRESSIONS_FILE} ")
      endif()
    endif()
    set(CMAKE_MEMORYCHECK_COMMAND_OPTIONS
        "${CMAKE_MEMORYCHECK_COMMAND_OPTIONS}"
        CACHE STRING "valgrind options")
  endif()

endfunction()

# ------------------------------------------------------------------------------------------------ #
# Toggle a compiler flag based on a bool
#
# Examples: toggle_compiler_flag( GCC_ENABLE_ALL_WARNINGS "-Weffc++" "CXX" "DEBUG" )
# ------------------------------------------------------------------------------------------------ #
macro(toggle_compiler_flag switch compiler_flag compiler_flag_var_names build_modes)

  # generate names that are safe for CMake RegEx MATCHES commands
  string(REPLACE "+" "x" safe_compiler_flag ${compiler_flag})

  # Loop over types of variables to check: CMAKE_C_FLAGS, CMAKE_CXX_FLAGS, etc.
  foreach(comp ${compiler_flag_var_names})

    # sanity check
    if(NOT ${comp} STREQUAL "C"
       AND NOT ${comp} STREQUAL "CXX"
       AND NOT ${comp} STREQUAL "Fortran"
       AND NOT ${comp} STREQUAL "CUDA"
       AND NOT ${comp} STREQUAL "EXE_LINKER"
       AND NOT ${comp} STREQUAL "SHARED_LINKER")
      message(
        FATAL_ERROR
          "When calling " "toggle_compiler_flag(switch, compiler_flag, compiler_flag_var_names), "
          "compiler_flag_var_names must be set to one or more of these valid name: C;CXX;Fortran;"
          "CUDA;EXE_LINKER;SHARED_LINKER.")
    endif()

    string(REPLACE "+" "x" safe_CMAKE_${comp}_FLAGS "${CMAKE_${comp}_FLAGS}")

    if("${build_modes}x" STREQUAL "x") # set flags for all build modes

      if(${switch})
        if(NOT "${safe_CMAKE_${comp}_FLAGS}" MATCHES "${safe_compiler_flag}")
          set(CMAKE_${comp}_FLAGS
              "${CMAKE_${comp}_FLAGS} ${compiler_flag} "
              CACHE STRING "compiler flags" FORCE)
        endif()
      else()
        if("${safe_CMAKE_${comp}_FLAGS}" MATCHES "${safe_compiler_flag}")
          string(REPLACE "${compiler_flag}" "" CMAKE_${comp}_FLAGS ${CMAKE_${comp}_FLAGS})
          set(CMAKE_${comp}_FLAGS
              "${CMAKE_${comp}_FLAGS}"
              CACHE STRING "compiler flags" FORCE)
        endif()
      endif()

    else() # build_modes listed

      foreach(bm ${build_modes})

        string(REPLACE "+" "x" safe_CMAKE_${comp}_FLAGS_${bm} ${CMAKE_${comp}_FLAGS_${bm}})

        if(${switch})
          if(NOT "${safe_CMAKE_${comp}_FLAGS_${bm}}" MATCHES "${safe_compiler_flag}")
            set(CMAKE_${comp}_FLAGS_${bm}
                "${CMAKE_${comp}_FLAGS_${bm}} ${compiler_flag} "
                CACHE STRING "compiler flags" FORCE)
          endif()
        else()
          if("${safe_CMAKE_${comp}_FLAGS_${bm}}" MATCHES "${safe_compiler_flag}")
            string(REPLACE "${compiler_flag}" "" CMAKE_${comp}_FLAGS_${bm}
                           ${CMAKE_${comp}_FLAGS_${bm}})
            set(CMAKE_${comp}_FLAGS_${bm}
                "${CMAKE_${comp}_FLAGS_${bm}}"
                CACHE STRING "compiler flags" FORCE)
          endif()
        endif()

      endforeach()

    endif()

  endforeach()
endmacro()

# ------------------------------------------------------------------------------------------------ #
# End config/compiler_env.cmake
# ------------------------------------------------------------------------------------------------ #
