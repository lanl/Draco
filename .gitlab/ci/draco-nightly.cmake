#--------------------------------------------*-cmake-*---------------------------------------------#
# file   .gitlab/ci/draco-nightly.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   Tuesday, Jun 02, 2020, 11:44 am
# brief  CTest regression script for DRACO
# note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #
# Ref: http://www.cmake.org/Wiki/CMake_Scripting_Of_CTest

# Avoid adding '-i' when using ninja
if(POLICY CMP0061)
  cmake_policy(SET CMP0061 NEW)
endif()

set(CTEST_SOURCE_DIRECTORY "$ENV{DRACO_SOURCE_DIR}")
set(CTEST_BINARY_DIRECTORY "$ENV{DRACO_BINARY_DIR}")
set(CTEST_PROJECT_NAME "$ENV{PROJECT}")

set(CTEST_SUBMIT_URL "https://rtt.lanl.gov/cdash3/submit.php?project=$ENV{PROJECT}")
set(CTEST_CURL_OPTIONS CURLOPT_SSL_VERIFYPEER_OFF CURLOPT_SSL_VERIFYHOST_OFF)
set(CTEST_NIGHTLY_START_TIME "22:00:00 MDT")

if(WIN32)
  set(CTEST_CMAKE_GENERATOR $ENV{GENERATOR})
else()
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
endif()
set(CTEST_UPDATE_COMMAND "git")
set(CTEST_GIT_UPDATE_CUSTOM "${CMAKE_COMMAND}" "-E" "echo" "Skipping git update (no-op).")

message(
"==> Starting a regression build...
CTEST_SOURCE_DIRECTORY = ${CTEST_SOURCE_DIRECTORY}
CTEST_BINARY_DIRECTORY = ${CTEST_BINARY_DIRECTORY}
CTEST_PROJECT_NAME     = ${CTEST_PROJECT_NAME}
CTEST_NIGHTLY_START_TIME = ${CTEST_NIGHTLY_START_TIME}
CTEST_CMAKE_GENERATOR  = ${CTEST_CMAKE_GENERATOR}
")

# ------------------------------------------------------------------------------------------------ #
# Collect information about this system
# - MPI_PHYSCIAL_CORES
# ------------------------------------------------------------------------------------------------ #
cmake_host_system_information(RESULT MPI_PHYSICAL_CORES QUERY NUMBER_OF_PHYSICAL_CORES)
# message("==> MPI_PHYSICAL_CORES = ${MPI_PHYSICAL_CORES}")
# message("==> MAXLOAD            = $ENV{MAXLOAD}")

if( (${MPI_PHYSICAL_CORES} LESS 2) AND ($ENV{MAXLOAD} GREATER 1))
  # for power9, the cmake command returns 1, so revert to the shell script value.
  set(MPI_PHYSICAL_CORES $ENV{MAXLOAD})
endif()
if(WIN32)
  # For Windows, just use the specified ${MAXLOAD}
  set(MPI_PHYSICAL_CORES $ENV{MAXLOAD})
endif()

# ------------------------------------------------------------------------------------------------ #
# Collect information specified by the *-jobs.yml job configuration
# - ARCH (optional)
# - CI_PROJECT_DIR
# - CTEST_BUILD_NAME
# - CTEST_BUILD_TYPE = {Debug, Release}
# - CTEST_MODE = {Nightly, Experimental}
# - CTEST_NPROC - Number of cores to use when running tests.
# - SITE_ID
#
# - EXTRA_CMAKE_ARGS
# - EXTRA_CTEST_ARGS
# - CTEST_MEMORYCHECK_SUPPRESSIONS_FILE
# ------------------------------------------------------------------------------------------------ #
file(TO_CMAKE_PATH "${CI_PROJECT_DIR}" CI_PROJECT_DIR)
set(CTEST_BUILD_NAME "$ENV{CTEST_BUILD_NAME}")
set(CTEST_BUILD_TYPE "$ENV{CMAKE_BUILD_TYPE}")
set(CTEST_MODE "$ENV{CTEST_MODE}")

set(CTEST_SITE "$ENV{SITE_ID}")
if(DEFINED ENV{ARCH})
  string(APPEND CTEST_SITE "-$ENV{ARCH}")
endif()
string(APPEND CTEST_SITE "-gr") # indicating job was initiated by a gitlab runner.

# separate_arguments(EXTRA_CMAKE_ARGS UNIX_COMMAND "$ENV{EXTRA_CMAKE_ARGS}") # leave as string
set(EXTRA_CMAKE_ARGS "$ENV{EXTRA_CMAKE_ARGS}")
separate_arguments(EXTRA_CTEST_ARGS UNIX_COMMAND "$ENV{EXTRA_CTEST_ARGS}")

if(DEFINED ENV{CTEST_MEMORYCHECK_SUPPRESSIONS_FILE})
  if(EXISTS "$ENV{CTEST_MEMORYCHECK_SUPPRESSIONS_FILE}")
    set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE "$ENV{CTEST_MEMORYCHECK_SUPPRESSIONS_FILE}")
    set(MEMORYCHECK_SUPPRESSIONS_FILE "$ENV{CTEST_MEMORYCHECK_SUPPRESSIONS_FILE}")
  else()
    message(FATAL_ERROR "CTEST_MEMORYCHECK_SUPPRESSIONS_FILE set to invalid path.")
  endif()
endif()

message(
  " ==> Job specification from yml
CI_PROJECT_DIR   = ${CI_PROJECT_DIR}
CTEST_BUILD_NAME = ${CTEST_BUILD_NAME}
CTEST_BUILD_TYPE = ${CTEST_BUILD_TYPE}
CTEST_MODE       = ${CTEST_MODE}
CTEST_SITE       = ${CTEST_SITE}
EXTRA_CMAKE_ARGS = ${EXTRA_CMAKE_ARGS}
EXTRA_CTEST_ARGS = ${EXTRA_CTEST_ARGS}
CTEST_NPROC      = $ENV{CTEST_NPROC}
MPI_PHYSICAL_CORES = ${MPI_PHYSICAL_CORES}
")

# ------------------------------------------------------------------------------------------------ #
# Options that control the build (but not set by yaml)
# ------------------------------------------------------------------------------------------------ #

# Verbosity controls
set(CTEST_OUTPUT_ON_FAILURE ON)
set(CTEST_USE_LAUNCHERS FALSE)
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 50)

# Build/test parallelism
set(CMAKE_BUILD_PARALLEL_LEVEL ${MPI_PHYSICAL_CORES})

if(DEFINED ENV{CTEST_NPROC})
  set(CTEST_PARALLEL_LEVEL $ENV{CTEST_NPROC})
else()
  set(CTEST_PARALLEL_LEVEL ${MPI_PHYSICAL_CORES})
endif()

string(TOUPPER ${CTEST_BUILD_TYPE} UPPER_CTEST_BUILD_TYPE)
if(${UPPER_CTEST_BUILD_TYPE} MATCHES DEBUG)
  set(CTEST_BUILD_CONFIGURATION "Debug")
elseif(${UPPER_CTEST_BUILD_TYPE} MATCHES RELWITHDEBINFO)
  set(CTEST_BUILD_CONFIGURATION "RelWithDebInfo")
elseif(${UPPER_CTEST_BUILD_TYPE} MATCHES MINSIZEREL)
  set(CTEST_BUILD_CONFIGURATION "MinSizeRel")
elseif(${UPPER_CTEST_BUILD_TYPE} MATCHES RELEASE)
  set(CTEST_BUILD_CONFIGURATION "Release")
else()
  message(FATAL_ERROR "CTEST_BUILD_TYPE (= ${UPPER_CTEST_BUILD_TYPE}) not set or not valid.")
endif()
set( CTEST_CONFIGURATION_TYPE "${CTEST_BUILD_CONFIGURATION}" )

# CTEST_CUSTOM_COVERAGE_EXCLUDE is a list of regular expressions. Any file name that matches any of
# the regular expressions in the list is excluded from the reported coverage data.
list(
  APPEND
  CTEST_CUSTOM_COVERAGE_EXCLUDE
  # don't report on actual unit tests
  "src/.*/test/"
  "src/.*/ftest/"
  # terminal isn't our code. don't report lack of coverage
  "src/.*/terminal/.*"
  "src/.*/test_.*/")

# ------------------------------------------------------------------------------------------------ #
# Start the build (clone, prep the build directory)
# ------------------------------------------------------------------------------------------------ #

# message("==> CTEST_SCRIPT_ARG = ${CTEST_SCRIPT_ARG}")

if(${CTEST_SCRIPT_ARG} MATCHES Configure)
  if( EXISTS "${CTEST_BINARY_DIRECTORY}")
    ctest_empty_binary_directory("${CTEST_BINARY_DIRECTORY}")
  endif()
  message("ctest_start( ${CTEST_MODE} )")
  ctest_start(${CTEST_MODE})
else()
  message("ctest_start( ${CTEST_MODE} APPEND )")
  ctest_start(${CTEST_MODE} APPEND)
endif()

ctest_update(SOURCE ${CTEST_SOURCE_DIRECTORY})

# ------------------------------------------------------------------------------------------------ #
# Configure the project with cmake.
# ------------------------------------------------------------------------------------------------ #
if(${CTEST_SCRIPT_ARG} MATCHES Configure)

  # create the comnand used to configure the build:

  file(TO_CMAKE_PATH $ENV{CI_PROJECT_DIR} CI_PROJECT_DIR)
  set(CTEST_CONFIGURE_COMMAND "cmake")
  #if(DEFINED ENV{GENERATOR})
  #  string(APPEND CTEST_CONFIGURE_COMMAND " -G $ENV{GENERATOR}")
  #endif()
  string(APPEND CTEST_CONFIGURE_COMMAND " -DCMAKE_INSTALL_PREFIX=${CI_PROJECT_DIR}/install"
    " ${EXTRA_CMAKE_ARGS}")
  if(WIN32)
    file(TO_CMAKE_PATH "${CTEST_SOURCE_DIRECTORY}" CTEST_SOURCE_DIRECTORY)
    file(TO_CMAKE_PATH "$ENV{CMAKE_TOOLCHAIN_FILE}" CMAKE_TOOLCHAIN_FILE)
    string(APPEND CTEST_CONFIGURE_COMMAND
       " -A x64 -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}")
    string(APPEND CTEST_CONFIGURE_COMMAND " -DDRACO_LIBRARY_TYPE=SHARED")
  else()
    string(APPEND CTEST_CONFIGURE_COMMAND " -DCMAKE_BUILD_TYPE=${UPPER_CTEST_BUILD_TYPE}")
  endif()
  if(DEFINED ENV{AUTODOCDIR})
    string(APPEND CTEST_CONFIGURE_COMMAND " -DAUTODOCDIR=$ENV{AUTODOCDIR}")
  endif()
  if(DEFINED ENV{CODECOV} AND "$ENV{CODECOV}" MATCHES "ON")
    set(CODE_COVERAGE "ON")
    string(APPEND CTEST_CONFIGURE_COMMAND " -DCODE_COVERAGE=ON")
  endif()
  # string(APPEND CTEST_CONFIGURE_COMMAND " -DCMAKE_VERBOSE_MAKEFILE=ON")
  string(APPEND CTEST_CONFIGURE_COMMAND " ${CTEST_SOURCE_DIRECTORY}")

  message(
    "
CTEST_CONFIGURE_COMMAND = ${CTEST_CONFIGURE_COMMAND}
ctest_configure(
  BUILD        ${CTEST_BINARY_DIRECTORY}
  RETURN_VALUE configure_failure
  CAPTURE_CMAKE_ERROR ctest_configure_errors)
")
  ctest_configure(
    BUILD "${CTEST_BINARY_DIRECTORY}"
    RETURN_VALUE configure_failure
    CAPTURE_CMAKE_ERROR ctest_configure_errors)

  if(configure_failure)
    message("${ctest_configure_errors}")
    if(${CTEST_SCRIPT_ARG} MATCHES Submit)
      ctest_submit()
    endif()
    message(FATAL_ERROR "configuration error")
  endif()

endif()

# ------------------------------------------------------------------------------------------------ #
# Build it
# ------------------------------------------------------------------------------------------------ #
if(${CTEST_SCRIPT_ARG} MATCHES Build)

  if(WIN32)
    set(CTEST_BUILD_FLAGS "-m:${MPI_PHYSICAL_CORES}")
  endif()

  if(DEFINED ENV{MAKEFILE_FLAGS})
    set(CTB_FLAGS FLAGS "$ENV{MAKEFILE_FLAGS}")
  endif()

  if(DEFINED ENV{AUTODOCDIR})
    # build one unit test
    message(
      "
ctest_build( ${CTB_FLAGS}
  TARGET Ut_dsxx_tstAssert_exe
  RETURN_VALUE build_failure
  CAPTURE_CMAKE_ERROR ctest_build_errors)
")
    ctest_build( ${CTB_FLAGS}
      TARGET Ut_dsxx_tstAssert_exe
      RETURN_VALUE build_failure
      CAPTURE_CMAKE_ERROR ctest_build_errors)
    # build autodoc target.
    message(
      "
ctest_build( ${CTB_FLAGS}
  TARGET autodoc
  RETURN_VALUE build_failure
  CAPTURE_CMAKE_ERROR ctest_build_errors)
")
    ctest_build( ${CTB_FLAGS}
      TARGET autodoc
      RETURN_VALUE build_failure
      CAPTURE_CMAKE_ERROR ctest_build_errors)
    if(build_failure OR (NOT "${ctest_build_errors}" STREQUAL "0") )
      message("${ctest_build_errors}")
      if(${CTEST_SCRIPT_ARG} MATCHES Submit)
        ctest_submit()
      endif()
      message(FATAL_ERROR "build error")
    endif()

  else()
    message(
      "
ctest_build( ${CTB_FLAGS}
  RETURN_VALUE build_failure
  CAPTURE_CMAKE_ERROR ctest_build_errors)")
    ctest_build( ${CTB_FLAGS}
      RETURN_VALUE build_failure
      # cmake-3.21 -- PARALLEL_LEVEL ${CMAKE_BUILD_PARALLEL_LEVEL}
      CAPTURE_CMAKE_ERROR ctest_build_errors)
    if(build_failure OR (NOT "${ctest_build_errors}" STREQUAL "0") )
      if(${CTEST_SCRIPT_ARG} MATCHES Submit)
        ctest_submit()
      endif()
      message(FATAL_ERROR "build error")
      file(WRITE "${CTEST_BINARY_DIRECTORY}/build_error.txt"
        "ctest_build_errors = ${ctest_build_errors}")
    endif()
  endif()
  unset(CTB_FLAGS)
endif()

# ------------------------------------------------------------------------------------------------ #
# Run the tests
# ------------------------------------------------------------------------------------------------ #
if(${CTEST_SCRIPT_ARG} MATCHES Test)

  if(DEFINED ENV{AUTODOCDIR})
    # Run one test to avoid missing data on CDash
    message("
ctest_test( RETURN_VALUE test_failure INCLUDE dsxx_tstAssert ${EXTRA_CTEST_ARGS})
")
    ctest_test(RETURN_VALUE test_failure INCLUDE dsxx_tstAssert ${EXTRA_CTEST_ARGS})
  else()
    list(APPEND EXTRA_CTEST_ARGS TEST_LOAD ${MPI_PHYSICAL_CORES})
    if(DEFINED ENV{TEST_EXCLUSIONS})
      list(APPEND EXTRA_CTEST_ARGS EXCLUDE $ENV{TEST_EXCLUSIONS})
    endif()
    list(APPEND EXTRA_CTEST_ARGS PARALLEL_LEVEL ${CTEST_PARALLEL_LEVEL})

    message("
ctest_test( RETURN_VALUE test_failure ${EXTRA_CTEST_ARGS})
")
    ctest_test(RETURN_VALUE test_failure ${EXTRA_CTEST_ARGS}
      CAPTURE_CMAKE_ERROR ctest_test_errors)

    if(test_failure OR (NOT "${ctest_test_errors}" STREQUAL "0") )
      if(${CTEST_SCRIPT_ARG} MATCHES Submit)
        ctest_submit()
      endif()
      message(FATAL_ERROR "test errors")
      file(WRITE "${CTEST_BINARY_DIRECTORY}/test_error.txt"
        "ctest_test_errors = ${ctest_test_errors}")
    endif()

    if(DEFINED ENV{CODECOV} AND "$ENV{CODECOV}" MATCHES "ON")
      set(CODE_COVERAGE "ON")

      # Count lines of code and upload results. CLOC commands:
      if( EXISTS "/ccs/codes/radtran/bin/cloc" )

        # Process and save lines of code
        message( "Generating lines of code statistics...
/ccs/codes/radtran/bin/cloc
        --exclude-list-file=${CTEST_SOURCE_DIRECTORY}/.gitlab/ci/cloc-exclude.cfg
        --exclude-lang=Text,Postscript
        --categorize=cloc-categorize.log
        --counted=cloc-counted.log
        --ignored=cloc-ignored.log
        --progress-rate=0
        --report-file=lines-of-code.log
        --force-lang-def=${CTEST_SOURCE_DIRECTORY}/.gitlab/ci/cloc-lang.defs
        ${CTEST_SOURCE_DIRECTORY}
            ")
        execute_process(
          COMMAND /ccs/codes/radtran/bin/cloc
          --exclude-list-file=${CTEST_SOURCE_DIRECTORY}/.gitlab/ci/cloc-exclude.cfg
          --exclude-lang=Text,Postscript
          --categorize=cloc-categorize.log
          --counted=cloc-counted.log
          --ignored=cloc-ignored.log
          --progress-rate=0
          --report-file=lines-of-code.log
          --force-lang-def=${CTEST_SOURCE_DIRECTORY}/.gitlab/ci/cloc-lang.defs
          ${CTEST_SOURCE_DIRECTORY}
          #  --3
          #  --diff
          WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
          )
        message( "Lines of code data at ${CTEST_SOURCE_DIRECTORY}/lines-of-code.log")
        message( "Generating lines of code statistics (omitting test directories)
/ccs/codes/radtran/bin/cloc
        --exclude-list-file=${CTEST_SOURCE_DIRECTORY}/.gitlab/ci/cloc-exclude.cfg
        --exclude-lang=Text,Postscript
        --categorize=cloc-categorize-notest.log
        --counted=cloc-counted-notest.log
        --ignored=cloc-ignored-notest.log
        --progress-rate=0
        --report-file=lines-of-code-notest.log
        --force-lang-def=${CTEST_SOURCE_DIRECTORY}/.gitlab/ci/cloc-lang.defs
        ${CTEST_SOURCE_DIRECTORY}
            ")
        execute_process(
          COMMAND /ccs/codes/radtran/bin/cloc
          --exclude-list-file=${CTEST_SOURCE_DIRECTORY}/.gitlab/ci/cloc-exclude.cfg
          --exclude-lang=Text,Postscript
          --categorize=cloc-categorize.log
          --counted=cloc-counted.log
          --ignored=cloc-ignored.log
          --progress-rate=0
          --report-file=lines-of-code-notest.log
          --force-lang-def=${CTEST_SOURCE_DIRECTORY}/.gitlab/ci/cloc-lang.defs
          ${CTEST_SOURCE_DIRECTORY}
          #  --3
          #  --diff
          WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
          )
        message( "Lines of code data at ${CTEST_SOURCE_DIRECTORY}/lines-of-code.log")
      endif()

      # GCOV/LCOV report
      find_program(CTEST_COVERAGE_COMMAND NAMES gcov)

      message("
CTEST_COVERAGE_COMMAND = ${CTEST_COVERAGE_COMMAND}
ctest_build(APPEND TARGET covrep)
")
      ctest_build(APPEND TARGET covrep)
      message("
ctest_coverage( BUILD \"${CTEST_BINARY_DIRECTORY}\" APPEND )
")
      ctest_coverage( BUILD "${CTEST_BINARY_DIRECTORY}" APPEND )
      message("ctest_coverage_collect_gcov( TARBALL gcov.tgz
        TARBALL_COMPRESSION \"FROM_EXT\"
        QUIET
        GCOV_OPTIONS -b -p -l -x
        DELETE )")
      include(CTestCoverageCollectGCOV)
      ctest_coverage_collect_gcov(
        TARBALL gcov.tgz
        TARBALL_COMPRESSION "FROM_EXT"
        QUIET
        GCOV_OPTIONS -b -p -l -x
        DELETE )

    endif() # CODECOV

    message("
==> ENV{MEMCHECK_CONFIGURATION} = $ENV{MEMCHECK_CONFIGURATION}
==> ENV{MEMORYCHECK_TYPE})      = $ENV{MEMORYCHECK_TYPE}
")
    if(DEFINED ENV{MEMCHECK_CONFIGURATION} AND ("$ENV{MEMCHECK_CONFIGURATION}" STREQUAL "ON"))
      if(DEFINED ENV{MEMORYCHECK_TYPE})
        set(CTEST_MEMORYCHECK_TYPE $ENV{MEMORYCHECK_TYPE})
      else()
        find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
        string(CONCAT CTEST_MEMORYCHECK_COMMAND_OPTIONS
          "-q --tool=memcheck --trace-children=yes --leak-check=full "
          "--show-reachable=yes --num-callers=20 --gen-suppressions=all")
        if(ENV{MEMCHECK_COMMAND_OPTIONS})
          string(APPEND CTEST_MEMORYCHECK_COMMAND_OPTIONS " $ENV{MEMCHECK_COMMAND_OPTIONS}")
        endif()
      endif()

      set(CTEST_TEST_TIMEOUT 1200) # 1200 seconds = 20 minutes per test
      message("ctest_memcheck( ${EXTRA_CTEST_ARGS})")
      ctest_memcheck(${EXTRA_CTEST_ARGS})
    endif()
  endif() # AUTODOC

endif() # Test

message("==> Done with testing")

# ------------------------------------------------------------------------------------------------ #
# Submit the results to CDash
# ------------------------------------------------------------------------------------------------ #
if(${CTEST_SCRIPT_ARG} MATCHES Submit)
  message("==> Trying to submit results to CDash")

  # files to look for
  set( datafiles
    "${CTEST_BINARY_DIRECTORY}/coverage.txt"
    "${CTEST_SOURCE_DIRECTORY}/lines-of-code.log"
    "${CTEST_SOURCE_DIRECTORY}/lines-of-code-notest.log"
    "${CTEST_SOURCE_DIRECTORY}/cloc-categorize.log"
    "${CTEST_SOURCE_DIRECTORY}/cloc-counted.log"
    "${CTEST_SOURCE_DIRECTORY}/cloc-ignored.log")
  foreach( file ${datafiles} )
    if(EXISTS "${file}")
      list(APPEND files_to_upload "${file}")
    endif()
  endforeach()
  message("ctest_upload( FILES ${files_to_upload} )" )
  ctest_upload( FILES ${files_to_upload} CAPTURE_CMAKE_ERROR err)
  message( "return code = ${err}" )
  if(EXISTS "${CTEST_BINARY_DIRECTORY}/gcov.tar")
    message("
==> ctest_submit(CDASH_UPLOAD \"${CTEST_BINARY_DIRECTORY}/gcov.tar\"
      CDASH_UPLOAD_TYPE GcovTar)
")
    ctest_submit(CDASH_UPLOAD "${CTEST_BINARY_DIRECTORY}/gcov.tar"
      CDASH_UPLOAD_TYPE GcovTar)
  endif()

  message("ctest_submit()")
  ctest_submit()

  if(EXISTS "${CTEST_BINARY_DIRECTORY}/build_error.txt" OR
      EXISTS "${CTEST_BINARY_DIRECTORY}/test_error.txt")
    message(FATAL_ERROR "build and/or test failures found")
  endif()

endif()

# ------------------------------------------------------------------------------------------------ #
# end .gitlab/ci/draco-nightly.cmake
# ------------------------------------------------------------------------------------------------ #
