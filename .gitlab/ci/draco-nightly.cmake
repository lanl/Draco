#--------------------------------------------*-cmake-*---------------------------------------------#
# file   .gitlab/ci/draco-nightly.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   Tuesday, Jun 02, 2020, 11:44 am
# brief  CTest regression script for DRACO
# note   Copyright (C) 2021 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #
# Ref: http://www.cmake.org/Wiki/CMake_Scripting_Of_CTest

set(CTEST_SOURCE_DIRECTORY "$ENV{DRACO_SOURCE_DIR}")
set(CTEST_BINARY_DIRECTORY "$ENV{DRACO_BINARY_DIR}")
set(CTEST_PROJECT_NAME "$ENV{PROJECT}")

set(CTEST_SUBMIT_URL "https://rtt.lanl.gov/cdash3/submit.php?project=$ENV{PROJECT}")
set(CTEST_CURL_OPTIONS CURLOPT_SSL_VERIFYPEER_OFF CURLOPT_SSL_VERIFYHOST_OFF)
set(CTEST_NIGHTLY_START_TIME "22:00:00 MDT")

if(WIN32)
  set(CTEST_CMAKE_GENERATOR "Visual Studio 16 2019")
else()
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
endif()
set(CTEST_SITE "$ENV{SITE_ID}")
if(DEFINED ENV{ARCH})
  string(APPEND CTEST_SITE "-$ENV{ARCH}")
endif()
string(APPEND CTEST_SITE "-gr") # indicating job was initiated by a gitlab runner.

set(CTEST_BUILD_NAME "$ENV{CTEST_BUILD_NAME}")
set(CTEST_MODE "$ENV{CTEST_MODE}")
set(CTEST_OUTPUT_ON_FAILURE ON)

cmake_host_system_information(RESULT MPI_PHYSICAL_CORES QUERY NUMBER_OF_PHYSICAL_CORES)
set(CMAKE_BUILD_PARALLEL_LEVEL ${MPI_PHYSICAL_CORES})
if(DEFINED ENV{CTEST_NPROC})
  set(CTEST_PARALLEL_LEVEL $ENV{CTEST_NPROC})
else()
  set(CTEST_PARALLEL_LEVEL ${MPI_PHYSICAL_CORES})
endif()
set(CTEST_BUILD_TYPE $ENV{CMAKE_BUILD_TYPE})
# if(CTEST_BUILD_TYPE STREQUAL MemCheck)
#   set(CTEST_BUILD_TYPE Debug)
# endif()
string(APPEND EXTRA_CMAKE_ARGS " $ENV{EXTRA_CMAKE_ARGS}")

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
set( CTEST_CONFIGURATION_TYPE "${CTEST_BUILD_CONFIGURATION}" ) # mayby only needed by msvc?

file(TO_CMAKE_PATH "${CI_PROJECT_DIR}" CI_PROJECT_DIR)
string(CONCAT CTEST_CONFIGURE_COMMAND "cmake"
              " -DCMAKE_INSTALL_PREFIX=$ENV{CI_PROJECT_DIR}/install" " ${EXTRA_CMAKE_ARGS}")
if(WIN32)
  file(TO_CMAKE_PATH "${CTEST_SOURCE_DIRECTORY}" CTEST_SOURCE_DIRECTORY)
  file(TO_CMAKE_PATH "$ENV{CMAKE_TOOLCHAIN_FILE}" CMAKE_TOOLCHAIN_FILE)
  string(APPEND CTEST_CONFIGURE_COMMAND " -A x64 -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}")
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
string(APPEND CTEST_CONFIGURE_COMMAND " -DCMAKE_VERBOSE_MAKEFILE=ON")
string(APPEND CTEST_CONFIGURE_COMMAND " ${CTEST_SOURCE_DIRECTORY}")
if(DEFINED ENV{CTEST_MEMORYCHECK_SUPPRESSIONS_FILE})
  if(EXISTS "$ENV{CTEST_MEMORYCHECK_SUPPRESSIONS_FILE}")
    set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE "$ENV{CTEST_MEMORYCHECK_SUPPRESSIONS_FILE}")
    set(MEMORYCHECK_SUPPRESSIONS_FILE "$ENV{CTEST_MEMORYCHECK_SUPPRESSIONS_FILE}")
  else()
    message(FATAL_ERROR "CTEST_MEMORYCHECK_SUPPRESSIONS_FILE set to invalid path.")
  endif()
endif()
message(
  "
CTEST_PARALLEL_LEVEL = ${CTEST_PARALLEL_LEVEL}
MPI_PHYSICAL_CORES   = ${MPI_PHYSICAL_CORES}

CTEST_CMAKE_GENERATOR     = ${CTEST_CMAKE_GENERATOR}
CTEST_BUILD_CONFIGURATION = ${CTEST_BUILD_CONFIGURATION}
CTEST_CONFIGURE_COMMAND   = ${CTEST_CONFIGURE_COMMAND}
")

set(CTEST_USE_LAUNCHERS FALSE)
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 50)
set(CTEST_UPDATE_COMMAND "git")
set(CTEST_GIT_UPDATE_CUSTOM "${CMAKE_COMMAND}" "-E" "echo" "Skipping git update (no-op).")

#message("Parsing ${CTEST_SCRIPT_DIRECTORY}/CTestCustom.cmake")
#ctest_read_custom_files("${CTEST_SCRIPT_DIRECTORY}")

# CTEST_CUSTOM_COVERAGE_EXCLUDE is a list of regular expressions. Any file name that matches any of
# the regular expressions in the list is excluded from the reported coverage data.
list(APPEND CTEST_CUSTOM_COVERAGE_EXCLUDE
  # don't report on actual unit tests
  "src/.*/test/"
  "src/.*/ftest/"
  # terminal isn't our code. don't report lack of coverage
  "src/.*/terminal/.*"
  "src/api/test_.*/"
  "src/.*/test_.*/"
  )

if(${CTEST_SCRIPT_ARG} MATCHES Configure)
  ctest_empty_binary_directory("${CTEST_BINARY_DIRECTORY}")
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

  message(
    "
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
    ctest_submit()
    message(FATAL_ERROR "configuration error")
  endif()

endif()

# ------------------------------------------------------------------------------------------------ #
# Build it
# ------------------------------------------------------------------------------------------------ #
if(${CTEST_SCRIPT_ARG} MATCHES Build)

  if(NOT WIN32)
    set(CTEST_BUILD_FLAGS "-j ${CTEST_PARALLEL_LEVEL} -l ${CMAKE_BUILD_PARALLEL_LEVEL}")
  endif()

  if(DEFINED ENV{AUTODOCDIR})
    # build one unit test
    message(
      "
ctest_build(
  TARGET Ut_dsxx_tstAssert_exe
  FLAGS $ENV{BUILD_FLAGS}
  RETURN_VALUE build_failure
  CAPTURE_CMAKE_ERROR ctest_build_errors)
")
    ctest_build(
      TARGET Ut_dsxx_tstAssert_exe
      FLAGS "$ENV{BUILD_FLAGS}"
      RETURN_VALUE build_failure
      CAPTURE_CMAKE_ERROR ctest_build_errors)
    # build autodoc target.
    message(
      "
ctest_build(
  TARGET autodoc
  FLAGS $ENV{BUILD_FLAGS}
  RETURN_VALUE build_failure
  CAPTURE_CMAKE_ERROR ctest_build_errors)
")
    ctest_build(
      TARGET autodoc
      FLAGS "$ENV{BUILD_FLAGS}"
      RETURN_VALUE build_failure
      CAPTURE_CMAKE_ERROR ctest_build_errors)
    if(build_failure)
      message("${ctest_build_errors}")
      ctest_submit()
      message(FATAL_ERROR "build error")
    endif()
  else()
    message(
      "
ctest_build(
  FLAGS $ENV{BUILD_FLAGS}
  RETURN_VALUE build_failure
  PARALLEL_LEVEL ${CTEST_PARALLEL_LEVEL}
  CAPTURE_CMAKE_ERROR ctest_build_errors)")
    ctest_build(
      FLAGS "$ENV{BUILD_FLAGS}"
      RETURN_VALUE build_failure
      # cmake-3.21 -- PARALLEL_LEVEL ${CMAKE_BUILD_PARALLEL_LEVEL}
      CAPTURE_CMAKE_ERROR ctest_build_errors)

    if(build_failure)
      message("${ctest_build_errors}")
      ctest_submit()
      message(FATAL_ERROR "build error")
    endif()
  endif()

endif()

# ------------------------------------------------------------------------------------------------ #
# Run the tests
# ------------------------------------------------------------------------------------------------ #
if(${CTEST_SCRIPT_ARG} MATCHES Test)

  if(DEFINED ENV{AUTODOCDIR})
    # Run one test to avoid missing data on CDash
    message("
ctest_test( RETURN_VALUE test_failure INCLUDE dsxx_tstAssert)
")
    ctest_test(RETURN_VALUE test_failure INCLUDE dsxx_tstAssert)
  else()
    if(DEFINED ENV{TEST_EXCLUSIONS})
      set(CTEST_TEST_EXTRAS EXCLUDE $ENV{TEST_EXCLUSIONS})
    endif()
    list(APPEND CTEST_TEST_EXTRAS PARALLEL_LEVEL ${CTEST_PARALLEL_LEVEL})
    message("
ctest_test( RETURN_VALUE test_failure ${CTEST_TEST_EXTRAS})
")
    ctest_test(RETURN_VALUE test_failure ${CTEST_TEST_EXTRAS})

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
      ctest_coverage_collect_gcov( TARBALL gcov.tgz
        TARBALL_COMPRESSION "FROM_EXT"
        QUIET
        GCOV_OPTIONS -b -p -l -x
        DELETE )

    endif()

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
        #if( DEFINED CTEST_MEMORYCHECK_SUPPRESSIONS_FILE )
        #  string(APPEND CTEST_MEMORYCHECK_COMMAND_OPTIONS
        #    " --suppressions=${CTEST_MEMORYCHECK_SUPPRESSIONS_FILE}")
        #endif()
      endif()

      set(CTEST_TEST_TIMEOUT 1200) # 1200 seconds = 20 minutes per test
      # message("ctest_memcheck( ${CTEST_TEST_EXTRAS} INCLUDE_LABEL memcheck"
      #         " EXCLUDE_LABEL nomemcheck)")
      # ctest_memcheck(${CTEST_TEST_EXTRAS} INCLUDE_LABEL memcheck EXCLUDE_LABEL nomemcheck)
      message("ctest_memcheck( ${CTEST_TEST_EXTRAS} EXCLUDE_LABEL nomemcheck)")
      ctest_memcheck(${CTEST_TEST_EXTRAS} EXCLUDE_LABEL nomemcheck)
    endif()
  endif()

endif()
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

  if(test_failure)
    message(FATAL_ERROR "test failure")
  endif()

endif()

# ------------------------------------------------------------------------------------------------ #
# end .gitlab/ci/draco-nightly.cmake
# ------------------------------------------------------------------------------------------------ #
