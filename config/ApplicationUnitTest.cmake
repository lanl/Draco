#--------------------------------------------*-cmake-*---------------------------------------------#
# file  config/ApplicationUnitTest.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date  Monday, Nov 19, 2012, 16:21 pm
# brief Provide macros that aid in creating unit tests that run interactive user codes (i.e.: run a
#       binary that reads an input file and diff the resulting output file).
# note  Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#

# * Reference: https://re-git.lanl.gov/draco/draco/-/wikis/CMake-based_ApplicationUnitTest
# * Example: draco/src/diagnostics/test/tDracoInfo.cmake (and associated CMakeLists.txt).

# ------------------------------------------------------------------------------------------------ #
# Example from draco/src/diagnostics/test/tDracoInfo.cmake
# ------------------------------------------------------------------------------------------------ #

# Use config/ApplicationUnitTest.cmake test registration:
#
# ~~~
# include( ApplicationUnitTest )
# add_app_unit_test(
#    DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/tDracoInfo.cmake
#    APP $<TARGET_FILE_DIR:Exe_draco_info>/$<TARGET_FILE_NAME:Exe_draco_info>
#    LABELS nomemcheck )
# ~~~
#
# Optional Parameters:
#
# * GOLDFILE      - Compare the output from APP against this file.
# * STDINFILE     - APP expects interactive input, use data from this file.
# * WORKDIR       - APP must be run from this directory.
# * BUILDENV
# * FAIL_REGEX
# * LABELS        - E.g.: nomemcheck, nr, perfbench
# * PASS_REGEX
# * PE_LIST       - How may mpi ranks to use "1;2;4"
# * RESOURCE_LOCK - Prevent tests with the same string from running concurrently.
# * RUN_AFTER     - Run the named tests before this test is started.
# * TEST_ARGS     - optional papmeters that will be given to APP.

# The above will generate a test with data similar to this:
#
# ~~~
# add_test(
#   NAME ${ctestname_base}${argname}
#   COMMAND ${CMAKE_COMMAND}
#   -D APP=${aut_APP}
#   -D ARGVALUE=${argvalue}
#   -D WORKDIR=${aut_WORKDIR}
#   -D TESTNAME=${ctestname_base}${argname}
#   -D DRACO_CONFIG_DIR=${DRACO_CONFIG_DIR}
#   -D DRACO_INFO=$<TARGET_FILE_DIR:Exe_draco_info>/$<TARGET_FILE_NAME:Exe_draco_info>
#   -D STDINFILE=${aut_STDINFILE}
#   -D GOLDFILE=${aut_GOLDFILE}
#   -D RUN_CMD=${RUN_CMD}
#   -D numPE=${numPE}
#   -D PROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
#   -D PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
#   ${BUILDENV}
#   -P ${aut_DRIVER}
# )
# set_tests_properties( diagnostics_draco_info
#    PROPERTIES
#      PASS_REGULAR_EXPRESSION Passes
#      FAIL_REGULAR_EXPRESSION Fails
#      LABELS nomemcheck
#  )
# ~~~
# Variables defined above can be used in this script.

# ------------------------------------------------------------------------------------------------ #

# cmake-lint: disable=R0912,R0915

include_guard(GLOBAL)
set(VERBOSE_DEBUG OFF)

# String manipulation
function(JOIN VALUES GLUE OUTPUT)
  string(REGEX REPLACE "([^\\]|^);" "\\1${GLUE}" _TMP_STR "${VALUES}")
  string(REGEX REPLACE "[\\](.)" "\\1" _TMP_STR "${_TMP_STR}") # fixes escaping
  set(${OUTPUT}
      "${_TMP_STR}"
      PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------------------------ #
# Check values for $APP
#
# Requires: APP        - name of executable that will be run
#
# Returns:
#
# * APP       - cleaned up path to executable
# * STDINFILE - cleaned up path to intput file.
# * BINDIR    - Directory location of binary file
# * PROJECT_BINARY_DIR - Parent directory of BINDIR
# * GOLDFILE  - cleaned up path to gold standard file.
# * OUTFILE   - Output filename derived from GOLDFILE.
# * ERRFILE   - Output(error) filename derived from GOLDFILE.
#
# if VERBOSE is set, also echo input values.
#
# * APP       - path name for executable.
# * STDINFILE - optional input file
# * GOLDFILE  - optional gold standard file.
#
# ------------------------------------------------------------------------------------------------ #
macro(aut_setup)

  # Setup and sanity check
  if("${APP}x" STREQUAL "x")
    message(FATAL_ERROR "You must provide a value for APP.")
  endif()

  # Set paths...
  get_filename_component(APP ${APP} ABSOLUTE)
  if(STDINFILE)
    get_filename_component(STDINFILE ${STDINFILE} ABSOLUTE)
  endif()
  get_filename_component(BINDIR ${APP} PATH)
  if(GOLDFILE)
    get_filename_component(OUTFILE ${GOLDFILE} NAME_WE)
  else()
    get_filename_component(OUTFILE ${APP} NAME_WE)
  endif()
  set(ERRFILE "${CMAKE_CURRENT_BINARY_DIR}/${OUTFILE}.err")
  set(OUTFILE "${CMAKE_CURRENT_BINARY_DIR}/${OUTFILE}.out")

  if(NOT EXISTS "${APP}")
    message(FATAL_ERROR "Cannot find ${APP}")
  else()
    message("Testing ${APP}")
  endif()

  set(numpasses 0)
  set(numfails 0)

  if(VERBOSE_DEBUG)
    message("Running with the following parameters:")
    message(
      "   APP       = ${APP}
   BINDIR    = ${BINDIR}
   PROJECT_BINARY_DIR = ${PROJECT_BINARY_DIR}
   OUTFILE   = ${OUTFILE}
   ERRFILE   = ${ERRFILE}")
    if(STDINFILE)
      message("   STDINFILE = ${STDINFILE}")
    endif()
    if(GOLDFILE)
      message("   GOLDFILE = ${GOLDFILE}")
    endif()
  endif()

  # Look for auxillary applications that are often used by ApplicationUnitTest.cmake.
  find_program(exenumdiff numdiff)
  if(NOT EXISTS ${exenumdiff})
    message(FATAL_ERROR "Numdiff not found in PATH")
  endif()
  if(VERBOSE_DEBUG)
    message("   exenumdiff = ${exenumdiff}")
  endif()

endmacro()

# ------------------------------------------------------------------------------------------------ #
# PASSMSG/FAILMSG
# ------------------------------------------------------------------------------------------------ #
macro(PASSMSG msg)
  math(EXPR numpasses "${numpasses} + 1")
  message("Test Passes: ${msg}")
endmacro()

# Increment the number of test failures and print a warning
macro(FAILMSG msg)
  math(EXPR numfails "${numfails} + 1")
  message("Test Fails: ${msg}")
endmacro()

# Increment the number of test failures
macro(ITFAILS)
  math(EXPR numfails "${numfails} + 1")
endmacro()

# ------------------------------------------------------------------------------------------------ #
# REGISTRATION
# ------------------------------------------------------------------------------------------------ #
macro(aut_register_test)
  # Register the test...

  if(EXISTS ${Draco_SOURCE_DIR}/config/ApplicationUnitTest.cmake)
    # this is draco
    set(DRACO_CONFIG_DIR ${Draco_SOURCE_DIR}/config)
  endif()

  set(num_procs 1)
  if(numPE)
    set(num_procs ${numPE})
  endif()
  string(REPLACE ";" " " RUN_CMD "${RUN_CMD}")

  if(VERBOSE_DEBUG)
    message(
      "
  add_test(
    NAME ${ctestname_base}${argname}
    COMMAND ${CMAKE_COMMAND}
    -D APP=${aut_APP}
    -D ARGVALUE=${argvalue}
    -D WORKDIR=${aut_WORKDIR}
    -D TESTNAME=${ctestname_base}${argname}
    -D DRACO_CONFIG_DIR=${DRACO_CONFIG_DIR}
    -D DRACO_INFO=$<TARGET_FILE_DIR:Exe_draco_info>/$<TARGET_FILE_NAME:Exe_draco_info>
    -D STDINFILE=${aut_STDINFILE}
    -D GOLDFILE=${aut_GOLDFILE}
    -D RUN_CMD=${RUN_CMD}
    -D numPE=${numPE}
    -D MPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}
    -D MPI_CORES_PER_CPU=${MPI_CORES_PER_CPU}
    -D SITENAME=${SITENAME}
    -D PROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
    -D PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
    ${BUILDENV}
    -P ${aut_DRIVER}
    )
  set_tests_properties( ${ctestname_base}${argname}
    PROPERTIES
      PASS_REGULAR_EXPRESSION \"${aut_PASS_REGEX}\"
      FAIL_REGULAR_EXPRESSION \"${aut_FAIL_REGEX}\"
      PROCESSORS              \"${num_procs}\"
    )")
    if(DEFINED aut_RESOURCE_LOCK)
      message("  set_tests_properties( ${ctestname_base}${argname}
      PROPERTIES RESOURCE_LOCK \"${aut_RESOURCE_LOCK}\" )")
    endif()
    if(DEFINED aut_RUN_AFTER)
      message("  set_tests_properties( ${ctestname_base}${argname}
      PROPERTIES DEPENDS \"${aut_RUN_AFTER}\" )")
    endif()
    if(DEFINED aut_LABELS)
      message("  set_tests_properties( ${ctestname_base}${argname}
      PROPERTIES LABELS \"${aut_LABELS}\" )")
    endif()
  endif(VERBOSE_DEBUG)

  # Look for python, which is used to drive application unit tests
  if(NOT Python_Interpreter_FOUND)
    # python should have been found when vendor_libraries.cmake was run.
    message(FATAL_ERROR "Draco requires python. Python not found in PATH.")
  endif()

  # Check to see if driver file is python or CMake
  set(PYTHON_TEST TRUE)
  string(FIND ${aut_DRIVER} ".py" PYTHON_TEST)
  if(${PYTHON_TEST} EQUAL -1)
    set(PYTHON_TEST FALSE)
  endif()

  # Set arguments that don't change between python and CMake driver
  set(SHARED_ARGUMENTS
      -DAPP=${aut_APP}
      -DARGVALUE=_av_${argvalue}_av_
      -DWORKDIR=${aut_WORKDIR}
      -DTESTNAME=${ctestname_base}${argname}
      -DDRACO_CONFIG_DIR=${DRACO_CONFIG_DIR}
      -DSTDINFILE=${aut_STDINFILE}
      -DGOLDFILE=${aut_GOLDFILE}
      -DRUN_CMD=${RUN_CMD}
      -DnumPE=${numPE}
      -D
      MPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}
      -DMPI_CORES_PER_CPU=${MPI_CORES_PER_CPU}
      -DSITENAME=${SITENAME}
      -DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
      -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
      ${BUILDENV})
  if(TARGET Exe_draco_info)
    list(APPEND SHARED_ARGUMENTS
         -DDRACO_INFO=$<TARGET_FILE_DIR:Exe_draco_info>/$<TARGET_FILE_NAME:Exe_draco_info>)
  endif()

  # Add python version if python driver file is specified

  if(${PYTHON_TEST})
    add_test(NAME ${ctestname_base}${argname} COMMAND "${Python_EXECUTABLE}" -B ${aut_DRIVER}
                                                      ${SHARED_ARGUMENTS})
  else()
    add_test(NAME ${ctestname_base}${argname} COMMAND ${CMAKE_COMMAND} ${SHARED_ARGUMENTS} -P
                                                      ${aut_DRIVER})
  endif()

  set_tests_properties(
    ${ctestname_base}${argname}
    PROPERTIES PASS_REGULAR_EXPRESSION "${aut_PASS_REGEX}" FAIL_REGULAR_EXPRESSION
               "${aut_FAIL_REGEX}" PROCESSORS "${num_procs}")
  if(DEFINED aut_RESOURCE_LOCK)
    set_tests_properties(${ctestname_base}${argname} PROPERTIES RESOURCE_LOCK
                                                                "${aut_RESOURCE_LOCK}")
  endif()
  if(DEFINED aut_RUN_AFTER)
    set_tests_properties(${ctestname_base}${argname} PROPERTIES DEPENDS "${aut_RUN_AFTER}")
  endif()
  if(DEFINED aut_LABELS)
    set_tests_properties(${ctestname_base}${argname} PROPERTIES LABELS "${aut_LABELS}")
  endif()

  unset(num_procs)
endmacro()

# ------------------------------------------------------------------------------------------------ #
# See documentation at the top of this file.
macro(add_app_unit_test)

  # These become variables of the form ${aut_APP}, etc.
  cmake_parse_arguments(
    aut "NONE" "APP;DRIVER;GOLDFILE;STDINFILE;WORKDIR"
    "BUILDENV;FAIL_REGEX;LABELS;PASS_REGEX;PE_LIST;RESOURCE_LOCK;RUN_AFTER;TEST_ARGS" ${ARGV})

  #
  # Check required intputs and set default vaules:
  #
  if(NOT EXISTS ${aut_DRIVER})
    message(FATAL_ERROR "Could not find the cmake driver script = ${aut_DRIVER}.")
  endif()
  if(NOT DEFINED aut_APP)
    message(FATAL_ERROR "You must provide a value for APP.")
  endif()
  if(NOT DEFINED aut_WORKDIR)
    set(aut_WORKDIR ${CMAKE_CURRENT_BINARY_DIR})
  endif()
  if(NOT DEFINED aut_PASS_REGEX)
    set(aut_PASS_REGEX "PASSED")
  endif()
  if(NOT DEFINED aut_FAIL_REGEX)
    set(aut_FAIL_REGEX "FAILED;fails")
  endif()
  if(DEFINED aut_GOLDFILE)
    if(NOT EXISTS ${aut_GOLDFILE})
      message(FATAL_ERROR "File not found, GOLDFILE=${aut_GOLDFILE}.")
    endif()
  endif()
  if(DEFINED aut_STDINFILE)
    if(NOT EXISTS ${aut_STDINFILE})
      message(FATAL_ERROR "File not found, STDINFILE=${aut_STDINFILE}.")
    endif()
  endif()

  # Load some information from the build environment:
  unset(RUN_CMD)

  if(DEFINED aut_PE_LIST AND ${DRACO_C4} MATCHES "MPI")

    # Parallel tests
    set(RUN_CMD "${MPIEXEC_EXECUTABLE} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG}")

  else()

    # Scalar tests
    if("${MPIEXEC_EXECUTABLE}" MATCHES "jsrun" OR "${MPIEXEC_EXECUTABLE}" MATCHES "srun")
      set(RUN_CMD "${MPIEXEC_EXECUTABLE} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 1")
    endif()
  endif()

  # Prove extra build environment to the cmake-scripted unit test
  unset(BUILDENV)
  if(DEFINED aut_BUILDENV)
    # expect a semicolon delimited list of parameters. FOO=myvalue1;BAR=myvalue2 translate this to
    # -DF00=myvalue1 -DBAR=myvalue2
    foreach(item ${aut_BUILDENV})
      list(APPEND BUILDENV "-D${item}")
    endforeach()
  endif()

  # Create a name for the test
  get_filename_component(drivername ${aut_DRIVER} NAME_WE)
  get_filename_component(package_name ${aut_DRIVER} PATH)
  get_filename_component(package_name ${package_name} PATH)
  get_filename_component(package_name ${package_name} NAME)
  set(ctestname_base ${package_name}_${drivername})
  # Make the test name safe for regex
  string(REGEX REPLACE "[+]" "x" ctestname_base ${ctestname_base})

  unset(argvalue)
  unset(argname)
  if(DEFINED aut_TEST_ARGS)

    # Create a suffix for the testname and generate a string that can be provided to the runTests
    # macro.
    if(${DRACO_C4} MATCHES "MPI" AND DEFINED aut_PE_LIST)
      foreach(numPE ${aut_PE_LIST})
        foreach(argvalue ${aut_TEST_ARGS})
          if(${argvalue} STREQUAL "none")
            set(argvalue "_${numPE}")
          else()
            get_filename_component(safe_argvalue "${argvalue}" NAME)
            string(REGEX REPLACE "[-]" "" safe_argvalue ${safe_argvalue})
            string(REGEX REPLACE "[ +.]" "_" safe_argvalue ${safe_argvalue})
            set(argname "_${numPE}_${safe_argvalue}")
          endif()
          # Register the test...
          if(VERBOSE_DEBUG)
            message("aut_register_test(${ctestname_base}${argname})")
          endif()
          aut_register_test()
        endforeach()
      endforeach()
    else()
      foreach(argvalue ${aut_TEST_ARGS})
        if(${argvalue} STREQUAL "none")
          set(argvalue "")
        else()
          get_filename_component(safe_argvalue "${argvalue}" NAME)
          string(REGEX REPLACE "[-]" "" safe_argvalue ${safe_argvalue})
          string(REGEX REPLACE "[ +.]" "_" safe_argvalue ${safe_argvalue})
          set(argname "_${safe_argvalue}")
        endif()
        # Register the test...
        if(VERBOSE_DEBUG)
          message("aut_register_test(${ctestname_base}${argname})")
        endif()
        aut_register_test()
      endforeach()
    endif()

  else(DEFINED aut_TEST_ARGS)

    # Register the test...
    if(${DRACO_C4} MATCHES "MPI" AND DEFINED aut_PE_LIST)
      foreach(numPE ${aut_PE_LIST})
        set(argname "_${numPE}")
        if(VERBOSE_DEBUG)
          message("aut_register_test(${ctestname_base}${argname})")
        endif()
        aut_register_test()
      endforeach()
    else()
      if(VERBOSE_DEBUG)
        message("aut_register_test(${ctestname_base})")
      endif()
      aut_register_test()
    endif()

  endif(DEFINED aut_TEST_ARGS)

  # cleanup
  unset(DRIVER)
  unset(APP)
  unset(STDINFILE)
  unset(GOLDFILE)
  unset(TEST_ARGS)
  unset(numPE)

endmacro()

# ------------------------------------------------------------------------------------------------ #
# Run the tests
# ------------------------------------------------------------------------------------------------ #
macro(aut_runTests)

  # Run the application and capture the output.
  message(
    "
=============================================
=== ${TESTNAME}
=============================================
")

  # Print version information
  separate_arguments(RUN_CMD UNIX_COMMAND "${RUN_CMD}")
  set(runcmd ${RUN_CMD}) # plain string with spaces (used in Capsaicin)

  if(numPE)
    # Use 1 proc to run draco_info
    set(draco_info_numPE 1)
  endif()
  if(EXISTS ${DRACO_INFO})
    execute_process(
      COMMAND ${RUN_CMD} ${draco_info_numPE} ${DRACO_INFO} --version
      RESULT_VARIABLE testres
      OUTPUT_VARIABLE testout
      ERROR_VARIABLE testerror)
    if(NOT ${testres} STREQUAL "0")
      failmsg("Unable to run '${runcmd} ${draco_info_numPE} ${DRACO_INFO} --version'")
    else()
      message("${testout}")
    endif()
  endif()
  unset(draco_info_numPE)

  if(numPE)
    string(REPLACE ".out" "-${numPE}.out" OUTFILE ${OUTFILE})
    string(REPLACE ".err" "-${numPE}.err" ERRFILE ${ERRFILE})
  endif()
  string(REGEX REPLACE "_av_([A-z-]*)_av_" "\\1" ARGVALUE "${ARGVALUE}")
  if(ARGVALUE)
    string(REGEX REPLACE "[-]" "" safe_argvalue ${ARGVALUE})
    string(REPLACE ".out" "-${safe_argvalue}.out" OUTFILE ${OUTFILE})
    string(REPLACE ".err" "-${safe_argvalue}.err" ERRFILE ${ERRFILE})
  endif()

  if(DEFINED RUN_CMD)
    string(REPLACE ";" " " run_cmd_string "${RUN_CMD}")
    message(">>> Running: ${run_cmd_string} ${numPE}")
    message(">>>          ${APP}")
    if(DEFINED ARGVALUE)
      message(">>>          ${ARGVALUE}")
    endif()
  else()
    message(">>> Running: ${APP} ${ARGVALUE}")
  endif()
  if(EXISTS ${STDINFILE})
    set(INPUT_FILE "INPUT_FILE ${STDINFILE}")
    message(">>>          < ${STDINFILE}")
  endif()
  message(">>>          > ${OUTFILE}
")

  # Run the application capturing all output.

  separate_arguments(INPUT_FILE UNIX_COMMAND "${INPUT_FILE}")
  separate_arguments(ARGVALUE UNIX_COMMAND "${ARGVALUE}")

  execute_process(
    COMMAND ${RUN_CMD} ${numPE} ${APP} ${ARGVALUE}
    WORKING_DIRECTORY ${WORKDIR} ${INPUT_FILE}
    RESULT_VARIABLE testres
    OUTPUT_VARIABLE testout
    ERROR_VARIABLE testerror)

  # Convert the ARGVALUE from a list back into a space separated string.
  if(ARGVALUE)
    string(REGEX REPLACE ";" " " ARGVALUE ${ARGVALUE})
  endif()

  if(FALSE)
    # This block was too slow for some Capsaicin tests

    # Preserve blank lines by settting the variable "Esc" to the ASCII value 27 - basically
    # something which is unlikely to conflict with anything in the file contents.
    string(ASCII 27 Esc)
    string(REGEX REPLACE "\n" "${Esc};" testout "${testout}")
    unset(newout)
    foreach(line ${testout})
      if(NOT line MATCHES "Application [0-9]* resources: utime")
        string(REGEX REPLACE "${Esc}" "\n" line ${line})
        set(newout "${newout}${line}")
      endif()
    endforeach()
    set(testout ${newout})
  endif()

  # now write the cleaned up file
  file(WRITE ${OUTFILE} "${testout}")

  # Ensure there are no errors
  if(NOT "${testres}" STREQUAL "0")
    message(FATAL_ERROR "Test FAILED:" "last message written to stderr: '${testerror}"
                        "See ${testout} for full details.")
  else()
    message("${testout}")
    passmsg("Application ran to completion.")
  endif()

endmacro()

# ------------------------------------------------------------------------------------------------ #
# Set numdiff run command
# ------------------------------------------------------------------------------------------------ #
function(set_numdiff_run_cmd RUN_CMD numdiff_run_cmd)
  if("${RUN_CMD}x" STREQUAL "x")
    set(numdiff_run_cmd "")
  else()
    if(NOT "${RUN_CMD}" MATCHES ";")
      separate_arguments(numdiff_run_cmd UNIX_COMMAND "${RUN_CMD}")
    else()
      set(numdiff_run_cmd ${RUN_CMD})
    endif()
    if("${MPIEXEC_EXECUTABLE}" MATCHES "mpiexec")
      # For Cray environments, let numdiff run on the login node.
      set(numdiff_run_cmd "")
    elseif(numPE)
      # Use 1 processor for srun, ssh, etc.
      set(numdiff_run_cmd "${numdiff_run_cmd};1")
    endif()
  endif()
  set(numdiff_run_cmd
      "${numdiff_run_cmd}"
      PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------------------------ #
# Run numdiff
# ------------------------------------------------------------------------------------------------ #
macro(aut_numdiff)

  set_numdiff_run_cmd("${RUN_CMD}" numdiff_run_cmd)
  string(REPLACE ";" " " pretty_run_cmd "${numdiff_run_cmd}")
  message(
    "Comparing output to goldfile:
${pretty_run_cmd} ${exenumdiff} ${ARGV2} ${ARGV3} ${ARGV4} ${ARGV5} ${ARGV6} \\
   ${OUTFILE} \\
   ${GOLDFILE}")
  execute_process(
    COMMAND ${numdiff_run_cmd} ${exenumdiff} ${ARGV2} ${ARGV3} ${ARGV4} ${ARGV5} ${ARGV6} ${OUTFILE}
            ${GOLDFILE}
    RESULT_VARIABLE numdiffres
    OUTPUT_VARIABLE numdiffout
    ERROR_VARIABLE numdifferror)
  if(${numdiffres} STREQUAL 0)
    passmsg("gold matches out.
")
  else()
    failmsg("gold does not match out.
numdiff output = ${numdiffout}")
  endif()
  unset(numdiff_run_cmd)

endmacro()

# ------------------------------------------------------------------------------------------------ #
# Run numdiff given 2 files
# ------------------------------------------------------------------------------------------------ #
macro(aut_numdiff_2files file1 file2)

  # Sometimes we are looking at a shared filesystem from different nodes.  If a file doesn't exist,
  # touch the file.  If the touch creates the file, then running numdiff should fail.
  if(NOT EXISTS ${file1})
    execute_process(COMMAND "${CMAKE_COMMAND}" -E touch ${file1})
  endif()
  if(NOT EXISTS ${file2})
    execute_process(COMMAND "${CMAKE_COMMAND}" -E touch ${file2})
  endif()

  set_numdiff_run_cmd("${RUN_CMD}" numdiff_run_cmd)
  string(REPLACE ";" " " pretty_run_cmd "${numdiff_run_cmd}")
  message(
    "
Comparing files:
${pretty_run_cmd} ${exenumdiff} ${ARGV2} ${ARGV3} ${ARGV4} ${ARGV5} ${ARGV6} \\
   ${file1} \\
   ${file2}")

  execute_process(
    COMMAND ${numdiff_run_cmd} ${exenumdiff} ${ARGV2} ${ARGV3} ${ARGV4} ${ARGV5} ${ARGV6} ${file1}
            ${file2}
    RESULT_VARIABLE numdiffres
    OUTPUT_VARIABLE numdiffout
    ERROR_VARIABLE numdifferror)
  if(${numdiffres} STREQUAL 0)
    passmsg("files match!
")
  else()
    failmsg("files do not match.
numdiff output = ${numdiffout}")
  endif()
  unset(numdiff_run_cmd)

endmacro()

# ------------------------------------------------------------------------------------------------ #
# Run gdiff (Capsaicin)
#
# Usage: aut_gdiff(input_file)
#
# GDIFF and possibly PGDIFF must be provided when registering the test:
#
# ~~~
# add_app_unit_test(
#   DRIVER   ${CMAKE_CURRENT_SOURCE_DIR}/tstAnaheim.cmake
#   APP      $<TARGET_FILE_DIR:Exe_anaheim>/$<TARGET_FILE_NAME:Exe_anaheim>
#   BUILDENV "GDIFF=$<TARGET_FILE_DIR:Exe_gdiff>/$<TARGET_FILE_NAME:Exe_gdiff>;\
#             PGDIFF=$<TARGET_FILE_DIR:Exe_pgdiff>/$<TARGET_FILE_NAME:Exe_pgdiff>"
# ~~~
# ------------------------------------------------------------------------------------------------ #
macro(aut_gdiff infile)

  # Sanity checks
  #
  # 1. Must be able to find gdiff.  Location is specified via BUILDENV
  #    "GDIFF=$<TARGET_FILE_DIR:Exe_gdiff>/$<TARGET_FILE_NAME:Exe_gdiff>"
  # 2. Input file ($1) must exist

  if(NOT EXISTS ${GDIFF})
    failmsg("Could not find gdiff!  Did you list it when registering this test?
GDIFF = ${GDIFF}")
  endif()
  if(NOT EXISTS ${infile})
    failmsg("Could not find specified intput file (${infile})!
Did you list it when registering this test?")
  endif()
  if(numPE)
    if("${numPE}" GREATER "1" AND NOT EXISTS ${PGDIFF})
      failmsg("If numPE > 1, you must provide the path to PGDIFF!")
    endif()
  endif()

  # ----------------------------------------
  # Choose pgdiff or gdiff

  if(numPE AND "${numPE}" GREATER "1")
    set(pgdiff_gdiff ${RUN_CMD} ${numPE} ${PGDIFF})
  else()
    # Use 1 proc to run gdiff
    set(pgdiff_gdiff ${RUN_CMD} 1 ${GDIFF})
  endif()

  # ----------------------------------------
  # Run GDIFF or PGDIFF

  separate_arguments(pgdiff_gdiff UNIX_COMMAND "${pgdiff_gdiff}")
  # pretty print string
  string(REPLACE ";" " " pgdiff_gdiff_string "${pgdiff_gdiff}")
  message(
    "
Comparing output to goldfile via gdiff:
   ${pgdiff_gdiff_string} ${infile}
")
  execute_process(
    COMMAND ${pgdiff_gdiff} ${infile}
    RESULT_VARIABLE gdiffres
    OUTPUT_VARIABLE gdiffout
    ERROR_VARIABLE gdifferror)
  if(${gdiffres} STREQUAL 0)
    passmsg("gdiff returned 0.")
  else()
    failmsg("gdiff returned non-zero.")
    message("gdiff messages =
${gdiffout}")
  endif()

  # should be no occurance of "FAILED"
  string(FIND "${gdiffout}" "FAILED" POS1)
  # should be  at least one occurance of "passed"
  string(FIND "${gdiffout}" "passed" POS2)

  if(${POS1} GREATER 0 OR ${POS2} EQUAL 0)
    # found failures or no passes.
    failmsg("Failed identical file check ( ${POS1}, ${POS2} ).")
  else()
    passmsg("Passed identical file check ( ${POS1}, ${POS2} ).")
  endif()

  unset(pgdiff_gdiff)

endmacro()

# ------------------------------------------------------------------------------------------------ #
# Run the tests
# ------------------------------------------------------------------------------------------------ #
macro(aut_report)

  message("
*****************************************************************")
  if(${numpasses} GREATER 0 AND ${numfails} STREQUAL 0)
    message("**** ${TESTNAME}: PASSED.")
  else()
    message("**** ${TESTNAME}: FAILED.")
  endif()
  message("*****************************************************************
")

endmacro()

# ------------------------------------------------------------------------------------------------ #
# End ApplicationUnitTest.cmake
# ------------------------------------------------------------------------------------------------ #
