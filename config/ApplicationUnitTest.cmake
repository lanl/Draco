#-----------------------------*-cmake-*----------------------------------------#
# file   config/ApplicationUnitTest.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   Monday, Nov 19, 2012, 16:21 pm
# brief  Provide macros that aid in creating unit tests that run
#        interactive user codes (i.e.: run a binary that reads an
#        input file and diff the resulting output file).
# note   Copyright (C) 2012-2015, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id: CMakeLists.txt 6732 2012-09-05 22:28:18Z kellyt $
#------------------------------------------------------------------------------#

# Reference: https://rtt.lanl.gov/redmine/projects/draco/wiki/CMake-based_ApplicationUnitTest
# Example: draco/src/diagnostics/test/tDracoInfo.cmake (and associated CMakeLists.txt).

##---------------------------------------------------------------------------------------##
## Example from draco/src/diagnostics/test/tDracoInfo.cmake
##---------------------------------------------------------------------------------------##

# Use config/ApplicationUnitTest.cmake test registration:
#
# include( ApplicationUnitTest )
# add_app_unit_test(
#   DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/tDracoInfo.cmake
#   APP    $<TARGET_FILE_DIR:Exe_draco_info>/$<TARGET_FILE_NAME:Exe_draco_info>
#   LABELS nomemcheck )

# The above will generate a test with data similar to this:
#
# add_test(
#    NAME diagnostics_tDracoInfo
#    COMMAND /yellow/usr/projects/draco/vendors/cmake-3.2.2-Linux-x86_64/bin/cmake
#      -D APP              = $<TARGET_FILE_DIR:Exe_draco_info>/$<TARGET_FILE_NAME:Exe_draco_info>
#      -D WORKDIR          = /users/kellyt/build/ml/intel-mpid/d/src/diagnostics/test
#      -D TESTNAME         = diagnostics_tDracoInfo
#      -D DRACO_CONFIG_DIR = /users/kellyt/draco/config
#      -D DRACO_INFO       = /users/kellyt/build/ml/intel-mpid/d/src/diagnostics/draco_info
#      -D RUN_CMD          =
#      -P /users/kellyt/draco/src/diagnostics/test/tDracoInfo.cmake
#    )
# set_tests_properties( diagnostics_draco_info
#    PROPERTIES
#      PASS_REGULAR_EXPRESSION Passes
#      FAIL_REGULAR_EXPRESSION Fails
#      LABELS nomemcheck
#    )

# Variables defined above can be used in this script.

##---------------------------------------------------------------------------------------##

include( parse_arguments )

# set( VERBOSE_DEBUG OFF )

##---------------------------------------------------------------------------##
## Check values for $APP
##
## Requires:
##  APP        - name of executable that will be run
##
## Returns:
##   APP       - cleaned up path to executable
##   STDINFILE - cleaned up path to intput file.
##   BINDIR    - Directory location of binary file
##   PROJECT_BINARY_DIR - Parent directory of BINDIR
##   GOLDFILE  - cleaned up path to gold standard file.
##   OUTFILE   - Output filename derived from GOLDFILE.
##   ERRFILE   - Output(error) filename derived from GOLDFILE.
##
## if VERBOSE is set, also echo input values.
##   APP       - path name for executable.
##   STDINFILE - optional input file
##   GOLDFILE  - optional gold standard file.
##---------------------------------------------------------------------------##
macro( aut_setup )

  # Setup and sanity check
  if( "${APP}x" STREQUAL "x" )
    message( FATAL_ERROR "You must provide a value for APP." )
  endif()

  # Set paths...
  get_filename_component( APP ${APP} ABSOLUTE )
  if( STDINFILE )
    get_filename_component( STDINFILE ${STDINFILE} ABSOLUTE )
  endif()
  get_filename_component( BINDIR ${APP} PATH )
  get_filename_component( PROJECT_BINARY_DIR ${BINDIR} PATH )
  if( GOLDFILE )
    get_filename_component( OUTFILE ${GOLDFILE} NAME_WE )
  else()
    get_filename_component( OUTFILE ${APP} NAME_WE )
  endif()
  set( ERRFILE "${CMAKE_CURRENT_BINARY_DIR}/${OUTFILE}.err")
  set( OUTFILE "${CMAKE_CURRENT_BINARY_DIR}/${OUTFILE}.out")

  if( NOT EXISTS "${APP}" )
    message( FATAL_ERROR "Cannot find ${APP}")
  else()
    message("Testing ${APP}")
  endif()

  set( numpasses 0 )
  set( numfails  0 )

  if( VERBOSE_DEBUG )
    message("Running with the following parameters:")
    message("   APP       = ${APP}
   BINDIR    = ${BINDIR}
   PROJECT_BINARY_DIR = ${PROJECT_BINARY_DIR}
   OUTFILE   = ${OUTFILE}
   ERRFILE   = ${ERRFILE}")
    if( STDINFILE )
      message("   STDINFILE = ${STDINFILE}")
    endif()
    if( GOLDFILE )
      message("   GOLDFILE = ${GOLDFILE}" )
    endif()
# This will cause REGEX PASS/FAIL to fail!
#  message("   numpasses   = ${numpasses}
#   numfails  = ${numfails}")
  endif()

endmacro()

##---------------------------------------------------------------------------##
## PASSMSG/FAILMSG
##---------------------------------------------------------------------------##

macro(PASSMSG msg)
    math( EXPR numpasses "${numpasses} + 1" )
    message( "Test Passes: ${msg}")
endmacro()

macro(FAILMSG msg)
    math( EXPR numfails "${numfails} + 1" )
    message( "Test Fails: ${msg}")
endmacro()

macro(ITFAILS)
    math( EXPR numfails "${numfails} + 1" )
endmacro()

##---------------------------------------------------------------------------##
## REGISTRATION
##---------------------------------------------------------------------------##

macro( aut_register_test )
  # Register the test...

  if( EXISTS ${Draco_SOURCE_DIR}/config/ApplicationUnitTest.cmake )
    # this is draco
    set( DRACO_CONFIG_DIR ${Draco_SOURCE_DIR}/config )
  endif()

  if( "${numPE}x" STREQUAL "x" )
    set(num_procs 1)
  else()
    set(num_procs ${numPE} )
  endif()

  if( VERBOSE_DEBUG )
    message("
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
    ${BUILDENV}
    -P ${aut_DRIVER}
    )
  set_tests_properties( ${ctestname_base}${argname}
    PROPERTIES
    PASS_REGULAR_EXPRESSION \"${aut_PASS_REGEX}\"
    FAIL_REGULAR_EXPRESSION \"${aut_FAIL_REGEX}\"
    PROCESSORS              \"${num_procs}\"
    ${LABELS}
    )")
    if( DEFINED aut_RESOURCE_LOCK )
message("  set_tests_properties( ${ctestname_base}${argname}
      PROPERTIES RESOURCE_LOCK \"${aut_RESOURCE_LOCK}\" )" )
    endif()
  endif()

  add_test(
    NAME ${ctestname_base}${argname}
    COMMAND ${CMAKE_COMMAND}
    -DAPP=${aut_APP}
    -DARGVALUE=${argvalue}
    -DWORKDIR=${aut_WORKDIR}
    -DTESTNAME=${ctestname_base}${argname}
    -DDRACO_CONFIG_DIR=${DRACO_CONFIG_DIR}
    -DDRACO_INFO=$<TARGET_FILE_DIR:Exe_draco_info>/$<TARGET_FILE_NAME:Exe_draco_info>
    -DSTDINFILE=${aut_STDINFILE}
    -DGOLDFILE=${aut_GOLDFILE}
    -DRUN_CMD=${RUN_CMD}
    -DnumPE=${numPE}
    ${BUILDENV}
    -P ${aut_DRIVER}
    )
  set_tests_properties( ${ctestname_base}${argname}
    PROPERTIES
    PASS_REGULAR_EXPRESSION "${aut_PASS_REGEX}"
    FAIL_REGULAR_EXPRESSION "${aut_FAIL_REGEX}"
    PROCESSORS              "${num_procs}"
    ${LABELS}
    )
  if( DEFINED aut_RESOURCE_LOCK )
    set_tests_properties( ${ctestname_base}${argname}
      PROPERTIES
      RESOURCE_LOCK "${aut_RESOURCE_LOCK}" )
  endif()

  unset(num_procs)
endmacro()

##---------------------------------------------------------------------------------------##

macro( add_app_unit_test )

  # These become variables of the form ${addscalartests_SOURCES}, etc.
  parse_arguments(
    # prefix
    aut
    # list names
    "APP;BUILDENV;DRIVER;FAIL_REGEX;GOLDFILE;LABELS;PASS_REGEX;PE_LIST;RESOURCE_LOCK;STDINFILE;TEST_ARGS;WORKDIR"
    # RUN_AFTER"
    # option names
    "NONE"
    ${ARGV}
    )

  #
  # Check required intputs and set default vaules:
  #
  if( NOT EXISTS ${aut_DRIVER} )
    message( FATAL_ERROR "Could not find the cmake driver script = ${aut_DRIVER}." )
  endif()
  if( NOT DEFINED aut_APP )
    message( FATAL_ERROR "You must provide a value for APP." )
  endif()
  if( NOT DEFINED aut_WORKDIR )
    set( aut_WORKDIR ${CMAKE_CURRENT_BINARY_DIR} )
  endif()
  if( NOT DEFINED aut_PASS_REGEX )
    set( aut_PASS_REGEX "PASSED" )
  endif()
  if( NOT DEFINED aut_FAIL_REGEX )
    set( aut_FAIL_REGEX "FAILED;fails" )
  endif()
  if( DEFINED aut_LABELS )
    set( LABEL "LABELS ${aut_LABELS}" )
  endif()
  if( DEFINED aut_GOLDFILE )
    if( NOT EXISTS ${aut_GOLDFILE} )
      message( FATAL_ERROR "File not found, GOLDFILE=${aut_GOLDFILE}.")
    endif()
  endif()
  if( DEFINED aut_STDINFILE )
    if( NOT EXISTS ${aut_STDINFILE} )
      message( FATAL_ERROR "File not found, STDINFILE=${aut_STDINFILE}.")
    endif()
  endif()

  # Load some information from the build environment:
  unset( RUN_CMD )
  set( MIC_RUN_CMD  "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null $ENV{HOSTNAME}-mic0 ${Draco_BINARY_DIR}/config/run_test_on_mic.sh ${aut_WORKDIR}" )
  if( DEFINED aut_PE_LIST AND ${DRACO_C4} MATCHES "MPI" )

    # Parallel tests
    if( "${C4_MPICMD}" MATCHES "aprun" )
      set( RUN_CMD "aprun -n" )
    elseif( HAVE_MIC )
      set( RUN_CMD "${MIC_RUN_CMD} ${MPIEXEC} ${MPIEXEC_POSTFLAGS} ${MPIEXEC_NUMPROC_FLAG}" )
    else()
      set( RUN_CMD "${MPIEXEC} ${MPIEXEC_POSTFLAGS} ${MPIEXEC_NUMPROC_FLAG}")
    endif()

  else()

    # Scalar tests
    if( "${C4_MPICMD}" MATCHES "aprun" )
      set( RUN_CMD "aprun -n 1" )
    elseif( HAVE_MIC )
      set( RUN_CMD "${MIC_RUN_CMD}" )
    endif()
  endif()

  # Prove extra build environment to the cmake-scripted unit test
  unset( BUILDENV )
  if( DEFINED aut_BUILDENV )
    # expect a semicolon delimited list of parameters.
    # FOO=myvalue1;BAR=myvalue2
    # translate this to -DF00=myvalue1 -DBAR=myvalue2
    foreach(item ${aut_BUILDENV} )
      list(APPEND BUILDENV "-D${item}" )
    endforeach()
  endif()

  # Create a name for the test
  get_filename_component( drivername   ${aut_DRIVER} NAME_WE )
  get_filename_component( package_name ${aut_DRIVER} PATH )
  get_filename_component( package_name ${package_name} PATH )
  get_filename_component( package_name ${package_name} NAME )
  set( ctestname_base ${package_name}_${drivername} )

  unset( argvalue )
  unset( argname  )
  if( DEFINED aut_TEST_ARGS )

    # Create a suffix for the testname and generate a string that can
    # be provided to the runTests macro.
    if( ${DRACO_C4} MATCHES "MPI" AND DEFINED aut_PE_LIST )
      foreach( numPE ${aut_PE_LIST} )
        set( iarg "0" )
        foreach( argvalue ${aut_TEST_ARGS} )
          math( EXPR iarg "${iarg} + 1" )
          if( ${argvalue} STREQUAL "none" )
            set( argvalue "_${numPE}" )
          else()
            # set( argname "_${numPE}_arg${iarg}" )
            string( REGEX REPLACE "[-]" "" safe_argvalue ${argvalue} )
            string( REGEX REPLACE "[.]" "_" safe_argvalue ${safe_argvalue} )
            set( argname "_${numPE}_${safe_argvalue}" )
          endif()
          # Register the test...
          if( VERBOSE_DEBUG )
            message("aut_register_test(${ctestname_base}${argname})")
          endif()
          aut_register_test()
        endforeach()
      endforeach()
    else()
      set( iarg "0" )
      foreach( argvalue ${aut_TEST_ARGS} )
        math( EXPR iarg "${iarg} + 1" )
        if( ${argvalue} STREQUAL "none" )
          set( argvalue "" )
        else()
          #set( argname "_arg${iarg}" )
          string( REGEX REPLACE "[-]" "" safe_argvalue ${argvalue} )
          string( REGEX REPLACE "[.]" "_" safe_argvalue ${safe_argvalue} )
          set( argname "_${safe_argvalue}" )
        endif()
        # Register the test...
        if( VERBOSE_DEBUG )
          message("aut_register_test(${ctestname_base}${argname})")
        endif()
        aut_register_test()
      endforeach()
    endif()

  else()

    # Register the test...
    if( ${DRACO_C4} MATCHES "MPI" AND DEFINED aut_PE_LIST )
      foreach( numPE ${aut_PE_LIST} )
        set( argname "_${numPE}" )
        if( VERBOSE_DEBUG )
          message("aut_register_test(${ctestname_base}${argname})")
        endif()
        aut_register_test()
      endforeach()
    endif()

  endif()

  # cleanup
  unset( DRIVER )
  unset( APP )
  unset( STDINFILE )
  unset( GOLDFILE  )
  unset( TEST_ARGS )
  unset( numPE )
  unset( LABEL )

endmacro()

##---------------------------------------------------------------------------##
## Run the tests
##---------------------------------------------------------------------------##
macro( aut_runTests )

  # Run the application and capture the output.
  message("
=============================================
=== ${TESTNAME}
=============================================
")
# === CMake driven ApplicationUnitTest: ${TESTNAME}

  # Print version information
  set( runcmd ${RUN_CMD} ) # plain string with spaces.
  separate_arguments(RUN_CMD)
  if( numPE )
    # Use 1 proc to run draco_info
    set( draco_info_numPE 1 )
  endif()
  if( EXISTS ${DRACO_INFO} )
    execute_process(
      COMMAND ${RUN_CMD} ${draco_info_numPE} ${DRACO_INFO} --version
      RESULT_VARIABLE testres
      OUTPUT_VARIABLE testout
      ERROR_VARIABLE  testerror
      )
    if( NOT ${testres} STREQUAL "0" )
      FAILMSG("Unable to run '${RUN_CMD} ${DRACO_INFO} --version'")
    else()
      message("${testout}")
    endif()
  endif()
  unset( draco_info_numPE )

  if( numPE )
    string( REPLACE ".out" "-${numPE}.out" OUTFILE ${OUTFILE} )
    string( REPLACE ".err" "-${numPE}.err" ERRFILE ${ERRFILE} )
  endif()
  if( ARGVALUE )
    string( REGEX REPLACE "[-]" "" safe_argvalue ${ARGVALUE} )
    string( REPLACE ".out" "-${safe_argvalue}.out" OUTFILE ${OUTFILE} )
    string( REPLACE ".err" "-${safe_argvalue}.err" ERRFILE ${ERRFILE} )
  endif()

  if( DEFINED RUN_CMD )
    string( REPLACE ";" " " run_cmd_string "${RUN_CMD}" )
    message(">>> Running: ${run_cmd_string} ${numPE} ${APP}
             ${ARGVALUE}")
  else()
    message(">>> Running: ${APP} ${ARGVALUE}" )
  endif()
  if( EXISTS ${STDINFILE} )
    set( INPUT_FILE "INPUT_FILE ${STDINFILE}")
    message(">>>          < ${STDINFILE}")
  endif()
  message(">>>          > ${OUTFILE}
")

  # Run the application capturing all output.
  separate_arguments(INPUT_FILE)
  execute_process(
    COMMAND ${RUN_CMD} ${numPE} ${APP} ${ARGVALUE}
    WORKING_DIRECTORY ${WORKDIR}
    ${INPUT_FILE}
    RESULT_VARIABLE testres
    OUTPUT_VARIABLE testout
    ERROR_VARIABLE  testerror
    )

  # Capture all the output to log files:
  file( WRITE ${OUTFILE} ${testout} )
  # [2015-07-28 KT] not sure we need to dump stderr values right now
  # file( WRITE ${ERRFILE} ${testerr} )

  # Ensure there are no errors
  if( NOT "${testres}" STREQUAL "0" )
    message( FATAL_ERROR "Test FAILED:
     error message = ${testerror}")
  else()
    message("${testout}")
  endif()

endmacro()


##---------------------------------------------------------------------------##
## Run numdiff
##---------------------------------------------------------------------------##
macro( aut_numdiff )

  if( DEFINED GOLDFILE )
    ## Use numdiff to compare output
    find_program( exenumdiff numdiff )
    if( NOT EXISTS ${exenumdiff} )
      message( FATAL_ERROR "Numdiff not found in PATH")
    endif()
    if( VERBOSE_DEBUG )
      message("   exenumdiff = ${exenumdiff}" )
    endif()
  endif()

  message("Comparing output to goldfile:
${exenumdiff} \\
   ${OUTFILE} \\
   ${GOLDFILE}")
  execute_process(
    COMMAND ${exenumdiff} ${OUTFILE} ${GOLDFILE}
    RESULT_VARIABLE numdiffres
    OUTPUT_VARIABLE numdiffout
    ERROR_VARIABLE numdifferror
    )
  if( ${numdiffres} STREQUAL 0 )
    PASSMSG("gold matches out.
")
  else()
    FAILMSG("gold does not match out.
numdiff output = ${numdiffout}" )
  endif()

endmacro()

##---------------------------------------------------------------------------##
## Run gdiff (Capsaicin)
##
## Usage:
##    aut_gdiff(input_file)
##---------------------------------------------------------------------------##
macro( aut_gdiff infile )

  # Sanity checks
  # 1. Must be able to find gdiff.  Location is specified via
  #    BUILDENV  "GDIFF=$<TARGET_FILE_DIR:Exe_gdiff>/$<TARGET_FILE_NAME:Exe_gdiff>"
  # 2. Input file ($1) must exist

  if( NOT EXISTS ${GDIFF} )
    FAILMSG( "Could not find gdiff!  Did you list it when registering this test?" )
  endif()
  if( NOT EXISTS ${infile} )
    FAILMSG( "Could not find gdiff!  Did you list it when registering this test?" )
  endif()
  if( numPE )
    if( "${numPE}" GREATER "1" )
      FAILMSG( "You can only run gdiff for scalar or ! PE tests!")
    endif()
  endif()

  #----------------------------------------
  message("
Comparing output to goldfile via gdiff:
   ${GDIFF} guajillo.ndi.2D.gdiff
")
  execute_process(
    COMMAND ${GDIFF} ${infile}
    RESULT_VARIABLE gdiffres
    OUTPUT_VARIABLE gdiffout
    ERROR_VARIABLE  gdifferror
    )
  if( ${gdiffres} STREQUAL 0 )
    PASSMSG("gdiff returned 0.")
  else()
    FAILMSG("gdiff returned non-zero.")
  endif()

  # should be no occurance of "FAILED"
  string(FIND "${gdiffout}" "FAILED" POS1)
  # should be  at least one occurance of "passed"
  string(FIND "${gdiffout}" "passed" POS2)

  if( ${POS1} GREATER 0 OR ${POS2} EQUAL 0 )
    # found failures or no passes.
    FAILMSG( "Failed identical file check ( ${POS1}, ${POS2} )." )
  else()
    PASSMSG( "Passed identical file check ( ${POS1}, ${POS2} )." )
  endif()

endmacro()

##---------------------------------------------------------------------------##
## Run the tests
##---------------------------------------------------------------------------##
macro( aut_report )

  message("
*********************************************")
  if( ${numpasses} GREATER 0 AND ${numfails} STREQUAL 0 )
    message("**** ${TESTNAME}: PASSED.")
  else()
    message("**** ${TESTNAME}: FAILED.")
  endif()
  message("*********************************************
")
# message("numpasses = ${numpasses}, numfails = ${numfails}")

endmacro()

##---------------------------------------------------------------------------##
## end
##---------------------------------------------------------------------------##
