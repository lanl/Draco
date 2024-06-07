#-----------------------------*-python-*---------------------------------------#
# file   dsxx/test/test_fpe_trap.py
# author Kelly Thompson <kgt@lanl.gov>
# date   Monday, Nov 28, 2016, 16:40 pm
# brief  This is a Python script that is used to test the fpe_trap features in
#        dsxx.
# note   Copyright (C) 2024 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
import sys
import re
import platform

try:

  #----------------------------------------------------------------------------#
  # Boilerplate code to get the location of the application_unit_test.py in
  # draco/config
  re_draco_config_dir = re.compile("DRACO_CONFIG_DIR=([^\s]*)")
  draco_config_dir = ""
  for arg in sys.argv:
    if (re_draco_config_dir.search(arg)):
      draco_config_dir = re_draco_config_dir.findall(arg)[0]
  if not draco_config_dir:
    raise Exception("Draco config directory not found, exiting")

  # import unit test functions
  sys.path.append(draco_config_dir)
  from application_unit_test import UnitTest
  #----------------------------------------------------------------------------#

  # Setup test using sys.argv and run:
  tFpeTrap = UnitTest()
  tFpeTrap.aut_runTests(True)

  ##---------------------------------------------------------------------------##
  ## Check the output
  ##---------------------------------------------------------------------------##

  print("\nChecking the generated output files...\n")
  case_number = tFpeTrap.get_arg_value( 'ARGVALUE' )

  #--------------------------------------------------------------------------------------------------#
  # Case 0: attempts to calculate 1.0 + zero + sqrt(-1*-1)
  #         should not throw IEEE exception.
  if case_number == '0':

      string_found = tFpeTrap.output_contains("result = 2")
      if( string_found ):
          tFpeTrap.passmsg("Found expected result (2).")
      else:
          tFpeTrap.failmsg("Failed to find expected result (2).")


  #--------------------------------------------------------------------------------------------------#
  # Case 1: attempts to divide by zero
  elif case_number == '1':

    # Check the test type
    string_found = tFpeTrap.output_contains("trying a div_by_zero operation")
    if( string_found ):
      tFpeTrap.passmsg("Case 1 -> division by zero test")
    else:
      tFpeTrap.failmsg("Case 1 did not try the division by zero test")

    if any(platform.win32_ver()):
      # Signaling error: A SIGFPE was detected!
      string_found = tFpeTrap.output_contains("A SIGFPE was detected!")
    if( not string_found ):
      # Signaling error: SIGFPE (Floating point divide by zero)
      string_found = tFpeTrap.output_contains("Floating point divide by zero")
    if( not string_found ):
      # Signaling error: runtime error: division by zero
      string_found = tFpeTrap.output_contains("runtime error: division by zero")

    if( string_found ):
      tFpeTrap.passmsg("Caught SIGFPE (Floating point divide by zero)")
    else:
      tFpeTrap.failmsg("Failed to catch SIGFPE (Floating point divide by zero)")

      print("Standard out:")
      with open(tFpeTrap.outfile) as f:
        for line in f:
          print("%s" % line)

      print("Standard error:")
      with open(tFpeTrap.errfile) as f:
        for line in f:
          print("%s" % line)

  #--------------------------------------------------------------------------------------------------#
  # Case 2: attempts to evaluate sqrt(-1.0)
  elif case_number == '2':

    # Check the test type
    string_found = tFpeTrap.output_contains("trying to evaluate sqrt")
    if( string_found ):
      tFpeTrap.passmsg("Case 2 -> sqrt(-1.0) test case")
    else:
      tFpeTrap.failmsg("Case 2 did not try the sqrt(-1.0) case")

    # As of 2016-11-29:
    # - GCC throws the FE_INVALID IEEE signal (stderr file)
    # - Intel 17 throws a C++ exception (stdout file)
    #
    # Look for FE_INVALID first, if that isn't found, allow the test to pass if
    # the C++ exception is detected.

    # Signaling error: SIGFPE (Invalid floating point operation)
    string_found = tFpeTrap.error_contains("Invalid floating point operation")
    if( string_found ):
      tFpeTrap.passmsg("Caught SIGFPE (Invalid floating point operation)")
    else:
      # 2nd chance: also look in the stdout file.
      string_found = tFpeTrap.output_contains("Invalid floating point operation")
      if( string_found ):
        tFpeTrap.passmsg("Caught SIGFPE (Invalid floating point operation)")
      else:
        tFpeTrap.failmsg("Failed to catch SIGFPE (Invalid floating point operation)")

        print("Standard out:")
        with open(tFpeTrap.outfile) as f:
          for line in f:
            print("%s" % line)

        print("Standard error:")
        with open(tFpeTrap.errfile) as f:
          for line in f:
            print("%s" % line)

  #--------------------------------------------------------------------------------------------------#
  # Case 3: An overflow condition is generated.
  elif case_number == '3':

    # Check the test type
    string_found = tFpeTrap.output_contains("trying to cause an overflow condition")
    if( string_found ):
      tFpeTrap.passmsg("Case 3 -> overflow condition test")
    else:
      tFpeTrap.failmsg("Case 3 did not try the overflow condition test")

    if any(platform.win32_ver()):
      # Signaling error: A SIGFPE was detected!
      string_found = tFpeTrap.output_contains("A SIGFPE was detected!")
    else:
      # Signaling error: SIGFPE (Floating point overflow)
      string_found = tFpeTrap.error_contains("Floating point overflow")
      # 2nd try - look in the stdout stream
      if not string_found:
        string_found = tFpeTrap.output_contains("Floating point overflow")

    if( string_found ):
      tFpeTrap.passmsg("Caught SIGFPE (Floating point overflow)")
    else:
      tFpeTrap.failmsg("Failed to catch SIGFPE (Floating point overflow)")

      print("Standard out:")
      with open(tFpeTrap.outfile) as f:
        for line in f:
          print("%s" % line)

      print("Standard error:")
      with open(tFpeTrap.errfile) as f:
        for line in f:
          print("%s" % line)

  print(" ")

  ##---------------------------------------------------------------------------##
  ## Final report
  ##---------------------------------------------------------------------------##
  tFpeTrap.aut_report()

##----------------------------------------------------------------------------##
## Handle outstanding exceptions
##----------------------------------------------------------------------------##
except Exception:
  print("Caught exception: {0}  {1}".format( sys.exc_info()[0], \
    sys.exc_info()[1]))
  print("*****************************************************************")
  print("**** TEST FAILED.")
  print("*****************************************************************")

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
