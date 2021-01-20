#-----------------------------*-python-*---------------------------------------#
# file   c4/test/tstYthi.py
# author Kelly Thompson <kgt@lanl.gov>
# date   Tuesday, Jun 05, 2018, 17:10 pm
# brief  This is a Python script that is used to test c4/bin/ythi
# note   Copyright (C) 2018, Triad National Security, LLC.
#        All rights reserved.
#--------------------------------------------------------------------------------------------------#
import sys
import re

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

  ##---------------------------------------------------------------------------##
  ## Test Ythi
  ##---------------------------------------------------------------------------##

  # Setup test using sys.argv and run:
  tstYthi = UnitTest()
  tstYthi.aut_runTests()

  ##---------------------------------------------------------------------------##
  ## Check the output
  ##---------------------------------------------------------------------------##

  print("\nChecking the generated output file...\n")

  # These strings should be found:

  string_found = \
    tstYthi.output_contains("Thread 000, core affinity = ")
  if(string_found):
    tstYthi.passmsg("Found thread 0.report")
  else:
    tstYthi.failmsg("Did not find thread 0 report.")

  string_found = \
    tstYthi.output_contains("Rank 00000, Thread")
  if(string_found):
    tstYthi.passmsg("Found MPI rank 0.")
  else:
    tstYthi.failmsg("Did not find MPI rank 0.")

  # Diff the output vs a gold file.
  # tstYthi.aut_numdiff()

  ##---------------------------------------------------------------------------##
  ## Final report
  ##---------------------------------------------------------------------------##
  tstYthi.aut_report()

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
