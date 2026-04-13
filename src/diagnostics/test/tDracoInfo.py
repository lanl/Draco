#---*-python-*-------------------------------------------------------------------------------------#
# file   diagnostics/test/tDracoInfo.py
# author Alex Long <along@lanl.gov>
# date   Wednesday, September 14, 2016, 14:16 pm
# brief  This is a CTest script that is used to test bin/draco_info.
# note   Copyright (C) 2016-2026 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#

import sys
import re

try:
  #------------------------------------------------------------------------------------------------#
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
  #------------------------------------------------------------------------------------------------#

  # Setup test using sys.argv and run:
  tDracoInfo = UnitTest()
  tDracoInfo.aut_runTests()

  #------------------------------------------------------------------------------------------------#
  # Examine the output to determine if the test passed
  #------------------------------------------------------------------------------------------------#

  # Scan output for various strings
  found_copyright = \
    tDracoInfo.output_contains("Copyright (C)")
  found_system_type = \
    tDracoInfo.output_contains("System type")
  found_build_data = \
    tDracoInfo.output_contains("build date" )
  found_draco_team = \
    tDracoInfo.output_contains("Draco Team" )
  found_help_string = \
    tDracoInfo.output_contains("Options without arguments:")

  # There are several versions of this test
  # version command
  if tDracoInfo.check_arg_value("ARGVALUE", "--version"):
    if(found_build_data):
      tDracoInfo.passmsg( "Found build date")
    else:
      tDracoInfo.failmsg( "Did not find build date")

    if(found_copyright):
      tDracoInfo.failmsg( "Found copyright date")
    else:
      tDracoInfo.passmsg( "Did not find copyright date")

  # brief command
  elif tDracoInfo.check_arg_value("ARGVALUE", "--brief"):
      if(found_build_data):
        tDracoInfo.passmsg( "Found build date")
      else:
        tDracoInfo.failmsg( "Did not find build date")

      if(found_copyright):
        tDracoInfo.passmsg( "Found copyright date")
      else:
        tDracoInfo.failmsg( "Did not find copyright date")

      if(found_system_type):
        tDracoInfo.failmsg( "Found system type id")
      else:
        tDracoInfo.passmsg( "Did not find system type id")

  # author command
  elif tDracoInfo.check_arg_value("ARGVALUE", "--author"):
      if(found_draco_team):
        tDracoInfo.passmsg( "Found CAI-2 Draco Team")
      else:
        tDracoInfo.failmsg( "Did not find CAI-2 Draco Team")

  # author command
  elif tDracoInfo.check_arg_value("ARGVALUE", "--help") or \
      tDracoInfo.check_arg_value("ARGVALUE", "--foo"):
      if(found_help_string):
        tDracoInfo.passmsg( "Found expected help output")
      else:
        tDracoInfo.failmsg( "Did not find expected help output")

  # no arguments
  elif tDracoInfo.check_arg_value("ARGVALUE", ""):
    if(found_copyright):
      tDracoInfo.passmsg( "Found copyright date")
    else:
      tDracoInfo.failmsg( "Did not find copyright date")

    if(found_system_type):
      tDracoInfo.passmsg( "Found system type id")
    else:
      tDracoInfo.failmsg( "Did not find system type id")

  else:
    if(found_help_string):
      tDracoInfo.passmsg( "Found expected help output")
    else:
      tDracoInfo.failmsg( "Did not find expected help output")

  #------------------------------------------------------------------------------------------------#
  # Final report
  #------------------------------------------------------------------------------------------------#
  tDracoInfo.aut_report()

#--------------------------------------------------------------------------------------------------#
# Handle outstanding exceptions
#--------------------------------------------------------------------------------------------------#
except Exception:
  print("Caught exception: {0}  {1}".format( sys.exc_info()[0], \
    sys.exc_info()[1]))
  print("*****************************************************************")
  print("**** TEST FAILED.")
  print("*****************************************************************")

#--------------------------------------------------------------------------------------------------#
# End
#--------------------------------------------------------------------------------------------------#
