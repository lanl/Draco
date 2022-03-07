#-*-python-*---------------------------------------------------------------------------------------#
# file   viz/test/tstEnsight_Translator_Unstructured.py
# author Alex Long <along@lanl.gov>
# date   Wednesday, September 14, 2016, 14:16 pm
# brief  This is a Python script that is used to check the output from
#        viz/test/tstEnsight_Translator
# note   Copyright (C) 2018-2022 Triad National Security, LLC., All rights reserved.
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
    tEnsight_Translator = UnitTest()
    tEnsight_Translator.aut_runTests()
    numpe = tEnsight_Translator.numPE

    ##--------------------------------------------------------------------------##
    # ASCII files: Run numdiff on the following files:
    ##--------------------------------------------------------------------------##

    vars = ["geo", "Temperatures", "Pressure", "Velocity", "Densities"]

    for var in vars:
        # testproblem_ensight directory
        tEnsight_Translator.diff_two_files("PROJECT_BINARY_DIR",
                                           "unstr2d_testproblem_serial_{0}_ensight/{1}/data.0001".format(
                                               numpe, var), "PROJECT_SOURCE_DIR",
                                           "unstr2d_bench_serial/{0}.0001".format(var))
        tEnsight_Translator.diff_two_files("PROJECT_BINARY_DIR",
                                           "unstr2d_testproblem_parallel_{0}_ensight/{1}/data.0001".format(
                                               numpe, var), "PROJECT_SOURCE_DIR",
                                           "unstr2d_bench_parallel_{0}/{1}.0001".format(numpe, var))
        tEnsight_Translator.diff_two_files("PROJECT_BINARY_DIR",
                                           "unstr3d_testproblem_serial_{0}_ensight/{1}/data.0001".format(
                                               numpe, var), "PROJECT_SOURCE_DIR",
                                           "unstr3d_bench_serial/{0}.0001".format(var))
        tEnsight_Translator.diff_two_files("PROJECT_BINARY_DIR",
                                           "unstr3d_testproblem_parallel_{0}_ensight/{1}/data.0001".format(
                                               numpe, var), "PROJECT_SOURCE_DIR",
                                           "unstr3d_bench_parallel_{0}/{1}.0001".format(numpe, var))

    # diff binary files if on little endian system
    if tEnsight_Translator.is_little_endian():
        if any(platform.win32_ver()):
            diff_prog = "fc"
        else:
            diff_prog = "diff"
        print("\nSystem is little endian, diffing binary files\n")
        for var in vars:
            # testproblem_binary_ensight directory
            tEnsight_Translator.diff_two_files("PROJECT_BINARY_DIR",
                                               "unstr2d_testproblem_serial_{0}_binary_ensight/{1}/data.0001".format(
                                                   numpe, var),
                                               "PROJECT_SOURCE_DIR", "unstr2d_bench_serial/{0}.bin.0001".format(
                                                   var),
                                               diff_name=diff_prog)
            tEnsight_Translator.diff_two_files("PROJECT_BINARY_DIR",
                                               "unstr2d_testproblem_parallel_{0}_binary_ensight/{1}/data.0001".format(
                                                   numpe, var),
                                               "PROJECT_SOURCE_DIR", "unstr2d_bench_parallel_{0}/{1}.bin.0001".format(
                                                   numpe, var),
                                               diff_name=diff_prog)
            tEnsight_Translator.diff_two_files("PROJECT_BINARY_DIR",
                                               "unstr3d_testproblem_serial_{0}_binary_ensight/{1}/data.0001".format(
                                                   numpe, var),
                                               "PROJECT_SOURCE_DIR", "unstr3d_bench_serial/{0}.bin.0001".format(
                                                   var),
                                               diff_name=diff_prog)
            tEnsight_Translator.diff_two_files("PROJECT_BINARY_DIR",
                                               "unstr3d_testproblem_parallel_{0}_binary_ensight/{1}/data.0001".format(
                                                   numpe, var),
                                               "PROJECT_SOURCE_DIR", "unstr3d_bench_parallel_{0}/{1}.bin.0001".format(
                                                   numpe, var),
                                               diff_name=diff_prog)

    ##--------------------------------------------------------------------------##
    # Final report
    ##--------------------------------------------------------------------------##
    tEnsight_Translator.aut_report()

##----------------------------------------------------------------------------##
# Handle outstanding exceptions
##----------------------------------------------------------------------------##
except Exception:
    print("Caught exception: {0}  {1}".format(sys.exc_info()[0],
                                              sys.exc_info()[1]))
    print("*****************************************************************")
    print("**** TEST FAILED.")
    print("*****************************************************************")

##---------------------------------------------------------------------------##
# End
##---------------------------------------------------------------------------##
