//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   device/test/gpu_device_info.cc
 * \author Alex Long
 * \date   Thu Mar 21 15:28:48 2019
 * \brief  Simple test of the CUDA Runtime API through the GPU_Device object
 * \note   Copyright (C) 2019-2020 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Assert.hh"
#include "dsxx/DracoStrings.hh"
#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/SystemCall.hh"
#include "dsxx/path.hh"
#include "device/GPU_Device.hh"
#include "device/config.h"
#include <iostream>
#include <sstream>
#include <time.h>

//------------------------------------------------------------------------------------------------//
// query_device
//------------------------------------------------------------------------------------------------//

void query_device(rtt_dsxx::ScalarUnitTest &ut) {
  using namespace std;

  cout << "Starting gpu_hello_driver_api::query_device()...\n" << endl;

  // Create a GPU_Device object.
  // Initialize the CUDA library and sets device and context handles.
  rtt_device::GPU_Device gpu;

  // Create and then print a summary of the devices found.
  std::ostringstream out;
  size_t const numDev(gpu.numDevicesAvailable());
  out << "GPU device summary:\n\n"
      << "   Number of devices found: " << numDev << "\n"
      << endl;
  for (size_t device = 0; device < numDev; ++device)
    gpu.printDeviceSummary(device, out);

  // Print the message to stdout
  cout << out.str();

  // Parse the output
  bool verbose(false);
  std::map<std::string, unsigned> wordCount = rtt_dsxx::get_word_count(out, verbose);

  FAIL_IF_NOT(wordCount[string("Device")] == numDev);
  // successful test output
  if (ut.numFails == 0)
    PASSMSG("gpu_device_info_test query_device test OK.");
  return;
}

//------------------------------------------------------------------------------------------------//
// Main
//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  using namespace std;

  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    query_device(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of gpu_device_info.cc
//------------------------------------------------------------------------------------------------//
