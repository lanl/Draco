//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   predict/test/tstnnetwork_interface.cc
 * \author Mathew Cleveland
 * \date   Nov. 10th 2020
 * \brief  KDE function tests
 * \note   Copyright (C) 2021-2021 Triad National Security, LLC., All rights reserved.
 *          */
//------------------------------------------------------------------------------------------------//

#include "../nnetwork_interface.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include <numeric>
#include <string>

using namespace rtt_dsxx;
using namespace rtt_c4;
using namespace rtt_predict;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
//
//
void test_replication(ParallelUnitTest &ut) {

#ifdef LIBTORCH
  std::string nn_file_name = "kde.pt";
  const int input_dim{100};
  const int output_dim{1};
  nnetwork_interface net(nn_file_name);
  if (!net.valid())
    ITFAILS;
  std::vector<float> input(input_dim, 1.0f / static_cast<float>(input_dim));
  std::vector<float> result = net.predict(input, input_dim, output_dim);
  if (!rtt_dsxx::soft_equiv(double(result[0]), 0.4891, 1.0e-4))
    ITFAILS;
  if (ut.numFails == 0) {
    PASSMSG("KDE checks pass");
  } else {
    FAILMSG("KDE checks failed");
  }
#else
  try {
    nnetwork_interface net("junk");
    FAILMSG("Failed to catch invalid NN backend error");
  } catch (...) {
    PASSMSG("Expected assertion caught: Not compiled with valid neural network backed");
  }
  if (ut.numFails == 0) {
    PASSMSG("KDE checks pass");
  } else {
    FAILMSG("KDE checks failed");
  }
#endif
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ParallelUnitTest ut(argc, argv, release);
  try {
    // >>> UNIT TESTS
    test_replication(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstnnetwork_interface.cc
//------------------------------------------------------------------------------------------------//
