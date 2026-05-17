//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   predict/test/tstnnetwork_interface.cc
 * \author Mathew Cleveland
 * \brief  KDE function tests
 * \note   Copyright (C) 2022-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "predict/nnetwork_interface.hh"
#include "c4/ParallelUnitTest.hh"
#include <string>
#include <vector>

// cSpell:ignore tstnnetwork

namespace {

using namespace rtt_dsxx;
using namespace rtt_c4;
using namespace rtt_predict;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_replication(ParallelUnitTest &ut) {

  std::string const nn_file_name = "kde";
  const int input_dim{100};
  const int output_dim{1};
  nnetwork_interface net(nn_file_name);
  FAIL_IF_NOT(net.valid());
  std::vector<float> input(input_dim, 1.0F / static_cast<float>(input_dim));
  std::vector<float> result = net.predict(input, input_dim, output_dim);
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(double(result[0]), 0.4872, 1.0e-4));
  if (ut.numFails == 0) {
    PASSMSG("KDE checks pass");
  } else {
    FAILMSG("KDE checks failed");
  }
}

} // namespace

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
