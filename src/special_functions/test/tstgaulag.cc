//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/test/tstgaulag.cc
 * \author Kent Budge
 * \date   Tue Sep 27 12:49:39 2005
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/DracoMath.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "special_functions/gaulag.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_special_functions;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstgaulag(UnitTest &ut) {
  vector<double> x, w;
  gaulag(x, w, 0.0, 3);
  double sum = 0.0;
  for (unsigned i = 0; i < 3; ++i)
    sum += x[i] * square(square(x[i])) * w[i];

  if (!soft_equiv(sum, 120.0))
    ut.failure("gaulag NOT accurate");
  else
    ut.passes("gaulag accurate");
  return;
}

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tstgaulag(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstgaulag.cc
//------------------------------------------------------------------------------------------------//
