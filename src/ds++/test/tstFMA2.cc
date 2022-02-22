//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstFMA.cc
 * \author Kelly Thompson <kgt@lanl.gov>
 * \date   Wed Nov 10 09:35:09 2010
 * \brief  Test functions defined in ds++/FMA.cc
 * \note   Copyright (C) 2017-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include <sstream>

#define FP_ACCURATE_FMA 0
#include "ds++/FMA.hh"

//------------------------------------------------------------------------------------------------//
// Test 2: Don't use FMA because it is too slow.
void test_fma2(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n>>> Begin test 2..." << std::endl;

  double const a(1.0e-16), b(1.0e16), c(-1.0);

  double const result = a * b + c;
  double const macro_fma_result = FMA(a, b, c);

  std::cout.precision(17);
  std::cout << "\nresult a*b+c      = " << result << "\n"
            << "result fma(a,b,c) = " << fma(a, b, c) << "\n"
            << "result FMA(a,b,c) = " << macro_fma_result << "\n"
            << std::endl;

#ifdef HAVE_HARDWARE_FMA

  std::cout << "Hardware FMA is available.\n";

#if FP_ACCURATE_FMA > 0
  FAILMSG("This test requires FP_ACCURATE_FMA=0.");
#endif

#ifdef FP_FAST_FMA
  std::cout << "\t#define FMA(a,b,c) fma(a,b,c)\n" << std::endl;

  if (rtt_dsxx::soft_equiv(fma(a, b, c), macro_fma_result, std::numeric_limits<double>::epsilon()))
    PASSMSG("With FP_ACCURATE_FMA=0, fma(a,b,c) == FMA(a,b,c).");
  else
    FAILMSG("With FP_ACCURATE_FMA=0, fma(a,b,c) != FMA(a,b,c).");

#else
  std::cout << "\t#define FMA(a,b,c) ((a)*(b)+(c))\n" << std::endl;

  if (rtt_dsxx::soft_equiv(result, macro_fma_result, std::numeric_limits<double>::epsilon()))
    PASSMSG("With FP_ACCURATE_FMA=0, a*b+c == FMA(a,b,c).");
  else
    FAILMSG("With FP_ACCURATE_FMA=0, a*b+c != FMA(a,b,c).");

#endif

#else
  std::cout << "Hardware FMA is not available.\n"

#if FP_ACCURATE_FMA > 0
            << "\t#define FMA(a,b,c) fma(a,b,c)\n"
#else
            << "\t#define FMA(a,b,c) ((a)*(b)+(c))\n"
#endif
            << std::endl;

  if (rtt_dsxx::soft_equiv(result, macro_fma_result, std::numeric_limits<double>::epsilon()))
    PASSMSG("With FP_ACCURATE_FMA=0, a*b+c == FMA(a,b,c).");
  else
    FAILMSG("With FP_ACCURATE_FMA=0, a*b+c != FMA(a,b,c).");

#endif

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_fma2(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstFMA.cc
//------------------------------------------------------------------------------------------------//
