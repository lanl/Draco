//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstDracoMath.cc
 * \author Kent G. Budge
 * \date   Wed Nov 10 09:35:09 2010
 * \brief  Test functions defined in ds++/draco_math.hh.
 * \note   Copyright (C) 2013-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "ds++/DracoMath.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstconj(rtt_dsxx::UnitTest &ut) {
  if (rtt_dsxx::soft_equiv(rtt_dsxx::conj(3.5), 3.5))
    PASSMSG("conj(double) is correct");
  else
    FAILMSG("conj(double) is NOT correct");

  std::complex<double> c(2.7, -1.4);
  if (rtt_dsxx::soft_equiv((rtt_dsxx::conj(c) * c).real(), rtt_dsxx::square(abs(c))))
    PASSMSG("conj(std::complex) is correct");
  else
    FAILMSG("conj(std::complex) is NOT correct");
  return;
}

//------------------------------------------------------------------------------------------------//
void tstcube(rtt_dsxx::UnitTest &ut) {
  if (rtt_dsxx::soft_equiv(rtt_dsxx::cube(2.0), 8.0))
    PASSMSG("rtt_dsxx::square function returned correct double");
  else
    FAILMSG("rtt_dsxx::square function did NOT return correct double.");
  return;
}

//------------------------------------------------------------------------------------------------//
void tstpythag(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::pythag;
  if (rtt_dsxx::soft_equiv(pythag(3.0e307, 4.0e307), 5.0e307))
    PASSMSG("pythag correct");
  else
    FAILMSG("pythag NOT correct");
  if (rtt_dsxx::soft_equiv(pythag(4.0e307, 3.0e307), 5.0e307))
    PASSMSG("pythag correct");
  else
    FAILMSG("pythag NOT correct");
  if (rtt_dsxx::soft_equiv(pythag(0.0, 0.0), 0.0))
    PASSMSG("pythag correct");
  else
    FAILMSG("pythag NOT correct");
  return;
}

//------------------------------------------------------------------------------------------------//
void tstsign(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::sign;
  if (!rtt_dsxx::soft_equiv(sign(3.2, 5.6), 3.2))
    FAILMSG("sign: FAILED");
  else
    PASSMSG("sign: passed");
  if (!rtt_dsxx::soft_equiv(sign(4.1, -0.3), -4.1))
    FAILMSG("sign: FAILED");
  else
    PASSMSG("sign: passed");
  return;
}

//------------------------------------------------------------------------------------------------//
void tstsquare(rtt_dsxx::UnitTest &ut) {
  if (rtt_dsxx::soft_equiv(rtt_dsxx::square(3.0), 9.0))
    PASSMSG("square function returned correct double");
  else
    FAILMSG("square function did NOT return correct double.");
  return;
}
//------------------------------------------------------------------------------------------------//
void tstceilintdiv(rtt_dsxx::UnitTest &ut) {
  using namespace rtt_dsxx;
  uint32_t const nf = ut.numFails;
  FAIL_IF_NOT(rtt_dsxx::ceil_int_division(1, 2) == 1);
  FAIL_IF_NOT(ceil_int_division(2, 2) == 1);
  FAIL_IF_NOT(ceil_int_division(0, 2) == 0);
  FAIL_IF_NOT(ceil_int_division(1, 200) == 1);
  FAIL_IF_NOT(ceil_int_division(-1, 2) == 0);
  if (ut.numFails == nf)
    PASSMSG("Fast ceiling integer division checks ok.");
  else
    FAILMSG("Fast ceiling integer division checks fail.");
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tstconj(ut);
    tstcube(ut);
    tstpythag(ut);
    tstsign(ut);
    tstsquare(ut);
    tstceilintdiv(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstDracoMath.cc
//------------------------------------------------------------------------------------------------//
