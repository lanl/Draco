//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/test/tstField_Traits.cc
 * \author Kent Budge
 * \date   Tue Aug 26 12:18:55 2008
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "ds++/Field_Traits.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"

using namespace std;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstFT(UnitTest &ut) {
  if (Field_Traits<complex<double>>::zero() == 0.0)
    PASSMSG("complex zero good");
  else
    FAILMSG("complex zero NOT good");
  if (Field_Traits<complex<double>>::one() == 1.0)
    PASSMSG("complex zero good");
  else
    FAILMSG("complex zero NOT good");
  double const x = 3.7;
  if (rtt_dsxx::soft_equiv(value(x), 3.7))
    PASSMSG("complex zero good");
  else
    FAILMSG("complex zero NOT good");

  double const eps = std::numeric_limits<double>::epsilon();
  double const mrv = std::numeric_limits<double>::min();

  if (rtt_dsxx::soft_equiv(Field_Traits<double const>::zero(), 0.0, mrv))
    PASSMSG("double zero good");
  else
    FAILMSG("double zero NOT good");
  if (rtt_dsxx::soft_equiv(Field_Traits<double const>::one(), 1.0, eps))
    PASSMSG("double zero good");
  else
    FAILMSG("double zero NOT good");
  return;
}

//------------------------------------------------------------------------------------------------//
struct unlabeled {
  int i;
};

struct labeled {
  unlabeled s;
  int j;

  operator unlabeled &() { return s; }
};

namespace rtt_dsxx {

template <> class Field_Traits<labeled> {
public:
  using unlabeled_type = unlabeled;
};
} // namespace rtt_dsxx

bool operator==(unlabeled const &a, labeled const &b) { return a.i == b.s.i; }

void tstvalue(UnitTest &ut) {
  double x = 3;
  double const cx = 4;

  if (rtt_dsxx::soft_equiv(x, value(x)) && rtt_dsxx::soft_equiv(cx, value(cx)))
    PASSMSG("value strips double correctly");
  else
    FAILMSG("value does NOT strip double correctly");

  labeled s = {{1}, 2};

  if (value(s) == s)
    PASSMSG("value strips struct correctly");
  else
    FAILMSG("value does NOT strip struct correctly");

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstFT(ut);
    tstvalue(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstField_Traits.cc
//------------------------------------------------------------------------------------------------//
