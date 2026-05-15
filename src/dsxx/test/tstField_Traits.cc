//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   dsxx/test/tstField_Traits.cc
 * \author Kent Budge
 * \note   Copyright (C) 2024-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Field_Traits.hh"
#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include <complex>
#include <limits>

namespace {

using namespace std;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstFT(UnitTest &ut) {
  FAIL_IF_NOT(Field_Traits<complex<double>>::zero() == 0.0);
  FAIL_IF_NOT(Field_Traits<complex<double>>::one() == 1.0);
  double const x = 3.7;
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(value(x), 3.7));

  double constexpr eps = std::numeric_limits<double>::epsilon();
  double constexpr mrv = std::numeric_limits<double>::min();

  FAIL_IF_NOT(rtt_dsxx::soft_equiv(Field_Traits<double const>::zero(), 0.0, mrv));
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(Field_Traits<double const>::one(), 1.0, eps));
  PASSMSG("done with tstFT");
  return;
}

//------------------------------------------------------------------------------------------------//
struct unlabeled {
  int i;
};

struct labeled {
  unlabeled s;
  int j;
  operator unlabeled &() { return s; } // NOLINT [hicpp-explicit-conversions]
};

} // namespace

//------------------------------------------------------------------------------------------------//
namespace rtt_dsxx {

template <> class Field_Traits<labeled> {
public:
  using unlabeled_type = unlabeled;
};
} // namespace rtt_dsxx

namespace {

bool operator==(unlabeled const &a, labeled const &b) { return a.i == b.s.i; }

//------------------------------------------------------------------------------------------------//
void tstvalue(UnitTest &ut) {
  double x = 3;
  double const cx = 4;
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(x, value(x)) && rtt_dsxx::soft_equiv(cx, value(cx)));
  labeled s = {.s = {1}, .j = 2};
  FAIL_IF_NOT(value(s) == s);
  PASSMSG("Done with tstvalue");
  return;
}

} // namespace

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
