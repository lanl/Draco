//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/test/tstUnit.cc
 * \author Kent G. Budge
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "parser/Unit.hh"
#include <iostream>
#include <ostream>
#include <sstream>

namespace {

using namespace std;
using namespace rtt_parser;
using namespace rtt_dsxx;

// NOLINTBEGIN(modernize-use-designated-initializers)

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void unit_test(UnitTest &ut) {
  Unit const tstC = {0, 0, 1, 1, 0, 0, 0, 0, 0, 1};
  if (C != tstC)
    FAILMSG("unit C does NOT have expected dimensions");
  else
    PASSMSG("unit C has expected dimensions");

  Unit const tstHz = {0, 0, -1, 0, 0, 0, 0, 0, 0, 1};
  if (Hz != tstHz)
    FAILMSG("unit Hz does NOT have expected dimensions");
  else
    PASSMSG("unit Hz has expected dimensions");

  Unit const tstN = {1, 1, -2, 0, 0, 0, 0, 0, 0, 1};
  if (N != tstN)
    FAILMSG("unit N does NOT have expected dimensions");
  else
    PASSMSG("unit N has expected dimensions");

  Unit const tstJ = {2, 1, -2, 0, 0, 0, 0, 0, 0, 1};
  if (J != tstJ)
    FAILMSG("unit J does NOT have expected dimensions");
  else
    PASSMSG("unit J has expected dimensions");

  {
    ostringstream buffer;
    buffer << tstJ;
    if (buffer.str() == "1 m^2-kg-s^-2")
      PASSMSG("correct text representation of J");
    else {
      FAILMSG("NOT correct text representation of J");
      std::cout << "\tfound buffer.str() = " << buffer.str() << "\n\texpected '1 m^2-kg-s^-2'\n";
    }
  }

  Unit const tstinch = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0.0254};
  if (inch == tstinch)
    PASSMSG("unit inch has expected dimensions");
  else
    FAILMSG("unit inch does NOT have expected dimensions");

  // Test divisor
  {
    Unit const five_cm = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0.05};
    Unit const one_cm = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0.01};
    if (five_cm / 5.0 == one_cm) {
      PASSMSG("Units: Division by a constant works.");
    } else {
      ostringstream buffer;
      buffer << "Units: Division by a constant fails.\n\t"
             << "five_cm/5.0 = " << five_cm / 5.0 << "\n\t"
             << "one_cm = " << one_cm << "\n\t"
             << "test five_cm/5.0 == one_cm ?\n";
      FAILMSG(buffer.str());
    }
  }

  {
    Unit const t1 = {0, 0, 1, 1, 0, 0, 0, 0, 0, 1};
    Unit const t2 = {0.1, 0, 1, 1, 0, 0, 0, 0, 0, 1};
    if (t1 == t2)
      FAILMSG("failed to detect length difference");
    if (is_compatible(t1, t2))
      FAILMSG("failed to detect length difference");
  }
  {
    Unit const t1 = {0, 0, 1, 1, 0, 0, 0, 0, 0, 1};
    Unit const t2 = {0, 0.1, 1, 1, 0, 0, 0, 0, 0, 1};
    if (t1 == t2)
      FAILMSG("failed to detect mass difference");
    if (is_compatible(t1, t2))
      FAILMSG("failed to detect mass difference");
  }
  {
    Unit const t1 = {0, 0, 1, 1, 0, 0, 0, 0, 0, 1};
    Unit const t2 = {0, 0, 1.1, 1, 0, 0, 0, 0, 0, 1};
    if (t1 == t2)
      FAILMSG("failed to detect time difference");
    if (is_compatible(t1, t2))
      FAILMSG("failed to detect time difference");
  }
  {
    Unit const t1 = {0, 0, 1, 1, 0, 0, 0, 0, 0, 1};
    Unit const t2 = {0, 0, 1, 1.1, 0, 0, 0, 0, 0, 1};
    if (t1 == t2)
      FAILMSG("failed to detect current difference");
    if (is_compatible(t1, t2))
      FAILMSG("failed to detect current difference");
  }
  {
    Unit const t1 = {0, 0, 1, 1, 0, 0, 0, 0, 0, 1};
    Unit const t2 = {0, 0, 1, 1, 0.1, 0, 0, 0, 0, 1};
    if (t1 == t2)
      FAILMSG("failed to detect temperature difference");
    if (is_compatible(t1, t2))
      FAILMSG("failed to detect temperature difference");
  }
  {
    Unit const t1 = {0, 0, 1, 1, 0, 0, 0, 0, 0, 1};
    Unit const t2 = {0, 0, 1, 1, 0, 0.1, 0, 0, 0, 1};
    if (t1 == t2)
      FAILMSG("failed to detect mole difference");
    if (is_compatible(t1, t2))
      FAILMSG("failed to detect mole difference");
  }
  {
    Unit const t1 = {0, 0, 1, 1, 0, 0, 0, 0, 0, 1};
    Unit const t2 = {0, 0, 1, 1, 0, 0, 0.1, 0, 0, 1};
    if (t1 == t2)
      FAILMSG("failed to detect luminance difference");
    if (is_compatible(t1, t2))
      FAILMSG("failed to detect luminance difference");
  }
  {
    Unit const t1 = {0, 0, 1, 1, 0, 0, 0, 0, 0, 1};
    Unit const t2 = {0, 0, 1, 1, 0, 0, 0, 0.1, 0, 1};
    if (t1 == t2)
      FAILMSG("failed to detect rad difference");
    if (is_compatible(t1, t2))
      FAILMSG("failed to detect rad difference");
  }
  {
    Unit const t1 = {0, 0, 1, 1, 0, 0, 0, 0, 0, 1};
    Unit const t2 = {0, 0, 1, 1, 0, 0, 0, 0, 0.1, 1};
    if (t1 == t2)
      FAILMSG("failed to detect solid angle difference");
    if (is_compatible(t1, t2))
      FAILMSG("failed to detect solid angle difference");
  }
  {
    Unit const t1 = {0, 0, 1, 1, 0, 0, 0, 0, 0, 1};
    Unit const t2 = {0, 0, 1, 1, 0, 0, 0, 0, 0, 1.1};
    if (t1 == t2)
      FAILMSG("failed to detect conversion factor difference");
    if (!is_compatible(t1, t2))
      FAILMSG("failed to ignore conversion factor difference");
  }
  if (soft_equiv(conversion_factor(J, CGS), erg.conv))
    PASSMSG("Conversion to CGS good");
  else
    FAILMSG("Conversion to CGS NOT good");
}

// NOLINTEND(modernize-use-designated-initializers)

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    unit_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstUnit.cc
//------------------------------------------------------------------------------------------------//
