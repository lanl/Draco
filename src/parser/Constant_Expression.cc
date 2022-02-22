//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/Constant_Expression.cc
 * \author Kent Budge
 * \date   Wed Jul 26 07:53:32 2006
 * \brief  Definition of methods of class Constant_Expression
 * \note   Copyright (C) 2006-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Constant_Expression.hh"
#include <cmath>
#include <limits>

namespace {
using namespace std;

//------------------------------------------------------------------------------------------------//
void upper(ostream &out, char const *const label, double const power, bool &first,
           unsigned &icount) {
  if (power > 0.0) {
    if (!first) {
      out << '*';
    }
    first = false;
    if (rtt_dsxx::soft_equiv(power, trunc(power), std::numeric_limits<double>::epsilon())) {
      out << label;
      for (unsigned i = 1; i < power; ++i) {
        out << '*' << label;
      }
    } else {
      out << "pow(" << label << "," << power << ")";
    }
  } else if (power < 0.0) {
    icount++;
  }
}

//------------------------------------------------------------------------------------------------//
void lower(ostream &out, char const *const label, double const power, bool &first) {
  if (power < 0.0) {
    if (!first) {
      out << '*';
    }
    first = false;
    auto ipower = static_cast<unsigned>(-power);
    out << "pow(" << label << "," << ipower << ")";
  }
}

} // namespace

namespace rtt_parser {

//------------------------------------------------------------------------------------------------//
void write_c(Unit const &u, ostream &out) {
  double const eps = std::numeric_limits<double>::epsilon();
  double const mrv = std::numeric_limits<double>::min();

  unsigned count = !rtt_dsxx::soft_equiv(u.conv, 1.0, eps);
  double p(0.0);
  if (!rtt_dsxx::soft_equiv(u.m, 0.0, mrv)) {
    count++;
    p = u.m;
  }
  if (!rtt_dsxx::soft_equiv(u.kg, 0.0, mrv)) {
    count++;
    p = u.kg;
  }
  if (!rtt_dsxx::soft_equiv(u.s, 0.0, mrv)) {
    count++;
    p = u.s;
  }
  if (!rtt_dsxx::soft_equiv(u.A, 0.0, mrv)) {
    count++;
    p = u.A;
  }
  if (!rtt_dsxx::soft_equiv(u.K, 0.0, mrv)) {
    count++;
    p = u.K;
  }
  if (!rtt_dsxx::soft_equiv(u.mol, 0.0, mrv)) {
    count++;
    p = u.mol;
  }
  if (!rtt_dsxx::soft_equiv(u.cd, 0.0, mrv)) {
    count++;
    p = u.cd;
  }
  if (!rtt_dsxx::soft_equiv(u.rad, 0.0, mrv)) {
    count++;
    p = u.rad;
  }
  if (!rtt_dsxx::soft_equiv(u.sr, 0.0, mrv)) {
    count++;
    p = u.sr;
  }

  Require(count != 0); // should not come here if dimensionless

  if (count == 1) {
    if (p < 0 && rtt_dsxx::soft_equiv(ceil(p), p, eps)) {
      out << '(';
    }
  }

  bool first = true;
  if (!rtt_dsxx::soft_equiv(u.conv, 1.0, eps)) {
    out << u.conv;
    first = false;
  }
  unsigned icount = 0;
  upper(out, "m", u.m, first, icount);
  upper(out, "kg", u.kg, first, icount);
  upper(out, "s", u.s, first, icount);
  upper(out, "A", u.A, first, icount);
  upper(out, "K", u.K, first, icount);
  upper(out, "mol", u.mol, first, icount);
  upper(out, "cd", u.cd, first, icount);
  upper(out, "rad", u.rad, first, icount);
  upper(out, "sr", u.sr, first, icount);
  if (first) {
    out << '1';
  }
  first = false;
  if (icount > 0) {
    out << '/';
    if (icount > 1) {
      out << '(';
    }
    lower(out, "m", u.m, first);
    lower(out, "kg", u.kg, first);
    lower(out, "s", u.s, first);
    lower(out, "A", u.A, first);
    lower(out, "K", u.K, first);
    lower(out, "mol", u.mol, first);
    lower(out, "cd", u.cd, first);
    lower(out, "rad", u.rad, first);
    lower(out, "sr", u.sr, first);
    if (icount > 1) {
      out << ')';
    }
  }
  if (count == 1) {
    if (p < 0 && rtt_dsxx::soft_equiv(ceil(p), p, eps)) {
      out << ')';
    }
  }
}

} // end namespace rtt_parser

//------------------------------------------------------------------------------------------------//
// end of parser/Constant_Expression.cc
//------------------------------------------------------------------------------------------------//
