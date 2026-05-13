//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_analytic/test/tstPseudo_Line_Analytic_MultigroupOpacity.cc
 * \author Kent G. Budge
 * \note   Copyright (C) 2011-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/UnitTest.hh"
#include "c4/ParallelUnitTest.hh"
#include "cdi/OpacityCommon.hh"
#include "cdi_analytic/Pseudo_Line_Analytic_MultigroupOpacity.hh"
#include "parser/String_Token_Stream.hh"
#include "parser/Unit.hh"
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace {

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_cdi_analytic;
using namespace rtt_parser;
using PLAMO = Pseudo_Line_Analytic_MultigroupOpacity;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

PLAMO::Averaging const NONE = PLAMO::NONE;
PLAMO::Averaging const ROSSELAND = PLAMO::ROSSELAND;
PLAMO::Averaging const PLANCK = PLAMO::PLANCK;

void tstPLAMO(UnitTest &ut) {
  unsigned const NG = 12288; // 10;
  int const number_of_lines = 20;
  int const number_of_edges = 10;
  unsigned constexpr seed = 1;
  double constexpr Tref = 1.0;
  double constexpr Tpow = 0.0;

  std::shared_ptr<Expression const> continuum;
  {
    map<string, pair<unsigned, Unit>> variables;
    variables["x"] = pair<unsigned, Unit>(0, raw);

    String_Token_Stream expr("1.0e-2 + 20/(x+1)^3 + 1e-4*x*x*x*x");
    continuum = Expression::parse(1, variables, expr);
  }

  double const peak = 1e1;
  double const width = 0.02; // keV
  double const edge_ratio = 10.0;
  double const emax = 10.0; // keV
  double const emin = 0.0;

  vector<double> group_bounds(NG + 1);
  for (unsigned i = 0; i <= NG; i++) {
    group_bounds[i] = i * (emax - emin) / NG + emin;
  }

  {
    PLAMO const model(group_bounds, rtt_cdi::ABSORPTION, continuum, number_of_lines, peak, width,
                      number_of_edges, edge_ratio, Tref, Tpow, emin, emax, NONE, 0, seed);

    PASSMSG("Created Pseudo_Line_Analytic_MultigroupOpacity");

    vector<double> sigma = model.getOpacity(1.0, 1.0);

    ofstream out("pseudo_none.dat");
    for (unsigned gg = 0; gg < NG; ++gg) {
      out << ((gg + 0.5) * (emax - emin) / NG + emin) << ' ' << sigma[gg] << '\n';
    }

    PASSMSG("Calculated Pseudo_Line_Analytic_MultigroupOpacity opacity");
  }

  {
    PLAMO const model(group_bounds, rtt_cdi::ABSORPTION, continuum, number_of_lines, peak, width,
                      number_of_edges, edge_ratio, Tref, Tpow, emin, emax, ROSSELAND, 0, seed);

    PASSMSG("Created Pseudo_Line_Analytic_MultigroupOpacity");

    vector<double> sigma = model.getOpacity(1.0, 1.0);

    ofstream out("pseudo_rosseland.dat");
    for (unsigned gg = 0; gg < NG; ++gg) {
      out << ((gg + 0.5) * (emax - emin) / NG + emin) << ' ' << sigma[gg] << '\n';
    }

    PASSMSG("Calculated Pseudo_Line_Analytic_MultigroupOpacity opacity");
  }

  {
    PLAMO const model(group_bounds, rtt_cdi::ABSORPTION, continuum, number_of_lines, peak, width,
                      number_of_edges, edge_ratio, Tref, Tpow, emin, emax, PLANCK, 0, seed);

    PASSMSG("Created Pseudo_Line_Analytic_MultigroupOpacity");

    vector<double> sigma = model.getOpacity(1.0, 1.0);

    ofstream out("pseudo_planck.dat");
    for (unsigned gg = 0; gg < NG; ++gg) {
      out << ((gg + 0.5) * (emax - emin) / NG + emin) << ' ' << sigma[gg] << '\n';
    }

    PASSMSG("Calculated Pseudo_Line_Analytic_MultigroupOpacity opacity");
  }

  {
    PLAMO const model(group_bounds, rtt_cdi::ABSORPTION, continuum, number_of_lines, peak, width,
                      number_of_edges, edge_ratio, Tref, Tpow, emin, emax, ROSSELAND, 1, seed);

    PASSMSG("Created Pseudo_Line_Analytic_MultigroupOpacity");

    vector<double> sigma = model.getOpacity(1.0, 1.0);

    ofstream out("pseudo_rosseland_d.dat");
    for (unsigned fg = 0; fg < NG; ++fg) {
      out << ((fg + 0.5) * (emax - emin) / NG + emin) << ' ' << sigma[fg] << '\n';
    }

    PASSMSG("Calculated Pseudo_Line_Analytic_MultigroupOpacity opacity");
  }

  {
    PLAMO const model(group_bounds, rtt_cdi::ABSORPTION, continuum, number_of_lines, peak, width,
                      number_of_edges, edge_ratio, Tref, Tpow, emin, emax, PLANCK, 1, seed);

    PASSMSG("Created Pseudo_Line_Analytic_MultigroupOpacity");

    vector<double> sigma = model.getOpacity(1.0, 1.0);

    ofstream out("pseudo_planck_d.dat");
    for (unsigned gg = 0; gg < NG; ++gg) {
      out << ((gg + 0.5) * (emax - emin) / NG + emin) << ' ' << sigma[gg] << '\n';
    }

    PASSMSG("Calculated Pseudo_Line_Analytic_MultigroupOpacity opacity");
  }
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, release);
  try {
    tstPLAMO(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstPseudo_Line_Analytic_MultigroupOpacity.cc
//------------------------------------------------------------------------------------------------//
