//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_analytic/test/tstPseudo_Line_Analytic_MultigroupOpacity.cc
 * \author Kent G. Budge
 * \date   Tue Apr  5 09:01:03 2011
 * \note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "cdi_analytic/Pseudo_Line_Analytic_MultigroupOpacity.hh"
#include "ds++/Release.hh"
#include "parser/Constant_Expression.hh"
#include "parser/String_Token_Stream.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_cdi_analytic;
using namespace rtt_parser;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

Pseudo_Line_Analytic_MultigroupOpacity::Averaging const NONE =
    Pseudo_Line_Analytic_MultigroupOpacity::NONE;

Pseudo_Line_Analytic_MultigroupOpacity::Averaging const ROSSELAND =
    Pseudo_Line_Analytic_MultigroupOpacity::ROSSELAND;

Pseudo_Line_Analytic_MultigroupOpacity::Averaging const PLANCK =
    Pseudo_Line_Analytic_MultigroupOpacity::PLANCK;

void tstPseudo_Line_Analytic_MultigroupOpacity(UnitTest &ut) {
  unsigned const NG = 12288; // 10;
  int const number_of_lines = 20;
  int const number_of_edges = 10;
  unsigned seed = 1;

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
    Pseudo_Line_Analytic_MultigroupOpacity model(group_bounds, rtt_cdi::ABSORPTION, continuum,
                                                 number_of_lines, peak, width, number_of_edges,
                                                 edge_ratio, 1.0, 0.0, emin, emax, NONE, 0, seed);

    ut.passes("Created Pseudo_Line_Analytic_MultigroupOpacity");

    vector<double> sigma = model.getOpacity(1.0, 1.0);

    ofstream out("pseudo_none.dat");
    for (unsigned gg = 0; gg < NG; ++gg) {
      out << ((gg + 0.5) * (emax - emin) / NG + emin) << ' ' << sigma[gg] << endl;
    }

    ut.passes("Calculated Pseudo_Line_Analytic_MultigroupOpacity opacity");
  }

  {
    Pseudo_Line_Analytic_MultigroupOpacity model(
        group_bounds, rtt_cdi::ABSORPTION, continuum, number_of_lines, peak, width, number_of_edges,
        edge_ratio, 1.0, 0.0, emin, emax, ROSSELAND, 0, seed);

    ut.passes("Created Pseudo_Line_Analytic_MultigroupOpacity");

    vector<double> sigma = model.getOpacity(1.0, 1.0);

    ofstream out("pseudo_rosseland.dat");
    for (unsigned gg = 0; gg < NG; ++gg) {
      out << ((gg + 0.5) * (emax - emin) / NG + emin) << ' ' << sigma[gg] << endl;
    }

    ut.passes("Calculated Pseudo_Line_Analytic_MultigroupOpacity opacity");
  }

  {
    Pseudo_Line_Analytic_MultigroupOpacity model(group_bounds, rtt_cdi::ABSORPTION, continuum,
                                                 number_of_lines, peak, width, number_of_edges,
                                                 edge_ratio, 1.0, 0.0, emin, emax, PLANCK, 0, seed);

    ut.passes("Created Pseudo_Line_Analytic_MultigroupOpacity");

    vector<double> sigma = model.getOpacity(1.0, 1.0);

    ofstream out("pseudo_planck.dat");
    for (unsigned gg = 0; gg < NG; ++gg) {
      out << ((gg + 0.5) * (emax - emin) / NG + emin) << ' ' << sigma[gg] << endl;
    }

    ut.passes("Calculated Pseudo_Line_Analytic_MultigroupOpacity opacity");
  }

  {
    Pseudo_Line_Analytic_MultigroupOpacity model(
        group_bounds, rtt_cdi::ABSORPTION, continuum, number_of_lines, peak, width, number_of_edges,
        edge_ratio, 1.0, 0.0, emin, emax, ROSSELAND, 1, seed);

    ut.passes("Created Pseudo_Line_Analytic_MultigroupOpacity");

    vector<double> sigma = model.getOpacity(1.0, 1.0);

    ofstream out("pseudo_rosseland_d.dat");
    for (unsigned fg = 0; fg < NG; ++fg) {
      out << ((fg + 0.5) * (emax - emin) / NG + emin) << ' ' << sigma[fg] << endl;
    }

    ut.passes("Calculated Pseudo_Line_Analytic_MultigroupOpacity opacity");
  }

  {
    Pseudo_Line_Analytic_MultigroupOpacity model(group_bounds, rtt_cdi::ABSORPTION, continuum,
                                                 number_of_lines, peak, width, number_of_edges,
                                                 edge_ratio, 1.0, 0.0, emin, emax, PLANCK, 1, seed);

    ut.passes("Created Pseudo_Line_Analytic_MultigroupOpacity");

    vector<double> sigma = model.getOpacity(1.0, 1.0);

    ofstream out("pseudo_planck_d.dat");
    for (unsigned gg = 0; gg < NG; ++gg) {
      out << ((gg + 0.5) * (emax - emin) / NG + emin) << ' ' << sigma[gg] << endl;
    }

    ut.passes("Calculated Pseudo_Line_Analytic_MultigroupOpacity opacity");
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, release);
  try {
    tstPseudo_Line_Analytic_MultigroupOpacity(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstPseudo_Line_Analytic_MultigroupOpacity.cc
//------------------------------------------------------------------------------------------------//
