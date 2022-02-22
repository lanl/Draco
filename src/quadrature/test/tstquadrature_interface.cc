//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/test/tstquadrature_interface.cc
 * \author Jae Chang
 * \date   Tue Feb  3 13:08:49 2012
 * \brief  Quadrature Interface test
 * \note   Copyright (C) 2015-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "quadrature/Quadrature.hh"
#include "quadrature/Quadrature_Interface.hh"
#include <array>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_quadrature;

//------------------------------------------------------------------------------------------------//
// Unit tests
//------------------------------------------------------------------------------------------------//

void quadrature_interface_test(rtt_dsxx::UnitTest &ut) {

  int dimension = 1;
  int type = 0;
  int order = 4;
  int azimuthal_order = 2;
  int geometry = 0;
  std::vector<double> mu(30, 0.0);
  std::vector<double> eta(30, 0.0);
  std::vector<double> xi(30, 0.0);
  std::vector<double> weights(30, 0.0);

  quadrature_data quad;

  init_quadrature(quad);

  // testing Gauss_Legendre

  quad.dimension = dimension;
  quad.type = type;
  quad.order = order;
  quad.azimuthal_order = azimuthal_order;
  quad.geometry = geometry;
  quad.mu = &mu[0];
  quad.weights = &weights[0];

  get_quadrature(quad);

  std::array<double, 4> gs_wt = {0.1739274226, 0.3260725774, 0.3260725774, 0.1739274226};
  std::array<double, 4> gs_mu = {-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116};

  for (unsigned i = 0; i < 4; ++i) {
    if (!soft_equiv(quad.mu[i], gs_mu[i], 1e-8)) {
      ut.failure("incorrect mu");
      return;
    }
    if (!soft_equiv(quad.weights[i], gs_wt[i], 1e-8)) {
      ut.failure("incorrect wt");
      return;
    }
  }
  ut.passes("Gauss_Legendre passed");

  // testing Lobato
  init_quadrature(quad);

  quad.dimension = 1;
  quad.type = 1;
  quad.order = 4;
  quad.geometry = 0;
  quad.mu = &mu[0];
  quad.weights = &weights[0];

  get_quadrature(quad);

  std::array<double, 4> lobato_mu = {-1, -0.4472135955, 0.4472135955, 1};
  std::array<double, 4> lobato_wt = {0.08333333333, 0.4166666667, 0.4166666667, 0.08333333333};

  for (unsigned i = 0; i < 4; ++i) {
    if (!soft_equiv(quad.mu[i], lobato_mu[i], 1e-8)) {
      ut.failure("incorrect mu");
      return;
    }
    if (!soft_equiv(quad.weights[i], lobato_wt[i], 1e-8)) {
      ut.failure("incorrect wt");
      return;
    }
  }
  ut.passes("Lobato passed");

  // testing Level Symmetric
  init_quadrature(quad);
  quad.dimension = 2;
  quad.type = 0;
  quad.order = 4;
  quad.geometry = 1;
  quad.mu = &mu[0];
  quad.eta = &eta[0];
  quad.xi = &xi[0];
  quad.weights = &weights[0];

  get_quadrature(quad);

  std::array<double, 16> ls_mu = {-0.4950046922, -0.3500211746, 0.3500211746,  -0.936741788,
                                  -0.8688903007, -0.3500211746, 0.3500211746,  0.8688903007,
                                  -0.936741788,  -0.8688903007, -0.3500211746, 0.3500211746,
                                  0.8688903007,  -0.4950046922, -0.3500211746, 0.3500211746};
  std::array<double, 16> ls_eta = {-0.8688903007, -0.8688903007, -0.8688903007, -0.3500211746,
                                   -0.3500211746, -0.3500211746, -0.3500211746, -0.3500211746,
                                   0.3500211746,  0.3500211746,  0.3500211746,  0.3500211746,
                                   0.3500211746,  0.8688903007,  0.8688903007,  0.8688903007};
  std::array<double, 16> ls_wt = {0,
                                  0.08333333333,
                                  0.08333333333,
                                  0,
                                  0.08333333333,
                                  0.08333333333,
                                  0.08333333333,
                                  0.08333333333,
                                  0,
                                  0.08333333333,
                                  0.08333333333,
                                  0.08333333333,
                                  0.08333333333,
                                  0,
                                  0.08333333333,
                                  0.08333333333};

  for (unsigned i = 0; i < 16; ++i) {
    if (!soft_equiv(quad.mu[i], ls_mu[i], 1e-8)) {
      ut.failure("incorrect mu");
      return;
    }
    if (!soft_equiv(quad.eta[i], ls_eta[i], 1e-8)) {
      ut.failure("incorrect eta");
      return;
    }
    if (!soft_equiv(quad.weights[i], ls_wt[i], 1e-8)) {
      ut.failure("incorrect wt");
      return;
    }
  }
  ut.passes("Level Symmetric passed");

  // testing Tri Chebyshev Legendre
  init_quadrature(quad);
  quad.dimension = 2;
  quad.type = 1;
  quad.order = 4;
  quad.geometry = 0;
  quad.mu = &mu[0];
  quad.eta = &eta[0];
  quad.xi = &xi[0];
  quad.weights = &weights[0];

  get_quadrature(quad);

  std::array<double, 12> tcl_mu = {-0.3598878562, 0.3598878562, -0.8688461434, 0.8688461434,
                                   -0.3594747925, 0.3594747925, -0.3594747925, 0.3594747925,
                                   -0.8688461434, 0.8688461434, -0.3598878562, 0.3598878562};
  std::array<double, 12> tcl_eta = {-0.8688461434, -0.8688461434, -0.3598878562, -0.3598878562,
                                    -0.3594747925, -0.3594747925, 0.3594747925,  0.3594747925,
                                    0.3598878562,  0.3598878562,  0.8688461434,  0.8688461434};
  std::array<double, 12> tcl_wt = {0.08151814436, 0.08151814436, 0.08151814436, 0.08151814436,
                                   0.08696371128, 0.08696371128, 0.08696371128, 0.08696371128,
                                   0.08151814436, 0.08151814436, 0.08151814436, 0.08151814436};

  for (unsigned i = 0; i < 12; ++i) {
    if (!soft_equiv(quad.mu[i], tcl_mu[i], 1e-8)) {
      ut.failure("incorrect mu");
      return;
    }
    if (!soft_equiv(quad.eta[i], tcl_eta[i], 1e-8)) {
      ut.failure("incorrect eta");
      return;
    }
    if (!soft_equiv(quad.weights[i], tcl_wt[i], 1e-8)) {
      ut.failure("incorrect wt");
      return;
    }
  }
  ut.passes("Tri Chebyshev Legendre passed");

  // testing Square Chebyshev Legendre
  init_quadrature(quad);
  quad.dimension = 2;
  quad.type = 2;
  quad.order = 4;
  quad.geometry = 0;
  quad.mu = &mu[0];
  quad.eta = &eta[0];
  quad.xi = &xi[0];

  quad.weights = &weights[0];

  get_quadrature(quad);

  std::array<double, 16> scl_mu = {-0.3598878562, 0.3598878562, -0.1945463558, 0.1945463558,
                                   -0.8688461434, 0.8688461434, -0.4696764507, 0.4696764507,
                                   -0.4696764507, 0.4696764507, -0.8688461434, 0.8688461434,
                                   -0.1945463558, 0.1945463558, -0.3598878562, 0.3598878562};
  std::array<double, 16> scl_eta = {-0.8688461434, -0.8688461434, -0.4696764507, -0.4696764507,
                                    -0.3598878562, -0.3598878562, -0.1945463558, -0.1945463558,
                                    0.1945463558,  0.1945463558,  0.3598878562,  0.3598878562,
                                    0.4696764507,  0.4696764507,  0.8688461434,  0.8688461434};
  std::array<double, 16> scl_wt = {0.08151814436, 0.08151814436, 0.04348185564, 0.04348185564,
                                   0.08151814436, 0.08151814436, 0.04348185564, 0.04348185564,
                                   0.04348185564, 0.04348185564, 0.08151814436, 0.08151814436,
                                   0.04348185564, 0.04348185564, 0.08151814436, 0.08151814436};

  for (unsigned i = 0; i < 16; ++i) {
    if (!soft_equiv(quad.mu[i], scl_mu[i], 1e-8)) {
      ut.failure("incorrect mu");
      return;
    }
    if (!soft_equiv(quad.eta[i], scl_eta[i], 1e-8)) {
      ut.failure("incorrect eta");
      return;
    }
    if (!soft_equiv(quad.weights[i], scl_wt[i], 1e-8)) {
      ut.failure("incorrect wt");
      return;
    }
  }
  ut.passes("Square Chebyshev Legendre passed");

  // testing Product Chebyshev Legendre
  init_quadrature(quad);
  quad.dimension = 2;
  quad.type = 3;
  quad.order = 2;
  quad.azimuthal_order = 4;
  quad.geometry = 0;
  quad.mu = &mu[0];
  quad.eta = &eta[0];
  quad.xi = &xi[0];
  quad.weights = &weights[0];

  get_quadrature(quad);

  std::array<double, 8> pcl_mu = {-0.3124597141, 0.3124597141, -0.7543444795, 0.7543444795,
                                  -0.7543444795, 0.7543444795, -0.3124597141, 0.3124597141};
  std::array<double, 8> pcl_eta = {-0.7543444795, -0.7543444795, -0.3124597141, -0.3124597141,
                                   0.3124597141,  0.3124597141,  0.7543444795,  0.7543444795};
  std::array<double, 9> pcl_wt = {0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125};

  for (unsigned i = 0; i < 8; ++i) {
    if (!soft_equiv(quad.mu[i], pcl_mu[i], 1e-8)) {
      ut.failure("incorrect mu");
      return;
    }
    if (!soft_equiv(quad.eta[i], pcl_eta[i], 1e-8)) {
      ut.failure("incorrect eta");
      return;
    }
    if (!soft_equiv(quad.weights[i], pcl_wt[i], 1e-8)) {
      ut.failure("incorrect wt");
      return;
    }
  }
  ut.passes("Product Chebyshev Legendre passed");

  if (ut.numFails == 0)
    PASSMSG("test successfull.");

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    quadrature_interface_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstquadrature_interface.cc
//------------------------------------------------------------------------------------------------//
