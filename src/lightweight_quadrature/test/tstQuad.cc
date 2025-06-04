//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \brief  Test 1D and 2D quadratures
 * \note   Copyright (C) 2025 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
//      Quadrature Sets       //
//------------------------------------------------------------------------------------------------//
//  1D Polar   //  1D Azimuthal  //
//------------------------------------------------------------------------------------------------//
//    gauss    //  chebyshev  //
// doublegauss //             //
//   lobatto   //  az_qrs45   //
//    radau    //  az_qra45   //
// doubleradau //  az_qrs30   //
//             //  az_qra30   //
//  pol_qrs45  //  az_qrs90   //
//  pol_qra45  //  az_qra90   //
//  pol_qrs90  //             //
//------------------------------------------------------------------------------------------------//
//          product           //
//------------------------------------------------------------------------------------------------//
//          triangle          //
//------------------------------------------------------------------------------------------------//

#include <cmath>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"

#include "lightweight_quadrature/quad_1D.hh"
#include "lightweight_quadrature/quad_2D.hh"

using namespace rtt_lightweight_quadrature;

void print_1D_quad(std::vector<std::vector<double>> const &quad) {
  printf("\n 1D %19s %19s \n", "points ", "weights ");
  double sumwt{0.0};

  for (size_t i = 0; i < quad[0].size(); ++i) {
    printf("%5ld %19.10e %19.10e\n", i, quad[0][i], quad[1][i]);
    sumwt += quad[1][i];
  }
  printf("  normalization %19.10e\n\n", sumwt);
}

//-------------------------------------------------------------------------------------------------//

static constexpr double tol = 1.0e-12;

static constexpr int num = 6;

// integrate sin^m(x) cos^n(x) sin(x) functions
auto sin_cos_pol_powers(int N, std::vector<std::vector<double>> &quad) {
  std::vector<double> ints((num + 1) * (num + 2) / 2, 0);

  for (int i = 0; i < N; i++) {
    const double x = quad[0][i];
    const double w = quad[1][i];

    for (int m = 0; m < num + 1; ++m) // sin^m(x) sin(x)
    {
      const double f = w * std::pow(1 - x * x, m / 2.0);

      for (int n = 0; n < num / 2 + 1; ++n) // cos^(2n)(x) (even powers)
      {
        ints[n + m * (num / 2 + 1)] += f * std::pow(x, 2 * n);
      }
    }
  };

  return ints;
}

// integrate sin^m(x) cos^n(x) functions
auto sin_cos_az_powers(int N, std::vector<std::vector<double>> &quad) {
  std::vector<double> ints((num + 1) * (num + 2) / 2, 0);

  for (int i = 0; i < N; i++) {
    const double x = quad[0][i];
    const double w = quad[1][i];

    for (int m = 0; m < num + 1; ++m) // sin^m(x)
    {
      const double f = w * std::pow(1 - x * x, m / 2.0);

      for (int n = 0; n < num / 2 + 1; ++n) // cos^(2n)(x) (even powers)
      {
        ints[n + m * (num / 2 + 1)] += f * std::pow(x, 2 * n);
      }
    }
  };

  return ints;
}

// integrate cos^n(x) functions
auto cos_powers(int N, std::vector<std::vector<double>> &quad) {
  std::vector<double> ints(num + 1, 0);

  for (int n = 0; n < N; n++) {
    const double x = quad[0][n];
    const double w = quad[1][n];

    double acc = w;
    for (int j = 0; j < num + 1; ++j) {
      ints[j] += acc;
      acc *= x;
    }
  };

  return ints;
}

void test_pol_cos(int N, std::vector<std::vector<double>> &quad, rtt_dsxx::UnitTest &ut) {
  printf("polar full-range [0,pi]\n");

  // calculate integrals
  auto ints = cos_powers(N, quad);

  for (int j = 0; j < num + 1; ++j) {
    const double v = (j % 2 ? 0.0 : 2.0 / (j + 1));
    printf("integral %d %.10e, should be %.10e\n", j, ints[j], v);
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(ints[j], v, tol));
  }
}

void test_pol_cos_half(int N, std::vector<std::vector<double>> &quad, rtt_dsxx::UnitTest &ut) {
  printf("polar half-range [0,pi/2]\n");

  // calculate integrals
  auto ints = cos_powers(N / 2, quad);

  for (int j = 0; j < num + 1; ++j) {
    const double v = 1.0 / (j + 1);
    printf("integral %d %.10e, should be %.10e\n", j, ints[j], v);
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(ints[j], v, tol));
  }
}

//-------------------------------------------------------------------------------------------------//

void test_az_cos(int N, std::vector<std::vector<double>> &quad, rtt_dsxx::UnitTest &ut) {
  printf("azimuthal full-range [0,pi]\n");

  // calculate integrals
  auto ints = cos_powers(N, quad);

  auto sqi = [](double &&n) { return gsl_sf_gamma(n + 0.5) / gsl_sf_gamma(n + 1); };

  for (int j = 0; j < num + 1; ++j) {
    const double v = (j % 2 ? 0.0 : sqi(static_cast<double>(j) / 2)) * std::sqrt(cPI);
    printf("integral %d %.10e, should be %.10e\n", j, ints[j], v);
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(ints[j], v, tol));
  }
}

void test_az_cos_half(int N, std::vector<std::vector<double>> &quad, rtt_dsxx::UnitTest &ut) {
  printf("azimuthal half-range [0,pi/2]\n");

  // calculate integrals
  auto ints = cos_powers(N / 2, quad);

  auto sqi = [](double &&n) { return gsl_sf_gamma(n + 0.5) / gsl_sf_gamma(n + 1) / 2.0; };

  for (int j = 0; j < num + 1; ++j) {
    const double v = sqi(static_cast<double>(j) / 2) * std::sqrt(cPI);
    printf("integral %d %.10e, should be %.10e\n", j, ints[j], v);
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(ints[j], v, tol));
  }
}

void test_az_cos_third(int N, std::vector<std::vector<double>> &quad, rtt_dsxx::UnitTest &ut) {
  printf("azimuthal half-range [0,pi/2]\n");

  // calculate integrals
  auto ints = cos_powers(N / 3, quad);

  auto sqi = [](double &&n) { return gsl_sf_gamma(n + 0.5) / gsl_sf_gamma(n + 1) / 2.0; };

  for (int j = 0; j < num + 1; ++j) {
    const double v = sqi(static_cast<double>(j) / 2) * std::sqrt(cPI) / 3;
    printf("integral %d %.10e, should be %.10e\n", j, ints[j], v);
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(ints[j], v, tol));
  }
}

void test_az_sin_cos(int N, std::vector<std::vector<double>> &quad, rtt_dsxx::UnitTest &ut) {
  printf("azimuthal full-range [0,pi]\n");

  // calculate integrals
  auto ints = sin_cos_az_powers(N, quad);

  auto sqi = [](double &&m, double &&n) {
    return gsl_sf_gamma(n + 0.5) * gsl_sf_gamma(m + 0.5) / gsl_sf_gamma(n + m + 1);
  };

  for (int m = 0; m < num + 1; ++m) {
    for (int n = 0; n < num / 2 + 1; ++n) {
      const double v = sqi(static_cast<double>(m) / 2, static_cast<double>(n));
      printf("integral %d %.10e, should be %.10e\n", n + m * (num / 2 + 1),
             ints[n + m * (num / 2 + 1)], v);
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(ints[n + m * (num / 2 + 1)], v, tol));
    }
  }
}

void test_az_sin_cos_half(int N, std::vector<std::vector<double>> &quad, rtt_dsxx::UnitTest &ut) {
  printf("azimuthal half-range [0,pi/2]\n");

  // calculate integrals
  auto ints = sin_cos_az_powers(N / 2, quad);

  auto sqi = [](double &&m, double &&n) {
    return gsl_sf_gamma(n + 0.5) * gsl_sf_gamma(m + 0.5) / gsl_sf_gamma(n + m + 1);
  };

  for (int m = 0; m < num + 1; ++m) {
    for (int n = 0; n < num / 2 + 1; ++n) {
      const double v = sqi(static_cast<double>(m) / 2, static_cast<double>(n)) / 2;
      printf("integral %d %.10e, should be %.10e\n", n + m * (num / 2 + 1),
             ints[n + m * (num / 2 + 1)], v);
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(ints[n + m * (num / 2 + 1)], v, tol));
    }
  }
}

void test_az_sin_cos_third(int N, std::vector<std::vector<double>> &quad, rtt_dsxx::UnitTest &ut) {
  printf("azimuthal third-range [0,pi/3]\n");

  // calculate integrals
  auto ints = sin_cos_az_powers(N / 3, quad);

  auto sqi = [](double &&m, double &&n) {
    return std::pow(3, m / 2 + 1.5) * std::pow(2, -(m + 1)) *
           gsl_sf_hyperg_2F1(m / 2 + 0.5, 0.5 - n / 2, m / 2 + 1.5, 0.75) / 3 / (m + 1);
  };

  for (int m = 0; m < num + 1; ++m) {
    for (int n = 0; n < num / 2 + 1; ++n) {
      const double v = sqi(static_cast<double>(m), 2.0 * n);
      printf("integral %d %.10e, should be %.10e\n", n + m * (num / 2 + 1),
             ints[n + m * (num / 2 + 1)], v);
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(ints[n + m * (num / 2 + 1)], v, tol));
    }
  }
}

void test_pol_sin_cos(int N, std::vector<std::vector<double>> &quad, rtt_dsxx::UnitTest &ut) {
  printf("polar full-range [0,pi]\n");

  // calculate integrals
  auto ints = sin_cos_pol_powers(N, quad);

  auto sqi = [](double &&m, double &&n) {
    return gsl_sf_gamma(n + 0.5) * gsl_sf_gamma(m + 1) / gsl_sf_gamma(n + m + 1.5);
  };

  for (int m = 0; m < num + 1; ++m) {
    for (int n = 0; n < num / 2 + 1; ++n) {
      const double v = sqi(static_cast<double>(m) / 2, static_cast<double>(n));
      printf("integral %d %.10e, should be %.10e\n", n + m * (num / 2 + 1),
             ints[n + m * (num / 2 + 1)], v);
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(ints[n + m * (num / 2 + 1)], v, tol));
    }
  }
}

void test_pol_sin_cos_half(int N, std::vector<std::vector<double>> &quad, rtt_dsxx::UnitTest &ut) {
  printf("polar half-range [0,pi/2]\n");

  // calculate integrals
  auto ints = sin_cos_pol_powers(N / 2, quad);

  auto sqi = [](double &&m, double &&n) {
    return gsl_sf_gamma(n + 0.5) * gsl_sf_gamma(m + 1) / gsl_sf_gamma(n + m + 1.5) / 2;
  };

  for (int m = 0; m < num + 1; ++m) {
    for (int n = 0; n < num / 2 + 1; ++n) {
      const double v = sqi(static_cast<double>(m) / 2, static_cast<double>(n));
      printf("integral %d %.10e, should be %.10e\n", n + m * (num / 2 + 1),
             ints[n + m * (num / 2 + 1)], v);
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(ints[n + m * (num / 2 + 1)], v, tol));
    }
  }
}

//-------------------------------------------------------------------------------------------------//

// 1D quadrature tests
void quad_test_1D(int const N, rtt_dsxx::UnitTest &ut) {
  std::vector<std::vector<double>> quad;

  // Classical full-range and double-range

  std::cout << "\nGauss-Chebyshev\n";
  quad = chebyshev(N);
  print_1D_quad(quad);
  test_az_cos(N, quad, ut);

  std::cout << "\nGauss-Legendre\n";
  quad = gauss(N);
  print_1D_quad(quad);
  test_pol_cos(N, quad, ut);

  std::cout << "\nDouble Gauss-Legendre\n";
  quad = doublegauss(N);
  print_1D_quad(quad);
  test_pol_cos(N, quad, ut);
  test_pol_cos_half(N, quad, ut);

  std::cout << "\nGauss-Radau\n";
  quad = radau(N);
  print_1D_quad(quad);
  test_pol_cos(N, quad, ut);

  std::cout << "\nDouble Radau\n";
  quad = doubleradau(N);
  print_1D_quad(quad);
  test_pol_cos(N, quad, ut);
  test_pol_cos_half(N, quad, ut);

  std::cout << "\nGauss-Lobatto\n";
  quad = lobatto(N);
  print_1D_quad(quad);
  test_pol_cos(N, quad, ut);

  // Non-classical quadruple-range
  // Spence (QRS) and Abu-Shumays (QRA)

  std::cout << "\nAzimuthal QRS90\n";
  quad = az_qrs90(N);
  print_1D_quad(quad);
  test_az_cos(N, quad, ut);
  test_az_sin_cos(N, quad, ut);
  test_az_sin_cos_half(N, quad, ut);

  std::cout << "\nAzimuthal2D QRA90\n";
  quad = az2D_qra90(N);
  print_1D_quad(quad);
  test_az_cos(N, quad, ut);
  test_az_cos_half(N, quad, ut);

  std::cout << "\nAzimuthal QRS45\n";
  quad = az_qrs45(N);
  print_1D_quad(quad);
  test_az_cos(N, quad, ut);
  test_az_sin_cos(N, quad, ut);
  test_az_sin_cos_half(N, quad, ut);

  std::cout << "\nAzimuthal QRA45\n";
  quad = az_qra45(N);
  print_1D_quad(quad);
  test_az_cos(N, quad, ut);
  test_az_cos_half(N, quad, ut);
  test_az_sin_cos(N, quad, ut);

  if (!(N % 3)) {
    std::cout << "\nAzimuthal QRS30\n";
    quad = az_qrs30(N);
    print_1D_quad(quad);
    test_az_cos(N, quad, ut);
    test_az_sin_cos_third(N, quad, ut);

    std::cout << "\nAzimuthal QRA30\n";
    quad = az_qra30(N);
    print_1D_quad(quad);
    test_az_cos(N, quad, ut);
    test_az_sin_cos_third(N, quad, ut);
  }

  std::cout << "\nPolar QRS90\n";
  quad = pol_qrs90(N);
  print_1D_quad(quad);
  test_pol_cos(N, quad, ut);
  test_pol_sin_cos(N, quad, ut);
  test_pol_sin_cos_half(N, quad, ut);

  std::cout << "\nPolar QRS45\n";
  quad = pol_qrs45(N);
  print_1D_quad(quad);
  test_pol_cos(N, quad, ut);
  test_pol_cos_half(N, quad, ut);
  test_pol_sin_cos_half(N, quad, ut);
}

//-------------------------------------------------------------------------------------------------//

// 2D quadrature tests
void test_sph(quad_2D const &quad, bool const prod_type, rtt_dsxx::UnitTest &ut) {
  std::vector<double> ox, oy, oz, wt;
  quad.create_quadrature(ox, oy, oz, wt);

  double sum0 = 0;
  double sum1 = 0;
  double sum2 = 0;
  double sum3 = 0;
  double sum4 = 0;
  double sum5 = 0;

  for (size_t i = 0; i < wt.size(); ++i) {
    //printf(" %lu %19.10e  %19.10e  %19.10e  %19.10e\n", i, ox[i], oy[i], oz[i], wt[i]);

    sum0 += wt[i];
    sum1 += wt[i] * ox[i] * ox[i];
    sum2 += wt[i] * oy[i] * oy[i];
    sum3 += wt[i] * oz[i] * oz[i];
    sum4 += wt[i] * oy[i] * oy[i] * oz[i] * oz[i];
    sum5 += wt[i] * ox[i] * ox[i] * oy[i] * oy[i] * oz[i] * oz[i];
  }

  printf("integral 0 %.10e, should be %.10e\n", sum0, 4 * cPI);
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(sum0, 4 * cPI, tol));

  printf("integral 1 %.10e, should be %.10e\n", sum1, 4 * cPI / 3);
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(sum1, 4 * cPI / 3, tol));

  printf("integral 2 %.10e, should be %.10e\n", sum2, 4 * cPI / 3);
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(sum2, 4 * cPI / 3, tol));

  printf("integral 3 %.10e, should be %.10e\n", sum3, 4 * cPI / 3);
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(sum3, 4 * cPI / 3, tol));

  printf("integral 4 %.10e, should be %.10e\n", sum4, 4 * cPI / 15);
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(sum4, 4 * cPI / 15, tol));

  if (prod_type) // triangular sets can't integrate this funciton exactly
  {
    printf("integral 5 %.10e, should be %.10e\n", sum5, 4 * cPI / 105);
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(sum5, 4 * cPI / 105, tol));
  }
}

void quad_test_2D(int const polar_order, int const azimuthal_order, rtt_dsxx::UnitTest &ut) {
  int quad_type;

  // product quadratures
  quad_type = 0;

  for (int i = 0; i < 4; ++i) // polar types
  {
    for (int j = 0; j < 7; ++j) // azimuthal types
    {
      std::cout << "\nProduct qudrature type (" << i << ", " << j << ")\n";
      auto qp = quad_2D(polar_order, azimuthal_order, i, j, quad_type);
      test_sph(qp, true, ut);
    }
  }

  // triangular quadratures
  quad_type = 1;

  // spence quadratures are accurate only for certain combinations of polynomials of varying
  // degree in sin and cos so they are not tested here

  // a triangular set can be created with az_qrs30 and az_qra30 because the first level
  // requires and N=2 and N must be evenly divisble by 3 for those sets

  const std::vector<int> az_types = {0, 2, 5};

  for (int i = 0; i < 4; ++i) // polar types
  {
    for (auto j : az_types) // azimuthal types
    {
      std::cout << "\nTriangular quadrature type (" << i << ", " << j << ")\n";
      auto qt = quad_2D(polar_order, azimuthal_order, i, j, quad_type);
      test_sph(qt, false, ut);
    }
  }

  // level symmetric quadratures
  quad_type = 2;
  for (int i = 6; i <= 24; ++i, ++i) // start at S6 to integrate exactly
  {
    std::cout << "\nLevel Symmetric order " << i << "\n";
    auto ql = quad_2D(i, 0, 0, 0, quad_type); // level symmetric only uses polar order i
    test_sph(ql, false, ut);
  }

  // test s2 and s4 level symmetric
  for (int j = 2; j <= 4; ++j, ++j) {
    auto ql = quad_2D(2, 0, 0, 0, quad_type);
    std::vector<double> ox, oy, oz, wt;
    ql.create_quadrature(ox, oy, oz, wt);

    double sum0 = 0;
    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;

    for (size_t i = 0; i < wt.size(); ++i) {
      sum0 += wt[i];
      sum1 += wt[i] * ox[i] * ox[i];
      sum2 += wt[i] * oy[i] * oy[i];
      sum3 += wt[i] * oz[i] * oz[i];
    }

    printf("integral 0 %.10e, should be %.10e\n", sum0, 4 * cPI);
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(sum0, 4 * cPI, tol));

    printf("integral 1 %.10e, should be %.10e\n", sum1, 4 * cPI / 3);
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(sum1, 4 * cPI / 3, tol));

    printf("integral 2 %.10e, should be %.10e\n", sum2, 4 * cPI / 3);
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(sum2, 4 * cPI / 3, tol));

    printf("integral 3 %.10e, should be %.10e\n", sum3, 4 * cPI / 3);
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(sum3, 4 * cPI / 3, tol));
  }
}

//-------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);

  ut.passes("lightweight_quadrature");

  try {
    // some partial-range quadratures require high quadrature order to integrate high-degree polynomials accurately

    quad_test_1D(36, ut);

    quad_test_2D(36, 36, ut);
  }

  UT_EPILOG(ut);
}
