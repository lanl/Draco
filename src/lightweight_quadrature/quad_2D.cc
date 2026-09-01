//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   lightweight_quadrature/quad_2D.cc
 * \brief  two-dimensional quadratures for integration over the unit sphere combining various polar
 *         and azimuthal quadratures into product or triangular arrangements, as well as good ol'
 *         level-symmetric quadrature.
 * \note   Copyright (C) 2025-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
//      1D polar           1D azimuthal
//------------------------------------------------------------------------------------------------//
//     gauss                chebyshev
//     double-gauss
//     lobatto              az_qrs45
//     radau                az_qra45
//     double-radau         az_qrs30
//                          az_qra30
//     pol_qrs90            az_qrs90
//     pol_qrs45            az_qra90
//
//------------------------------------------------------------------------------------------------//
//     recommended product quadratures (without closed 1D polar quadratures)
//      polar              azimuthal
//   (0) gauss           (0) chebyshev
//   (1) double-gauss
//
//   (2) qrs90           (1) az_qrs45
//   (3) qrs45           (2) az_qra45
//                       (3) az_qrs30
//                       (4) az_qra30
//                       (5) az_qrs90
//                       (6) az_qra90
//
//------------------------------------------------------------------------------------------------//
//   recommended triangular quadratures (based on integration accuracy)
//      polar              azimuthal
//   (0) gauss           (0) chebyshev
//   (1) double-gauss
//
//   (2) qrs90           (2) az_qra45
//   (3) qrs45           (5) az_qrs90
//
//------------------------------------------------------------------------------------------------//

#include <cassert>
#include <cmath>
#include <cstdio>
#include <functional>
#include <vector>

#include "quad_1D.hh"
#include "quad_2D.hh"
#include "quad_level_symmetric.hh"

namespace rtt_lightweight_quadrature {

quad_2D::quad_2D(const int polar_order, const int azimuthal_order, const int polar_type,
                 const int azimuthal_type, const int quad_type)
    : polar_order_(polar_order), azimuthal_order_(azimuthal_order), polar_type_(polar_type),
      azimuthal_type_(azimuthal_type), quad_type_(quad_type) {
  assert(polar_order >= 2);
  assert((quad_type == 0 && azimuthal_order >= 2) || quad_type > 0);
}

template <typename P, typename A>
auto quad_2D::product(P &polar, A &azimuthal, std::vector<double> &omega_x,
                      std::vector<double> &omega_y, std::vector<double> &omega_z,
                      std::vector<double> &weights) const {
  // polar quadrature
  auto pol = polar(polar_order_);
  std::vector<double> const &z = pol[0];
  std::vector<double> const &w_z = pol[1];

  // azimuthal quadrature
  auto az = azimuthal(azimuthal_order_);
  std::vector<double> const &x = az[0];
  std::vector<double> const &w_x = az[1];

  int const num = 2 * polar_order_ * azimuthal_order_;
  omega_x.reserve(num);
  omega_y.reserve(num);
  omega_z.reserve(num);
  weights.reserve(num);

  // product quadrature
  for (int m = 0; m < polar_order_; m++) {
    for (int n = 0; n < azimuthal_order_; n++) {
      const double mu = std::sqrt(1 - z[m] * z[m]) * x[n];
      const double eta = std::sqrt(1 - z[m] * z[m]) * std::sqrt(1 - x[n] * x[n]);
      const double xi = z[m];
      const double w = w_z[m] * w_x[n];

      omega_x.push_back(mu);
      omega_y.push_back(eta);
      omega_z.push_back(xi);
      weights.push_back(w);

      omega_x.push_back(mu);
      omega_y.push_back(-eta);
      omega_z.push_back(xi);
      weights.push_back(w);
    }
  }
}

// triangular arrangement
template <typename P, typename A>
auto quad_2D::triangular(P &polar, A &azimuthal, std::vector<double> &omega_x,
                         std::vector<double> &omega_y, std::vector<double> &omega_z,
                         std::vector<double> &weights) const {
  // polar quadrature
  auto pol = polar(polar_order_);
  std::vector<double> const &z = pol[0];
  std::vector<double> const &w_z = pol[1];

  int const num = (polar_order_ + 2) * polar_order_;

  omega_x.reserve(num);
  omega_y.reserve(num);
  omega_z.reserve(num);
  weights.reserve(num);

  // triangular quadrature
  for (int m = 0; m < polar_order_ / 2; m++) {
    // azimuthal quadrature
    auto az = azimuthal(2 * (m + 1));
    std::vector<double> const &x = az[0];
    std::vector<double> const &w_x = az[1];

    for (int n = 0; n < 2 * (m + 1); n++) {
      const double mu = std::sqrt(1 - z[m] * z[m]) * x[n];
      const double eta = std::sqrt(1 - z[m] * z[m]) * std::sqrt(1 - x[n] * x[n]);
      const double xi = z[m];
      const double w = w_z[m] * w_x[n];

      omega_x.push_back(mu);
      omega_y.push_back(eta);
      omega_z.push_back(xi);
      weights.push_back(w);

      omega_x.push_back(mu);
      omega_y.push_back(-eta);
      omega_z.push_back(xi);
      weights.push_back(w);

      omega_x.push_back(mu);
      omega_y.push_back(eta);
      omega_z.push_back(-xi);
      weights.push_back(w);

      omega_x.push_back(mu);
      omega_y.push_back(-eta);
      omega_z.push_back(-xi);
      weights.push_back(w);
    }
  }
}

// create quadrature set
void quad_2D::create_quadrature(std::vector<double> &omega_x, std::vector<double> &omega_y,
                                std::vector<double> &omega_z, std::vector<double> &weights) const {
  // polar quadrature
  std::function<std::vector<std::vector<double>>(int)> polar;

  switch (polar_type_) {
  // gauss-legendre
  case 0:
    polar = gauss;
    break;
  case 1:
    polar = doublegauss;
    break;
  case 2:
    polar = pol_qrs90;
    break;
  case 3:
    polar = pol_qrs45;
    break;
  default: {
    printf("polar type not found\n"); // LCOV_EXCL_LINE
    throw;                            // LCOV_EXCL_LINE
  }
  };

  // azimuthal quadrature
  std::function<std::vector<std::vector<double>>(int)> azimuthal;
  switch (azimuthal_type_) {
  // chebyshev
  case 0:
    azimuthal = chebyshev;
    break;
    // Spence
  case 1:
    azimuthal = az_qrs90;
    break;
  case 2:
    azimuthal = az_qrs45;
    break;
  case 3:
    azimuthal = az_qrs30;
    break;
    // Abu-Shumays
  case 4:
    azimuthal = az2D_qra90;
    break;
  case 5:
    azimuthal = az_qra45;
    break;
  case 6:
    azimuthal = az_qra30;
    break;
  default: {
    printf("azimuthal type not found\n"); // LCOV_EXCL_LINE
    throw;                                // LCOV_EXCL_LINE
  }
  };

  // type
  switch (quad_type_) {
  // product quadrature
  case 0:
    /*
      printf("creating product quadrature order (%d, %d) type (%d, %d) \n",
      polar_order_,
      azimuthal_order_,
      polar_type_,
      azimuthal_type_);
     */
    product(polar, azimuthal, omega_x, omega_y, omega_z, weights);
    break;
    // triangular quadrature
  case 1:
    /*
      printf("creating triangular quadrature order (%d) type (%d, %d) \n",
      polar_order_,
      polar_type_,
      azimuthal_type_);
     */
    triangular(polar, azimuthal, omega_x, omega_y, omega_z, weights);
    break;
    // level symmetric quadrature
  case 2:
    // printf("creating level symmetric quadrature order (%d) \n", polar_order_);
    level_symmetric(polar_order_, omega_x, omega_y, omega_z, weights);
    break;
  default:
    throw; // LCOV_EXCL_LINE
  };
}

} // namespace rtt_lightweight_quadrature
