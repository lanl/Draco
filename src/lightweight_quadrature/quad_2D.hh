//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   lightweight_quadrature/quad_2D.hh
 * \brief  two-dimensional quadratures for integration over the unit sphere combining various 
 * polar and azimuthal quadratures into product or triangular arrangements, as well as good ol'
 * level-symmetric quadrature.
 * \note   Copyright (C) 2025-2026 Triad National Security, LLC., All rights reserved.
 */
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

#ifndef QUAD_2D_HH
#define QUAD_2D_HH

#include <cassert>
#include <functional>

#include "quad_1D.hh"
#include "quad_level_symmetric.hh"

// set number of directions (need to check this against host code at some point)
// n is polar order, m is aziumuthal order
// geom == 0 is Cartesian, geom == 1 is axisymmetric

inline auto qdirs = [](int geom, int dim, int n, int m, int qtype) {
  int ndirs(0);
  switch (qtype) {
  case 0:
    if (geom == 0)
      ndirs = 2 * n * m / (dim == 3 ? 1 : 2);
    else {
      if (dim == 3)
        throw; // LCOV_EXCL_LINE
      ndirs = n * (m + 1);
    }
    break;
  case 1:
  case 2:
    if (geom == 0)
      ndirs = (n + 2) * n / (dim == 3 ? 1 : 2);
    else {
      if (dim == 3)
        throw; // LCOV_EXCL_LINE
      ndirs = n * n / 2 + (2 * n);
    }
    break;

  default:
    throw; // LCOV_EXCL_LINE
  }
  return ndirs;
};

namespace rtt_lightweight_quadrature {

class quad_2D {

public:
  quad_2D(const int polar_order, const int azimuthal_order, const int polar_type,
          const int azimuthal_type, const int quad_type);

private:
  int const polar_order_;
  int const azimuthal_order_;
  int const polar_type_;
  int const azimuthal_type_;
  int const quad_type_;

  // product arrangement
  template <typename P, typename A>
  auto product(P &polar, A &azimuthal, std::vector<double> &omega_x, std::vector<double> &omega_y,
               std::vector<double> &omega_z, std::vector<double> &weights) const;

  // triangular arrangement
  template <typename P, typename A>
  auto triangular(P &polar, A &azimuthal, std::vector<double> &omega_x,
                  std::vector<double> &omega_y, std::vector<double> &omega_z,
                  std::vector<double> &weights) const;

public:
  // construct quadrature set
  void create_quadrature(std::vector<double> &omega_x, std::vector<double> &omega_y,
                         std::vector<double> &omega_z, std::vector<double> &weights) const;
};

} // namespace rtt_lightweight_quadrature

#endif // QUAD_2D_HH
