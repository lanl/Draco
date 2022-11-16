//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/Interpolate.cc
 * \author Kendra P. Long, Kent G. Budge
 * \date   Wed Jan 22 15:18:23 MST 2003
 * \brief  Interpolation functions.
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Interpolate.hh"

namespace rtt_dsxx {
namespace interpolate {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do a 3D linear interpolation between vertices of a rectangular prism.
 *
 * Algorithm from wikipedia's Trilinear Interpolation article, hat tip to E.  Norris for the
 * reference.
 *
 * \param[in] x0   lower x coordinate of lattice
 * \param[in] x1   upper x coordinate of lattice
 * \param[in] y0   lower y coordinate of lattice
 * \param[in] y1   upper y coordinate of lattice
 * \param[in] z0   lower z coordinate of lattice
 * \param[in] z1   upper z coordinate of lattice
 * \param[in] f000 function at (x0,y0,z0)
 * \param[in] f100 function at (x1,y0,z0)
 * \param[in] f001 function at (x0,y0,z1)
 * \param[in] f101 function at (x1,y0,z1)
 * \param[in] f010 function at (x0,y1,z0)
 * \param[in] f110 function at (x1,y1,z0)
 * \param[in] f011 function at (x0,y1,z1)
 * \param[in] f111 function at (x1,y1,z1)
 * \param[in] x    x coordinate of interpolation point
 * \param[in] y    y coordinate of interpolation point
 * \param[in] z    z coordinate of interpolation point
 * \return The function value linearly interpolated to (x,y,z)
 */
double linear_3d(double const x0, double const x1, double const y0, double const y1,
                 double const z0, double const z1, double const f000, double const f100,
                 double const f001, double const f101, double const f010, double const f110,
                 double const f011, double const f111, double const x, double const y,
                 double const z) {
  Require(std::abs(x1 - x0) > std::numeric_limits<double>::epsilon());
  Require(std::abs(y1 - y0) > std::numeric_limits<double>::epsilon());
  Require(std::abs(z1 - z0) > std::numeric_limits<double>::epsilon());
  Require(x >= x0);
  Require(x <= x1);
  Require(y >= y0);
  Require(y <= y1);
  Require(z >= z0);
  Require(z <= z1);

  const double xd = (x - x0) / (x1 - x0);
  const double yd = (y - y0) / (y1 - y0);
  const double zd = (z - z0) / (z1 - z0);

  const double f00 = f000 * (1. - xd) + f100 * xd;
  const double f01 = f001 * (1. - xd) + f101 * xd;
  const double f10 = f010 * (1. - xd) + f110 * xd;
  const double f11 = f011 * (1. - xd) + f111 * xd;

  const double f0 = f00 * (1. - yd) + f10 * yd;
  const double f1 = f01 * (1. - yd) + f11 * yd;

  const double f = f0 * (1. - zd) + f1 * zd;

  return f;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute the Lagrange multipliers for a set of breakpoint regions,
 *        each with n_local points. Return the multipliers to to the calling function.
 *
 * \param[in] n_break The global number of breakpoints (each pair defines a sub-region)
 * \param[in] n_local The number of points per (sub-)region
 * \param[in] points A set of evaluation points
 * \return A vector of Lagrange multipliers for the given evaluation points
 */
std::vector<double> lagrange_multipliers(const size_t n_break, const size_t n_local,
                                         const std::vector<double> &points) {
  // make sure the data is the right size:
  Insist(points.size() == (n_break - 1) * n_local,
         "Points passed to interpolate::lagrange_multipliers are wrong size!");

  std::vector<double> mults(points.size(), 0.0);

  for (size_t i = 0; i < n_break - 1; ++i) {
    for (size_t j = 0; j < n_local; ++j) {
      const size_t k{j + n_local * i};
      mults[k] = 1.;
      for (size_t l = 0; l < n_local; ++l) {
        if (l != j) {
          const size_t k1{l + n_local * i};
          mults[k] = mults[k] / (points[k] - points[k1]);
        }
      }
    }
  }
  return mults;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Perform one-dimensional polynomial interpolation, given y values
 *         (data), x values (xs), multipliers (cxs), and the value to interpolate to (x)
 *
 * \param[in] data y values
 * \param[in] xs   x data points
 * \param[in] cxs  Lagrange multipliers for data points
 * \param[in] x    value to interpolate to.
 * \return The y value associated with x based on Lagrange polynomial interpolation).
 *
 */
double lagrange_1d(const std::vector<double> &data, const std::vector<double> &xs,
                   const std::vector<double> &cxs, const double x) {
  // make sure all the data is the same size:
  Insist(data.size() == xs.size(), "Sizes passed to interpolate::lagrange_1d do not match!");
  Insist(data.size() == cxs.size(), "Sizes passed to interpolate::lagrange_1d do not match!");

  const size_t nx_local{cxs.size()};

  // Here, we check to see if the data point lies directly on one of our
  // eval points. If it does, we just return the value at the data point.
  double phi = 1.0;
  for (size_t j = 0; j < nx_local; ++j) {
    // use the relative diff (scaled by x) if x is "big enough", or absolute diff if x is O(eps):
    double rel_dx =
        fabs(x) > std::numeric_limits<double>::epsilon() ? (fabs(x - xs[j]) / x) : fabs(x - xs[j]);
    if (rel_dx < 1.0e-6) {
      return data[j];
    } else {
      phi *= (x - xs[j]);
    }
  }
  double val = 0.;
  for (size_t j = 0; j < nx_local; ++j) {
    val += data[j] * (phi * cxs[j]) / (x - xs[j]);
  }

  return val;
}

} // namespace interpolate
} // namespace rtt_dsxx
