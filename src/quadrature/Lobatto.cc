//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Lobatto.cc
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  A class representing an interval Gauss-Legendre quadrature set.
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Lobatto.hh"
#include "gauleg.hh"
#include "ds++/DracoStrings.hh"
#include <numeric>

namespace rtt_quadrature {
using namespace std;
using rtt_dsxx::to_string;

//------------------------------------------------------------------------------------------------//
Lobatto::Lobatto(unsigned sn_order) : Interval_Quadrature(sn_order) {
  Require(sn_order > 0 && sn_order % 2 == 0);

  // base class data members
  mu_.resize(sn_order);
  wt_.resize(sn_order);

  // Mus and weights are worth precomputing

  using rtt_dsxx::soft_equiv;
  using rtt_units::PI;
  using std::cos;
  using std::numeric_limits;

  double const tolerance(100 * std::numeric_limits<double>::epsilon());

  // size the member data vectors
  unsigned const numGaussPoints = sn_order_;
  unsigned const N(numGaussPoints);

  double const mu1(-1.0); // minimum value for mu
  double const mu2(1.0);  // maximum value for mu

  // Number of Gauss points in the half range.
  // The roots are symmetric in the interval.  We only need to search for
  // half of them.
  unsigned const numHrGaussPoints(numGaussPoints / 2);

  mu_[0] = mu1;
  mu_[N - 1] = mu2;

  // Loop over the desired roots.
  for (size_t iroot = 0; iroot < numHrGaussPoints - 1; ++iroot) {
    // Approximate the i-th root.
    double z(cos(PI * (static_cast<double>(iroot) - 0.25) / (static_cast<double>(N - 2) + 0.5)));
    double z1;

    do // Use Newton's method to refine the value for the i-th root.
    {
      // P_{N-1}(z)
      double const pnm1(gsl_sf_legendre_Pl(N - 1, z));
      // P_{N}(z)
      double const pn(gsl_sf_legendre_Pl(N, z));
      // P_{N+1}(z)
      double const pnp1(gsl_sf_legendre_Pl(N + 1, z));

      // dP/dz _{N-1}(z)
      double const pp((N) * (z * pnm1 - pn) / (1.0 - z * z));
      // dP/dz _{N}(z)
      double const pp1((N + 1) * (z * pn - pnp1) / (1.0 - z * z));

      // d2P/dz2 _{N}(z)
      double const pdp(((N) * (z * pp + pnm1 - pp1) + 2 * z * pp) / (1.0 - z * z));

      // Do synthetic division to avoid roots already found

      double Dpp = pp;
      double Dpdp = pdp;
      for (unsigned i = 1; i <= iroot; ++i) {
        Dpp /= z + mu_[i]; // Use positive root
        Dpdp /= z + mu_[i];
      }
      for (unsigned i = 1; i <= iroot; ++i) {
        Dpdp -= Dpp / (z + mu_[i]);
      }

      // update

      z1 = z;
      z = z1 - Dpp / Dpdp;

    } while (!soft_equiv(z, z1, tolerance));

    // Roots will be in [-1,1], symmetric about the origin.
    mu_[iroot + 1] = -z;
    mu_[numGaussPoints - iroot - 2] = z;
  }

  // Loop over the quadrature points to compute weights.
  for (size_t m = 0; m < numHrGaussPoints; ++m) {
    double const z(mu_[m]);
    double const p(gsl_sf_legendre_Pl(N - 1, z));

    // Compute the associated weight and its symmetric counterpart.
    wt_[m] = 2.0 / N / (N - 1) / p / p;
    wt_[numGaussPoints - m - 1] = wt_[m];
  }

  Ensure(check_class_invariants());
  Ensure(this->sn_order() == sn_order);
}

//------------------------------------------------------------------------------------------------//
/* virtual */
string Lobatto::name() const { return "Lobatto"; }

//------------------------------------------------------------------------------------------------//
/* virtual */
string Lobatto::parse_name() const { return "axial"; }

//------------------------------------------------------------------------------------------------//
/* virtual */
unsigned Lobatto::number_of_levels() const { return sn_order_; }

//------------------------------------------------------------------------------------------------//
/* virtual */ string Lobatto::as_text(string const &indent) const {
  string Result =
      indent + "type = lobatto" + indent + "  order = " + to_string(sn_order()) + indent + "end";

  return Result;
}

//------------------------------------------------------------------------------------------------//
bool Lobatto::check_class_invariants() const { return sn_order() > 0 && sn_order() % 2 == 0; }

//------------------------------------------------------------------------------------------------//
/* virtual */
vector<Ordinate> Lobatto::create_level_ordinates_(double const norm) const {
  // Sanity Checks: none at present

  unsigned const numPoints(sn_order());

  double sumwt = 0.0;
  for (size_t i = 0; i < numPoints; ++i)
    sumwt += wt_[i];

  double c = norm / sumwt;

  // build the set of ordinates
  vector<Ordinate> Result(numPoints);
  for (size_t i = 0; i < numPoints; ++i) {
    // This is a 1D set.
    Result[i] = Ordinate(mu_[i], c * wt_[i]);
  }

  return Result;
}

//------------------------------------------------------------------------------------------------//
bool Lobatto::is_open_interval() const {
  // Lobatto is one of our few closed interval quadratures.
  return false;
}

} // end namespace rtt_quadrature

//------------------------------------------------------------------------------------------------//
// end of quadrature/Lobatto.cc
//------------------------------------------------------------------------------------------------//
