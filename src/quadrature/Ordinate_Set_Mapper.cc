//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Ordinate_Set_Mapper.cc
 * \author Allan Wollaber
 * \date   Mon Mar  7 10:42:56 EST 2016
 * \brief  Implementation file for the class
 *         rtt_quadrature::Ordinate_Set_Mapper.
 * \note   Copyright (C) 2016-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Ordinate_Set_Mapper.hh"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>

namespace {

using namespace std;
using namespace rtt_quadrature;

// convenience functions to check ordinates

#if DBC & 1

//------------------------------------------------------------------------------------------------//
bool check_4(Ordinate const &ordinate) {
  // In 1-D spherical geometry, the ordinates must be confined to the first
  // two octants.
  if (ordinate.eta() < 0 || ordinate.xi() < 0)
    return false;
  return true;
}

//------------------------------------------------------------------------------------------------//
bool check_2(Ordinate const &ordinate) {
  // In 2-D geometry, the ordinates must be confined to the first
  // four octants
  if (ordinate.xi() < 0)
    return false;
  return true;
}

#endif

//------------------------------------------------------------------------------------------------//
using dsp = std::pair<double, size_t>;
bool bigger_pair(const dsp &d1, const dsp &d2) { return (d1.first > d2.first); }

} // end anonymous namespace

namespace rtt_quadrature {

//------------------------------------------------------------------------------------------------//
bool Ordinate_Set_Mapper::check_class_invariants() const { return os_.check_class_invariants(); }

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Primary service method to perform weight mapping
 *
 * \param[in] ord_in Ordinate with associated weight to remap to Ordinate Set
 * \param[in] interp_in Selected interpolation scheme to use for remapping
 * \param[out] weights vector of output weights produced by remapping
 *
 *  The output weights preserve the incoming weight assuming that the incoming
 *  ordinate's weight is associated via a delta function in angle.
 */
void Ordinate_Set_Mapper::map_angle_into_ordinates(const Ordinate &ord_in,
                                                   const Interpolation_Type &interp_in,
                                                   std::vector<double> &weights) const {

  Require(os_.ordinates().size() == weights.size());
  Require(os_.dimension() == 2 ? check_2(ord_in) : true);
  Require(os_.dimension() == 1 ? check_4(ord_in) : true);
  // check norm == 1 in 2-D and 3-D
  Require(os_.dimension() >= 2 ? rtt_dsxx::soft_equiv(dot_product_functor_3D(ord_in)(ord_in), 1.0)
                               : true);
  // check norm <= 1 in 1-D
  Require(os_.dimension() == 1 ? dot_product_functor_1D(ord_in)(ord_in) <= 1.0 : true);

  // Vector of all ordinates in the ordinate set
  const vector<Ordinate> &ords(os_.ordinates());

  // Vector of dot products
  vector<double> dps(weights.size(), 0.0);
  // Perform all of the dot products between the incoming ordinate and the
  // ordinates in the Ordinate_Set
  if (os_.dimension() != 1) {
    dot_product_functor_3D dpf(ord_in);
    std::transform(ords.begin(), ords.end(), dps.begin(), dpf);
  } else {
    dot_product_functor_1D dpf(ord_in);
    std::transform(ords.begin(), ords.end(), dps.begin(), dpf);
  }

  // Remove the "starting directions" as valid ordinates in the quadrature, if
  // they exist
  if (os_.has_starting_directions()) {
    size_t i = 0;
    for (auto ord = ords.begin(); ord != ords.end(); ++ord, ++i) {
      // If the ordinate weight is zero, we found a "starting direction". We
      // remove it by saying its dot product is completely opposed to the
      // ordinate we passed in (set to -1).
      if (ord->wt() <= 0.0)
        dps[i] = -1.0;
    }
  }

  switch (interp_in) {
  case NEAREST_NEIGHBOR: {
    // Find the index with the largest dot product
    size_t max_e = std::max_element(dps.begin(), dps.end()) - dps.begin();

    // Put all of the associated weight into that element
    weights[max_e] = ord_in.wt() / ords[max_e].wt();
  } break;

  case NEAREST_THREE: {
    Require(dps.size() >= 3);

    // Vector of all ordinates in the ordinate set
    const std::vector<Ordinate> &ord(os_.ordinates());

    // Associate a container of indices with the dot products
    std::vector<std::pair<double, size_t>> dpsi(dps.size());
    size_t i(0);
    for (auto elem = dpsi.begin(); elem != dpsi.end(); ++elem, ++i) {
      elem->first = dps[i];
      elem->second = i;
    }
    // This sorts the dot products to find the largest 3 of them
    std::nth_element(dpsi.begin(), dpsi.begin() + 3, dpsi.end(), bigger_pair);

    // Assign weights based on normalization of nearest 3
    double w1(dpsi[0].first), w2(dpsi[1].first), w3(dpsi[2].first);
    size_t i1(dpsi[0].second), i2(dpsi[1].second), i3(dpsi[2].second);

    // Prevent adding energy into negative dot-product ordinates
    w1 = std::max(w1, 0.0);
    w2 = std::max(w2, 0.0);
    w3 = std::max(w3, 0.0);
    Check(w1 <= 1.0);
    Check(w2 <= 1.0);
    Check(w3 <= 1.0);

    // Allocate the weights using the inverse distance
    const double ord_tol(1e-6); //
    w1 = 1.0 - w1;
    w2 = 1.0 - w2;
    w3 = 1.0 - w3;

    // This block selects a single ordinate in the ordinate set if one of the
    // dot products is very near to 1.0 (to within ord_tol)
    if (w1 > ord_tol && w2 > ord_tol && w3 > ord_tol) {
      w1 = 1.0 / w1;
      w2 = 1.0 / w2;
      w3 = 1.0 / w3;
    } else if (w1 < ord_tol) {
      w1 = 1.0;
      w2 = 0.0;
      w3 = 0.0;
    } else if (w2 < ord_tol) {
      w1 = 0.0;
      w2 = 1.0;
      w3 = 0.0;
    } else if (w3 < ord_tol) {
      w1 = 0.0;
      w2 = 0.0;
      w3 = 1.0;
    }

    double wsum(w1 + w2 + w3);
    Check(wsum > 0.0);

    // Normalize the 3 weights
    w1 = w1 * ord_in.wt() / (ord[i1].wt() * wsum);
    w2 = w2 * ord_in.wt() / (ord[i2].wt() * wsum);
    w3 = w3 * ord_in.wt() / (ord[i3].wt() * wsum);

    weights[i1] = w1;
    weights[i2] = w2;
    weights[i3] = w3;
  } break;

  default:
    Insist(false, "Unimplemented interpolation type");
    break;
  }

  // Test for energy conservation by integrating over all angles i.e., summing
  // the quadrature weights
  Ensure(rtt_dsxx::soft_equiv(zeroth_moment(weights), ord_in.wt()));
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Simple private function to integrate the zeroth angular moment
 *
 * \param[in] weights vector in the ordinates
 *
 * \return double value for the zeroth moment
 */
double Ordinate_Set_Mapper::zeroth_moment(const std::vector<double> &weights) const {
  Require(weights.size() == os_.ordinates().size());

  // Vector of all ordinates in the ordinate set
  const vector<Ordinate> &ords(os_.ordinates());

  double phi(0.0);
  for (size_t i = 0; i < ords.size(); ++i) {
    phi += ords[i].wt() * weights[i];
  }

  return phi;
}

} // end namespace rtt_quadrature

//------------------------------------------------------------------------------------------------//
// end of quadrature/Ordinate_Set_Mapper.cc
//------------------------------------------------------------------------------------------------//
