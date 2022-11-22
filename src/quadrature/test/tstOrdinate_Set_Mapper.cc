//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/test/tstOrdinate_Set_Mapper.cc
 * \author Allan Wollaber
 * \date   Mon Mar  7 16:21:45 EST 2016
 * \brief  Ordinate Set Mapper tests
 * \note   Copyright (C) 2016-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "quadrature/Gauss_Legendre.hh"
#include "quadrature/Level_Symmetric.hh"
#include "quadrature/Ordinate_Set_Mapper.hh"
#include "quadrature/Product_Chebyshev_Legendre.hh"
#include <algorithm>
#include <numeric>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_quadrature;

//------------------------------------------------------------------------------------------------//
// Unit tests
//------------------------------------------------------------------------------------------------//

// -----------------------------------------------------------------------------
// This is used in an STL algorithm below to count the number of zeros
// -----------------------------------------------------------------------------
bool is_zero(double a) { return soft_equiv(a, 0.0); }

// -----------------------------------------------------------------------------
// A simple functor to be used in computing a bunch of 3D angle dot products between a given
// ordinate and all the ordinates in a container.
// -----------------------------------------------------------------------------
struct dot_product_functor_3D {
  explicit dot_product_functor_3D(const Ordinate &o_in) : o1(o_in) {}

  // Returns the dot product of the ordinate passed into the functor with the local ordinate
  double operator()(const Ordinate &o2) const {
    return o1.mu() * o2.mu() + o1.eta() * o2.eta() + o1.xi() * o2.xi();
  }

  const Ordinate o1;
};

// -----------------------------------------------------------------------------
// A simple functor to be used in computing a bunch of 1D angle dot products between a given
// ordinate and all the ordinates in a container.
// -----------------------------------------------------------------------------
struct dot_product_functor_1D {
  explicit dot_product_functor_1D(const Ordinate &o_in) : o1(o_in) {}

  // Returns the dot product of the ordinate passed into the functor with the local ordinate
  double operator()(const Ordinate &o2) const {
    // This uses a different approach, nominally slower than
    double phi1 = acos(o1.mu());
    double phi2 = acos(o2.mu());

    return cos(phi1 - phi2);
  }

  const Ordinate o1;
};

// -----------------------------------------------------------------------------
// Takes the zeroth angular moment of weights against the ordinate set
// -----------------------------------------------------------------------------
double zeroth_moment(const vector<double> &weights, const Ordinate_Set &os) {
  Require(weights.size() == os.ordinates().size());

  // Vector of all ordinates in the ordinate set
  const vector<Ordinate> &ords(os.ordinates());

  double phi(0.0);
  for (size_t i = 0; i < ords.size(); ++i) {
    phi += ords[i].wt() * weights[i];
  }

  return phi;
}

// -----------------------------------------------------------------------------
// A quick way to test the validity of several ordinates using the nearest neighbor interpolation
// scheme. This test harness does not account for the removal of "starting directions" in the
// Ordinate_Set.
// -----------------------------------------------------------------------------
void nearest_neighbor_test(rtt_dsxx::UnitTest &ut, const Ordinate &ord, const Ordinate_Set &os,
                           const vector<double> &wts, const size_t index) {
  Require(index < wts.size());
  Require(wts.size() == os.ordinates().size());

  const vector<Ordinate> &ordinates(os.ordinates());

  if (wts.size() != ordinates.size())
    ut.failure("Weight/size mismatch");

  // We should get exactly one nonzero entry in the weight vector
  size_t numzeros = count_if(wts.begin(), wts.end(), is_zero);
  if (numzeros != ordinates.size() - 1)
    ut.failure("Found multiple matches");

  // The sum of the weights should be the original ordinate weight
  double E = zeroth_moment(wts, os);
  if (!soft_equiv(E, ord.wt()))
    ut.failure("Weight summation mismatch");

  // We can confirm here that we actually did find the nearest neighbor, too
  vector<double> dps(wts.size(), 0.0);
  if (os.dimension() >= 2) {
    dot_product_functor_3D dpf(ord);
    std::transform(ordinates.begin(), ordinates.end(), dps.begin(), dpf);
  } else {
    dot_product_functor_1D dpf(ord);
    std::transform(ordinates.begin(), ordinates.end(), dps.begin(), dpf);
  }
  // Find the maximum dot product that we just calculated
  size_t max_e = std::max_element(dps.begin(), dps.end()) - dps.begin();

  // Does it correspond to the same location in the weight vector?
  size_t nz_e = std::max_element(wts.begin(), wts.end()) - wts.begin();
  if (max_e != nz_e)
    ut.failure("Nearest ordinate mismatch");

  // Was the weight correctly allocated?
  if (!soft_equiv(wts[nz_e], ord.wt() / ordinates[nz_e].wt()))
    ut.failure("Wrong weight in ordinate");

  // Was the correct location found?
  if (nz_e != index)
    ut.failure("wrong ordinate mapped");
}

// -----------------------------------------------------------------------------
// A quick way to test the validity of several ordinates using the nearest three ordinates
// interpolation scheme. This test harness does not account for the removal of "starting directions"
// in the Ordinate_Set.
// -----------------------------------------------------------------------------
void nearest_three_test(rtt_dsxx::UnitTest &ut, const Ordinate &ord, const Ordinate_Set &os,
                        const vector<double> &wts) {
  const vector<Ordinate> &ordinates(os.ordinates());

  if (wts.size() != ordinates.size())
    ut.failure("Weight/size mismatch");

  // We should get 1 or 3 nonzero entries in the weight vector
  size_t numzeros = count_if(wts.begin(), wts.end(), is_zero);
  if (numzeros != ordinates.size() - 3 && numzeros != ordinates.size() - 1)
    ut.failure("Should have found either 1 or 3 nonzero weights");

  // The sum of the weights should be the original ordinate weight
  double E = zeroth_moment(wts, os);
  if (!soft_equiv(E, ord.wt()))
    ut.failure("Weight summation mismatch");

  // We can confirm here that we actually did find the nearest neighbor, too
  vector<double> dps(wts.size(), 0.0);
  if (os.dimension() >= 2) {
    dot_product_functor_3D dpf(ord);
    std::transform(ordinates.begin(), ordinates.end(), dps.begin(), dpf);
  } else {
    dot_product_functor_1D dpf(ord);
    std::transform(ordinates.begin(), ordinates.end(), dps.begin(), dpf);
  }

  // Find the maximum 3 dot products and their locations
  vector<double> dpcopy(dps);
  std::partial_sort(dps.begin(), dps.begin() + 3, dps.end(), std::greater<>());
  size_t i1 = std::distance(dpcopy.begin(), std::find(dpcopy.begin(), dpcopy.end(), dps[0]));
  size_t i2 = std::distance(dpcopy.begin(), std::find(dpcopy.begin(), dpcopy.end(), dps[1]));
  size_t i3 = std::distance(dpcopy.begin(), std::find(dpcopy.begin(), dpcopy.end(), dps[2]));

  // These 3 indices should provide the ordered weight values in the test weight vector
  if (wts[i1] < wts[i2])
    ut.failure("Max weight is smaller than 2nd weight");
  if (wts[i1] < wts[i3])
    ut.failure("Max weight is smaller than 3rd weight");
  if (dpcopy[i2] > 0.0 && dpcopy[i3] > 0.0) {
    if (wts[i2] < wts[i3])
      ut.failure("Second and third weights are out of order");
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Main functions for testing; this one is 2D and uses nearest neighbor interp.
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void ordinate_set_2D_nn_mapper_test(rtt_dsxx::UnitTest &ut) {
  int N(2); // quadrature order
  rtt_mesh_element::Geometry geometry(rtt_mesh_element::Geometry::CARTESIAN);
  Level_Symmetric quadrature(N);

  std::shared_ptr<Ordinate_Set> os_LS2 =
      quadrature.create_ordinate_set(2, geometry,
                                     1.0,   // norm,
                                     false, //starting_directions?
                                     false, //extra directions?
                                     Ordinate_Set::LEVEL_ORDERED);

  Ordinate_Set_Mapper osm(*os_LS2);
  vector<Ordinate> ordinates(os_LS2->ordinates());

  {
    // Create a (1,1,1) angle in octant one
    double w(2.2);
    Ordinate o1(sqrt(3.0) / 3.0, sqrt(3.0) / 3.0, sqrt(3.0) / 3.0, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_NEIGHBOR, wt_distribute);

    // perform some sanity checks
    nearest_neighbor_test(ut, o1, *os_LS2, wt_distribute, 3);

    if (ut.numFails == 0)
      ut.passes("Level_Symmetric octant 1 map to nearest neighbor passes");
  }

  // Test angle in octant 2
  {
    // Create a (1,-1,1) angle in octant two
    double w(1.8);
    Ordinate o1(sqrt(3.0) / 3.0, -sqrt(3.0) / 3.0, sqrt(3.0) / 3.0, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_NEIGHBOR, wt_distribute);

    // perform some sanity checks
    nearest_neighbor_test(ut, o1, *os_LS2, wt_distribute, 1);

    if (ut.numFails == 0)
      ut.passes("Level_Symmetric octant 2 map to nearest neighbor passes");
  }

  // Test angle in octant 2; this time off the diagonal
  {
    double w(0.1);
    double dp2 = std::sqrt(1.0 + 1.0 + 0.1 * 0.1);
    Ordinate o1(1.0 / dp2, -1.0 / dp2, 0.1 / dp2, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_NEIGHBOR, wt_distribute);

    // perform some sanity checks
    nearest_neighbor_test(ut, o1, *os_LS2, wt_distribute, 1);

    if (ut.numFails == 0)
      ut.passes("Level_Symmetric octant 2 off-diagonal NN passes");
  }

  // Test angle in octant 5; this should fail in 2-D!
  if (ut.dbcOn() && !ut.dbcNothrow()) {
    double w(0.1);
    Ordinate o1(sqrt(3.0) / 3.0, sqrt(3.0) / 3.0, -sqrt(3.0) / 3.0, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    bool caught(false);

    try {
      osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_NEIGHBOR, wt_distribute);
    } catch (const rtt_dsxx::assertion &error) {
      ostringstream message;
      message << "Good, we caught the following exception\n" << error.what();
      PASSMSG(message.str());
      caught = true;
    }
    if (!caught) {
      ostringstream message;
      message << "Failed to catch an exception for out-of-range angle.";
      FAILMSG(message.str());
    }

    if (ut.numFails == 0)
      ut.passes("Caught an octant 5 angle supplied to a 2-D quadrature.");
  }

  if (ut.numFails == 0)
    ut.passes("2-D Level symmetric, nearest-neighbor mapping passes.");

  return;
}

//------------------------------------------------------------------------------------------------//
// Test the mapper for 1-D problems using nearest neighbor interpolation
//------------------------------------------------------------------------------------------------//
void ordinate_set_1D_nn_mapper_test(rtt_dsxx::UnitTest &ut) {
  int N(4); // quadrature order
  rtt_mesh_element::Geometry geometry(rtt_mesh_element::Geometry::CARTESIAN);
  Gauss_Legendre quadrature(N);

  std::shared_ptr<Ordinate_Set> os_LS4 =
      quadrature.create_ordinate_set(1, //1-D
                                     geometry,
                                     1.0,   // norm,
                                     false, //starting_directions?
                                     false, //extra directions?
                                     Ordinate_Set::LEVEL_ORDERED);

  Ordinate_Set_Mapper osm(*os_LS4);
  vector<Ordinate> ordinates(os_LS4->ordinates());

  {
    // Create an angle in octant one
    double w(2.2);
    Ordinate o1(sqrt(3.0) / 3.0, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_NEIGHBOR, wt_distribute);

    // perform some sanity checks
    nearest_neighbor_test(ut, o1, *os_LS4, wt_distribute, 2);

    if (ut.numFails == 0)
      ut.passes("Gauss_Legendre octant 1 map to nearest neighbor passes");
  }

  {
    // Create an angle near -1
    double w(2.2);
    Ordinate o1(-0.9, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_NEIGHBOR, wt_distribute);

    // perform some sanity checks
    nearest_neighbor_test(ut, o1, *os_LS4, wt_distribute, 0);

    if (ut.numFails == 0)
      ut.passes("Gauss_Legendre octant 2 map to nearest neighbor passes");
  }

  {
    // Create an angle near -1
    double w(2.2);
    Ordinate o1(-0.01, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_NEIGHBOR, wt_distribute);

    // perform some sanity checks
    nearest_neighbor_test(ut, o1, *os_LS4, wt_distribute, 1);

    if (ut.numFails == 0)
      ut.passes("Gauss_Legendre octant 2 map to nearest neighbor passes");
  }

  if (ut.dbcOn() && !ut.dbcNothrow()) {
    // Create a (-1,-1,1) angle in octant three
    double w(2.2);
    Ordinate o1(-sqrt(3.0) / 3.0, -sqrt(3.0) / 3.0, sqrt(3.0) / 3.0, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    // perform some sanity checks
    bool caught(false);
    try {
      osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_NEIGHBOR, wt_distribute);
    } catch (const rtt_dsxx::assertion &error) {
      ostringstream message;
      message << "Good, we caught the following exception\n" << error.what();
      PASSMSG(message.str());
      caught = true;
    }
    if (!caught) {
      ostringstream message;
      message << "Failed to catch an exception for out-of-range angle.";
      FAILMSG(message.str());
    }

    if (ut.numFails == 0)
      ut.passes("Caught an octant 3 angle supplied to a 1-D quadrature.");
  }

  if (ut.numFails == 0)
    ut.passes("1-D nearest-neighbor remapping tests all passed");
}

//------------------------------------------------------------------------------------------------//
// Test the mapper for 3D problems with nearest neighbor interpolation
//------------------------------------------------------------------------------------------------//
void ordinate_set_3D_nn_mapper_test(rtt_dsxx::UnitTest &ut) {
  int N(4); // quadrature order
  rtt_mesh_element::Geometry geometry(rtt_mesh_element::Geometry::CARTESIAN);
  Product_Chebyshev_Legendre quadrature(N, N);

  std::shared_ptr<Ordinate_Set> os_LS2 =
      quadrature.create_ordinate_set(3, geometry,
                                     1.0,   // norm,
                                     false, //starting_directions?
                                     false, //extra directions?
                                     Ordinate_Set::LEVEL_ORDERED);

  Ordinate_Set_Mapper osm(*os_LS2);
  vector<Ordinate> ordinates(os_LS2->ordinates());

  {
    // Create an angle in octant 1
    double w(2.2);
    double no1(sqrt(0.9 * 0.9 + 0.8 * 0.8 + 0.1 * 0.1));
    Ordinate o1(0.9 / no1, 0.8 / no1, 0.1 / no1, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_NEIGHBOR, wt_distribute);

    // perform some sanity checks
    nearest_neighbor_test(ut, o1, *os_LS2, wt_distribute, 23);

    if (ut.numFails == 0)
      ut.passes("Product quad octant 1 map to nearest neighbor passes");
  }

  {
    // Create an angle in octant 5
    double w(2.2);
    double no1(sqrt(0.9 * 0.9 + 0.8 * 0.8 + 0.1 * 0.1));
    Ordinate o1(0.9 / no1, 0.8 / no1, -0.1 / no1, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_NEIGHBOR, wt_distribute);

    // perform some sanity checks
    nearest_neighbor_test(ut, o1, *os_LS2, wt_distribute, 22);

    if (ut.numFails == 0)
      ut.passes("Product quad octant 5 map to nearest neighbor passes");
  }

  {
    // Create an angle in octant 6
    double w(2.2);
    double no1(sqrt(0.9 * 0.9 + 0.8 * 0.8 + 0.1 * 0.1));
    Ordinate o1(-0.9 / no1, 0.8 / no1, -0.1 / no1, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_NEIGHBOR, wt_distribute);

    // perform some sanity checks
    nearest_neighbor_test(ut, o1, *os_LS2, wt_distribute, 20);

    if (ut.numFails == 0)
      ut.passes("Product quad octant 6 map to nearest neighbor passes");
  }

  {
    // Create an angle in octant 7
    double w(2.2);
    double no1(sqrt(0.9 * 0.9 + 0.8 * 0.8 + 0.1 * 0.1));
    Ordinate o1(-0.9 / no1, -0.8 / no1, -0.1 / no1, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_NEIGHBOR, wt_distribute);

    // perform some sanity checks
    nearest_neighbor_test(ut, o1, *os_LS2, wt_distribute, 8);

    if (ut.numFails == 0)
      ut.passes("Product quad octant 7 map to nearest neighbor passes");
  }

  {
    // Create an angle in octant 8
    double w(2.2);
    double no1(sqrt(0.9 * 0.9 + 0.8 * 0.8 + 0.1 * 0.1));
    Ordinate o1(0.9 / no1, -0.8 / no1, -0.1 / no1, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_NEIGHBOR, wt_distribute);

    // perform some sanity checks
    nearest_neighbor_test(ut, o1, *os_LS2, wt_distribute, 10);

    if (ut.numFails == 0)
      ut.passes("Product quad octant 8 map to nearest neighbor passes");
  }

  if (ut.numFails == 0)
    ut.passes("3-D nearest-neighbor remapping tests all passed");
}

//------------------------------------------------------------------------------------------------//
// Test the mapper for Spherical problems using nearest neighbor interpolation
//------------------------------------------------------------------------------------------------//
void ordinate_set_1D_sph_nn_mapper_test(rtt_dsxx::UnitTest &ut) {
  int N(2); // quadrature order
  rtt_mesh_element::Geometry geometry(rtt_mesh_element::Geometry::SPHERICAL);
  Gauss_Legendre quadrature(N);

  // Note that this test uses "starting directions" in the Ordinate_Set These are necessary for SN
  // in curvilinear geometries, but these starting directions have zero weight and do not actually
  // contribute to the quadrature integration -- they should not be incluced as valid ordinates
  // during the remapping!
  std::shared_ptr<Ordinate_Set> os_LS2 =
      quadrature.create_ordinate_set(1, //1-D
                                     geometry,
                                     1.0,  // norm,
                                     true, //starting_directions?
                                     true, //extra directions?
                                     Ordinate_Set::LEVEL_ORDERED);

  Ordinate_Set_Mapper osm(*os_LS2);
  vector<Ordinate> ordinates(os_LS2->ordinates());

  {
    // Create an angle in octant one
    double w(2.2);
    Ordinate o1(0.99, w); // This should map to the 3rd ordinate
    // The first and last ordinates are "extra"
    vector<double> wts(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_NEIGHBOR, wts);

    if (wts.size() != ordinates.size())
      ut.failure("Weight/size mismatch");

    // We should get exactly one nonzero entry in the weight vector
    size_t numzeros = count_if(wts.begin(), wts.end(), is_zero);
    if (numzeros != ordinates.size() - 1)
      ut.failure("Found multiple matches");

    // The sum of the weights should be the original ordinate weight
    double E = zeroth_moment(wts, *os_LS2);
    if (!soft_equiv(E, o1.wt()))
      ut.failure("Weight summation mismatch");

    // Does it correspond to the same location in the weight vector?
    size_t nz_e = std::max_element(wts.begin(), wts.end()) - wts.begin();
    if (nz_e != 2)
      ut.failure("Nearest ordinate mismatch");

    if (ut.numFails == 0)
      ut.passes("Gauss_Legendre spherical nearest-neighbor passes");
  }
}

//------------------------------------------------------------------------------------------------//
// Nearest-three ordinate set mapping tests; this one is 1-D
//------------------------------------------------------------------------------------------------//
void ordinate_set_1D_nt_mapper_test(rtt_dsxx::UnitTest &ut) {
  int N(8); // quadrature order
  rtt_mesh_element::Geometry geometry(rtt_mesh_element::Geometry::CARTESIAN);
  Gauss_Legendre quadrature(N);

  std::shared_ptr<Ordinate_Set> os_LS4 =
      quadrature.create_ordinate_set(1, //1-D
                                     geometry,
                                     1.0,   // norm,
                                     false, //starting_directions?
                                     false, //extra directions?
                                     Ordinate_Set::LEVEL_ORDERED);

  Ordinate_Set_Mapper osm(*os_LS4);
  vector<Ordinate> ordinates(os_LS4->ordinates());

  {
    // Create an angle in octant one
    double w(1.0);
    Ordinate o1(1.0, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_THREE, wt_distribute);

    // perform some sanity checks
    nearest_three_test(ut, o1, *os_LS4, wt_distribute);

    if (ut.numFails == 0)
      ut.passes("Gauss_Legendre octant 1 map to nearest three passes");
  }

  {
    // Create an angle near -1
    double w(2.2);
    Ordinate o1(-0.9, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_THREE, wt_distribute);

    // perform some sanity checks
    nearest_three_test(ut, o1, *os_LS4, wt_distribute);

    if (ut.numFails == 0)
      ut.passes("Gauss_Legendre octant 2 map to nearest three passes");
  }

  {
    // Create a grazing angle
    double w(2.2);
    Ordinate o1(-0.01, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_THREE, wt_distribute);

    // perform some sanity checks
    nearest_three_test(ut, o1, *os_LS4, wt_distribute);

    if (ut.numFails == 0)
      ut.passes("Gauss_Legendre octant 2 map to nearest three passes");
  }

  {
    // Create an angle right on top of an ordinate
    double w(2.2);
    Ordinate o1(0.52553240991632899, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_THREE, wt_distribute);

    // perform some sanity checks
    nearest_three_test(ut, o1, *os_LS4, wt_distribute);

    if (ut.numFails == 0)
      ut.passes("Gauss_Legendre identical ord. to nearest three passes");
  }

  {
    // Create an angle verrrrrry close to an ordinate
    double w(2.2);
    Ordinate o1(0.53, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_THREE, wt_distribute);

    // perform some sanity checks
    nearest_three_test(ut, o1, *os_LS4, wt_distribute);

    if (ut.numFails == 0)
      ut.passes("Gauss_Legendre super close map to nearest three passes");
  }

  if (ut.dbcOn() && !ut.dbcNothrow()) {
    // Create a (-1,-1,1) angle in octant three
    double w(2.2);
    Ordinate o1(-sqrt(3.0) / 3.0, -sqrt(3.0) / 3.0, sqrt(3.0) / 3.0, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    // perform some sanity checks
    bool caught(false);
    try {
      osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_NEIGHBOR, wt_distribute);
    } catch (const rtt_dsxx::assertion &error) {
      ostringstream message;
      message << "Good, we caught the following exception\n" << error.what();
      PASSMSG(message.str());
      caught = true;
    }
    if (!caught) {
      ostringstream message;
      message << "Failed to catch an exception for out-of-range angle.";
      FAILMSG(message.str());
    }

    if (ut.numFails == 0)
      ut.passes("Caught an octant 3 angle supplied to a 1-D quadrature.");
  }

  if (ut.numFails == 0)
    ut.passes("1-D nearest-three remapping tests all passed");
}

// -----------------------------------------------------------------------------
// Main functions for testing; this one is 2D and uses nearest three interp.
// -----------------------------------------------------------------------------
void ordinate_set_2D_nt_mapper_test(rtt_dsxx::UnitTest &ut) {
  int N(6); // quadrature order
  rtt_mesh_element::Geometry geometry(rtt_mesh_element::Geometry::CARTESIAN);
  Level_Symmetric quadrature(N);

  std::shared_ptr<Ordinate_Set> os_LS2 =
      quadrature.create_ordinate_set(2, geometry,
                                     1.0,   // norm,
                                     false, //starting_directions?
                                     false, //extra directions?
                                     Ordinate_Set::LEVEL_ORDERED);

  Ordinate_Set_Mapper osm(*os_LS2);
  vector<Ordinate> ordinates(os_LS2->ordinates());

  {
    // Create a (1,1,1) angle in octant one
    double w(2.2);
    Ordinate o1(sqrt(3.0) / 3.0, sqrt(3.0) / 3.0, sqrt(3.0) / 3.0, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_THREE, wt_distribute);

    // perform some sanity checks
    nearest_three_test(ut, o1, *os_LS2, wt_distribute);

    if (ut.numFails == 0)
      ut.passes("Level_Symmetric octant 1 map to nearest three passes");
  }

  // Test angle in octant 2
  {
    // Create a (1,-1,1) angle in octant two
    double w(1.8);
    Ordinate o1(sqrt(3.0) / 3.0, -sqrt(3.0) / 3.0, sqrt(3.0) / 3.0, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_THREE, wt_distribute);

    // perform some sanity checks
    nearest_three_test(ut, o1, *os_LS2, wt_distribute);

    if (ut.numFails == 0)
      ut.passes("Level_Symmetric octant 2 map to nearest three passes");
  }

  // Test angle in octant 2; this time off the diagonal
  {
    double w(0.1);
    double dp2 = std::sqrt(1.0 + 1.0 + 0.1 * 0.1);
    Ordinate o1(1.0 / dp2, -1.0 / dp2, 0.1 / dp2, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_THREE, wt_distribute);

    // perform some sanity checks
    nearest_three_test(ut, o1, *os_LS2, wt_distribute);

    if (ut.numFails == 0)
      ut.passes("Level_Symmetric octant 2 off-diagonal NT passes");
  }

  // Test angle in octant 5; this should fail in 2-D!
  if (ut.dbcOn() && !ut.dbcNothrow()) {
    double w(0.1);
    Ordinate o1(sqrt(3.0) / 3.0, sqrt(3.0) / 3.0, -sqrt(3.0) / 3.0, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    bool caught(false);

    try {
      osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_THREE, wt_distribute);
    } catch (const rtt_dsxx::assertion &error) {
      ostringstream message;
      message << "Good, we caught the following exception\n" << error.what();
      PASSMSG(message.str());
      caught = true;
    }
    if (!caught) {
      ostringstream message;
      message << "Failed to catch an exception for out-of-range angle.";
      FAILMSG(message.str());
    }

    if (ut.numFails == 0)
      ut.passes("Caught an octant 5 angle supplied to a 2-D quadrature.");
  }

  if (ut.numFails == 0)
    ut.passes("2-D Level symmetric, nearest-three mapping passes.");

  return;
}

//------------------------------------------------------------------------------------------------//
// Test the mapper for 3D problems with nearest three interpolation
//------------------------------------------------------------------------------------------------//
void ordinate_set_3D_nt_mapper_test(rtt_dsxx::UnitTest &ut) {
  int N(4); // quadrature order
  rtt_mesh_element::Geometry geometry(rtt_mesh_element::Geometry::CARTESIAN);
  Product_Chebyshev_Legendre quadrature(N, N);

  std::shared_ptr<Ordinate_Set> os_LS2 =
      quadrature.create_ordinate_set(3, geometry,
                                     1.0,   // norm,
                                     false, //starting_directions?
                                     false, //extra directions?
                                     Ordinate_Set::LEVEL_ORDERED);

  Ordinate_Set_Mapper osm(*os_LS2);
  vector<Ordinate> ordinates(os_LS2->ordinates());

  {
    // Create an angle in octant 1
    double w(2.2);
    double no1(sqrt(0.9 * 0.9 + 0.8 * 0.8 + 0.1 * 0.1));
    Ordinate o1(0.9 / no1, 0.8 / no1, 0.1 / no1, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_THREE, wt_distribute);

    // perform some sanity checks
    nearest_three_test(ut, o1, *os_LS2, wt_distribute); // 23

    if (ut.numFails == 0)
      ut.passes("Product quad octant 1 map to nearest three passes");
  }

  {
    // Create an angle in octant 5
    double w(2.2);
    double no1(sqrt(0.9 * 0.9 + 0.8 * 0.8 + 0.1 * 0.1));
    Ordinate o1(0.9 / no1, 0.8 / no1, -0.1 / no1, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_THREE, wt_distribute);

    // perform some sanity checks
    nearest_three_test(ut, o1, *os_LS2, wt_distribute); //22

    if (ut.numFails == 0)
      ut.passes("Product quad octant 5 map to nearest three passes");
  }

  {
    // Create an angle in octant 6
    double w(2.2);
    double no1(sqrt(0.9 * 0.9 + 0.8 * 0.8 + 0.1 * 0.1));
    Ordinate o1(-0.9 / no1, 0.8 / no1, -0.1 / no1, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_THREE, wt_distribute);

    // perform some sanity checks
    nearest_three_test(ut, o1, *os_LS2, wt_distribute); // 20

    if (ut.numFails == 0)
      ut.passes("Product quad octant 6 map to nearest three passes");
  }

  {
    // Create an angle in octant 7
    double w(2.2);
    double no1(sqrt(0.9 * 0.9 + 0.8 * 0.8 + 0.1 * 0.1));
    Ordinate o1(-0.9 / no1, -0.8 / no1, -0.1 / no1, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_THREE, wt_distribute);

    // perform some sanity checks
    nearest_three_test(ut, o1, *os_LS2, wt_distribute); // 8

    if (ut.numFails == 0)
      ut.passes("Product quad octant 7 map to nearest three passes");
  }

  {
    // Create an angle in octant 8
    double w(2.2);
    double no1(sqrt(0.9 * 0.9 + 0.8 * 0.8 + 0.1 * 0.1));
    Ordinate o1(0.9 / no1, -0.8 / no1, -0.1 / no1, w);
    vector<double> wt_distribute(ordinates.size(), 0.0);

    osm.map_angle_into_ordinates(o1, Ordinate_Set_Mapper::NEAREST_THREE, wt_distribute);

    // perform some sanity checks
    nearest_three_test(ut, o1, *os_LS2, wt_distribute); // 10

    if (ut.numFails == 0)
      ut.passes("Product quad octant 8 map to nearest three passes");
  }

  if (ut.numFails == 0)
    ut.passes("3-D nearest-three remapping tests all passed");
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    // Perform the nearest neighbor mapping tests
    ordinate_set_1D_nn_mapper_test(ut);
    ordinate_set_2D_nn_mapper_test(ut);
    ordinate_set_3D_nn_mapper_test(ut);
    ordinate_set_1D_sph_nn_mapper_test(ut);

    // Perform the "nearest three" mapping tests
    ordinate_set_1D_nt_mapper_test(ut);
    ordinate_set_2D_nt_mapper_test(ut);
    ordinate_set_3D_nt_mapper_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstOrdinate_Set_Mapper.cc
//------------------------------------------------------------------------------------------------//
