//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ipcress/test/cdi_ipcress_test.hh
 * \author Thomas M. Evans
 * \date   Fri Oct 12 15:36:36 2001
 * \brief  cdi_ipcress test function headers.
 * \note   Copyright (C) 2011-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_ipcress_test_hh
#define rtt_cdi_ipcress_test_hh

#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include <sstream>

namespace rtt_cdi_ipcress_test {

//------------------------------------------------------------------------------------------------//
// COMPARISON FUNCTIONS USED IN IPCRESS OPACITY TESTS
//------------------------------------------------------------------------------------------------//

template <typename temperatureType, typename densityType, typename testValueType,
          typename opacityClassType>
bool opacityAccessorPassed(rtt_dsxx::ScalarUnitTest &ut, opacityClassType const spOpacity,
                           temperatureType const temperature, densityType const density,
                           testValueType const tabulatedValue) {
  using std::ostringstream;

  // Interpolate the multigroup opacities.
  testValueType grayOpacity = spOpacity->getOpacity(temperature, density);

  // Make sure that the interpolated value matches previous interpolations.

  if (rtt_dsxx::soft_equiv(grayOpacity, tabulatedValue)) {
    ostringstream message;
    message << spOpacity->getDataDescriptor() << " opacity computation was good for \n\t"
            << "\"" << spOpacity->getDataFilename() << "\" data.";
    ut.passes(message.str());
  } else {
    ostringstream message;
    message << spOpacity->getDataDescriptor() << " opacity value is out of spec. for \n\t"
            << "\"" << spOpacity->getDataFilename() << "\" data.";
    return ut.failure(message.str());
  }

  // If we get here then the test passed.
  return true;
}

//------------------------------------------------------------------------------------------------//

template <typename opacityClassType>
void testTemperatureGridAccessor(rtt_dsxx::ScalarUnitTest &ut, opacityClassType const spOpacity) {
  using std::ostringstream;

  // Read the temperature grid from the IPCRESS file.
  std::vector<double> temps = spOpacity->getTemperatureGrid();

  // Verify that the size of the temperature grid looks right.  If it is the right size then compare
  // the temperature grid data to the data specified when we created the IPCRESS file using TOPS.
  if (temps.size() == spOpacity->getNumTemperatures() && temps.size() == 3) {
    ostringstream message;
    message << "The number of temperature points found in the data\n\t"
            << "grid matches the number returned by the\n\t"
            << "getNumTemperatures() accessor.";
    ut.passes(message.str());

    // The grid specified by TOPS has 3 temperature points.
    std::vector<double> temps_ref = {0.1, 1.0, 10.0};
    Check(temps.size() == 3);

    // Compare the grids.
    if (rtt_dsxx::soft_equiv(temps, temps_ref))
      ut.passes("Temperature grid matches.");
    else
      ut.failure("Temperature grid did not match.");
  } else {
    ostringstream message;
    message << "The number of temperature points found in the data\n\t"
            << "grid does not match the number returned by the\n\t"
            << "getNumTemperatures() accessor. \n"
            << "Did not test the results returned by\n\t"
            << "getTemperatureGrid().";
    ut.failure(message.str());
  }
}

//------------------------------------------------------------------------------------------------//

template <typename opacityClassType>
void testDensityGridAccessor(rtt_dsxx::ScalarUnitTest &ut, opacityClassType const spOpacity) {
  using std::ostringstream;

  // Read the grid from the IPCRESS file.
  std::vector<double> density = spOpacity->getDensityGrid();

  // Verify that the size of the density grid looks right.  If it is the right
  // size then compare the density grid data to the data specified when we
  // created the IPCRESS file using TOPS.
  if (density.size() == 3 && density.size() == spOpacity->getNumDensities()) {
    ostringstream message;
    message << "The number of density points found in the data\n\t"
            << "grid matches the number returned by the\n\tgetNumDensities() accessor.";
    ut.passes(message.str());

    // The grid specified by TOPS has 3 density points
    std::vector<double> density_ref = {0.1, 0.5, 1.0};
    Check(density.size() == 3);

    // Compare the grids.
    if (rtt_dsxx::soft_equiv(density, density_ref))
      ut.passes("Density grid matches.");
    else
      ut.failure("Density grid did not match.");
  } else {
    ostringstream message;
    message << "The number of density points found in the data\n\t"
            << "grid does not match the number returned by the\n\t"
            << "getNumDensities() accessor. \n"
            << "Did not test the results returned by\n\tgetDensityGrid().";
    ut.failure(message.str());
  }
}

//------------------------------------------------------------------------------------------------//
template <typename opacityClassType>
void testEnergyBoundaryAccessor(rtt_dsxx::ScalarUnitTest &ut, opacityClassType const spOpacity) {
  using std::ostringstream;

  // Read the grid from the IPCRESS file.
  std::vector<double> ebounds = spOpacity->getGroupBoundaries();

  // Verify that the size of the group boundary grid looks right.  If it is the right size then
  // compare the energy groups grid data to the data specified when we created the IPCRESS file
  // using TOPS.
  if (ebounds.size() == 13 && ebounds.size() == spOpacity->getNumGroupBoundaries()) {
    ostringstream message;
    message << "The number of energy boundary points found in the data\n\t"
            << "grid matches the number returned by the\n\t"
            << "getNumGroupBoundaries() accessor.";
    ut.passes(message.str());

    // The grid specified by TOPS has 13 energy boundaries.
    std::vector<double> ebounds_ref = {0.01, 0.03, 0.07, 0.1,  0.3,  0.7,  1.0,
                                       3.0,  7.0,  10.0, 30.0, 70.0, 100.0};
    Check(ebounds.size() == 13);

    // Compare the grids.
    if (rtt_dsxx::soft_equiv(ebounds, ebounds_ref))
      ut.passes("Energy group boundary grid matches.");
    else
      ut.failure("Energy group boundary grid did not match.");
  } else {
    ostringstream message;
    message << "The number of energy boundary points found in the data\n\t"
            << "grid does not match the number returned by the\n\t"
            << "get NumGroupBoundaries() accessor. \n"
            << "Did not test the results returned by\n\tgetGroupBoundaries().";
    ut.failure(message.str());
  }
}

} // end namespace rtt_cdi_ipcress_test

#endif // rtt_cdi_ipcress_test_hh

//------------------------------------------------------------------------------------------------//
// end of cdi_ipcress/test/cdi_ipcress_test.hh
//------------------------------------------------------------------------------------------------//
