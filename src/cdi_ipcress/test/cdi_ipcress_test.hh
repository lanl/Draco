//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_ipcress/test/cdi_ipcress_test.hh
 * \author Thomas M. Evans
 * \date   Fri Oct 12 15:36:36 2001
 * \brief  cdi_ipcress test function headers.
 * \note   Copyright (C) 2011-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_ipcress_test_hh
#define rtt_cdi_ipcress_test_hh

#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include <sstream>

// cSpell:ignore IPCRESS

namespace rtt_cdi_ipcress_test {

//------------------------------------------------------------------------------------------------//
// COMPARISON FUNCTIONS USED IN IPCRESS OPACITY TESTS
//------------------------------------------------------------------------------------------------//

template <typename temperatureType, typename densityType, typename testValueType,
          typename opacityClassType>
bool opacityAccessorPassed(rtt_dsxx::ScalarUnitTest &ut, opacityClassType const spOpacity,
                           temperatureType const temperature, densityType const density,
                           testValueType const tabulatedValue) {

  // Interpolate the multigroup opacities.
  testValueType grayOpacity = spOpacity->getOpacity(temperature, density);

  // Make sure that the interpolated value matches previous interpolations.
  std::ostringstream message; // NOLINT(misc-const-correctness) // false-positive
  if (rtt_dsxx::soft_equiv(grayOpacity, tabulatedValue)) {
    message << spOpacity->getDataDescriptor() << " opacity computation was good for \n\t\""
            << spOpacity->getDataFilename() << "\" data.";
    PASSMSG(message.str());
  } else {
    message << spOpacity->getDataDescriptor() << " opacity value is out of spec. for \n\t\""
            << spOpacity->getDataFilename() << "\" data.";
    FAILMSG(message.str());
    return false;
  }

  // If we get here then the test passed.
  return true;
}

//------------------------------------------------------------------------------------------------//

template <typename opacityClassType>
void testTemperatureGridAccessor(rtt_dsxx::ScalarUnitTest &ut, opacityClassType const spOpacity) {

  // Read the temperature grid from the IPCRESS file.
  std::vector<double> const temps = spOpacity->getTemperatureGrid();

  // Verify that the size of the temperature grid looks right.  If it is the right size then compare
  // the temperature grid data to the data specified when we created the IPCRESS file using TOPS.
  if (temps.size() == spOpacity->getNumTemperatures() && temps.size() == 3) {
    PASSMSG("The number of temperature points found in the data grid matches the number returned"
            " by the getNumTemperatures() accessor.");

    // The grid specified by TOPS has 3 temperature points.
    std::vector<double> const temps_ref = {0.1, 1.0, 10.0};
    Check(temps.size() == 3);

    // Compare the grids.
    if (rtt_dsxx::soft_equiv(temps, temps_ref))
      PASSMSG("Temperature grid matches.");
    else
      FAILMSG("Temperature grid did not match.");
  } else {
    FAILMSG("The number of temperature points found in the data grid does not match the number "
            "returned by the getNumTemperatures() accessor.\nDid not test the results returned by "
            "getTemperatureGrid().");
  }
}

//------------------------------------------------------------------------------------------------//

template <typename opacityClassType>
void testDensityGridAccessor(rtt_dsxx::ScalarUnitTest &ut, opacityClassType const spOpacity) {

  // Read the grid from the IPCRESS file.
  std::vector<double> const density = spOpacity->getDensityGrid();

  // Verify that the size of the density grid looks right.  If it is the right size then compare the
  // density grid data to the data specified when we created the IPCRESS file using TOPS.
  if (density.size() == 3 && density.size() == spOpacity->getNumDensities()) {
    PASSMSG("The number of density points found in the data grid matches the number returned by"
            " the getNumDensities() accessor.");

    // The grid specified by TOPS has 3 density points
    std::vector<double> const density_ref = {0.1, 0.5, 1.0};
    Check(density.size() == 3);

    // Compare the grids.
    if (rtt_dsxx::soft_equiv(density, density_ref))
      PASSMSG("Density grid matches.");
    else
      FAILMSG("Density grid did not match.");
  } else {
    FAILMSG("The number of density points found in the data grid does not match the number "
            "returned by the getNumDensities() accessor.\nDid not test the results returned by "
            "getDensityGrid().");
  }
}

//------------------------------------------------------------------------------------------------//
template <typename opacityClassType>
void testEnergyBoundaryAccessor(rtt_dsxx::ScalarUnitTest &ut, opacityClassType const spOpacity) {

  // Read the grid from the IPCRESS file.
  std::vector<double> const eBounds = spOpacity->getGroupBoundaries();

  // Verify that the size of the group boundary grid looks right.  If it is the right size then
  // compare the energy groups grid data to the data specified when we created the IPCRESS file
  // using TOPS.
  if (eBounds.size() == 13 && eBounds.size() == spOpacity->getNumGroupBoundaries()) {
    PASSMSG("The number of energy boundary points found in the data grid matches the number "
            "returned by the getNumGroupBoundaries() accessor.");

    // The grid specified by TOPS has 13 energy boundaries.
    std::vector<double> const eBounds_ref = {0.01, 0.03, 0.07, 0.1,  0.3,  0.7,  1.0,
                                             3.0,  7.0,  10.0, 30.0, 70.0, 100.0};
    Check(eBounds.size() == 13);

    // Compare the grids.
    if (rtt_dsxx::soft_equiv(eBounds, eBounds_ref))
      PASSMSG("Energy group boundary grid matches.");
    else
      FAILMSG("Energy group boundary grid did not match.");
  } else {
    FAILMSG("The number of energy boundary points found in the data grid does not match the "
            "number returned by the get NumGroupBoundaries() accessor.\nDid not test the results "
            "returned by\n\tgetGroupBoundaries().");
  }
}

} // end namespace rtt_cdi_ipcress_test

#endif // rtt_cdi_ipcress_test_hh

//------------------------------------------------------------------------------------------------//
// end of cdi_ipcress/test/cdi_ipcress_test.hh
//------------------------------------------------------------------------------------------------//
