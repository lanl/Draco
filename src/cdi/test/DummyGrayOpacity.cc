//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi/test/DummyGrayOpacity.cc
 * \author Kelly Thompson
 * \date   Mon Jan 8 15:33:51 2001
 * \brief  DummyGrayOpacity class implementation file.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "DummyGrayOpacity.hh"

namespace rtt_cdi_test {

// ------------ //
// Constructors //
// ------------ //

/*!
 * \brief Constructor for DummyGrayOpacity object.
 *
 * The constructor assigns fixed values for all of the member data.  Every instance of this object
 * has the same member data.
 *
 * Temperatures = { 1.0, 2.0, 3.0 }
 * Densities    = { 0.1, 0.2 }
 */
DummyGrayOpacity::DummyGrayOpacity(rtt_cdi::Reaction reaction, rtt_cdi::Model model)
    : dataFilename("none"), dataDescriptor("DummyGrayOpacity"), energyPolicyDescriptor("Gray"),
      temperatureGrid(), densityGrid(), reaction_type(reaction), model_type(model) {
  // Set up the temperature and density grid.
  temperatureGrid.resize(numTemperatures);
  densityGrid.resize(numDensities);
  for (size_t i = 0; i < numTemperatures; ++i)
    temperatureGrid[i] = static_cast<double>(i + 1);
  for (size_t i = 0; i < numDensities; ++i)
    densityGrid[i] = 0.1 * static_cast<double>(i + 1);
}

// --------- //
// Accessors //
// --------- //

/*!
 * \brief Opacity accessor that returns a single opacity (or a vector of opacities for the
 *        multigroup EnergyPolicy) that corresponds to the provided temperature and density.
 *
 *     Opacity = temperature + density/1000
 */
double DummyGrayOpacity::getOpacity(double targetTemperature, double targetDensity) const {
  return targetTemperature + targetDensity / 1000.0;
}

/*!
 * \brief Opacity accessor that returns a vector of opacities that correspond to the provided vector
 *        of temperatures and a single density value.
 *
 *     Opacity[i] = temperature[i] + density/1000
 */
std::vector<double> DummyGrayOpacity::getOpacity(const std::vector<double> &targetTemperature,
                                                 double targetDensity) const {
  std::vector<double> grayOpacity(targetTemperature.size());
  for (size_t i = 0; i < targetTemperature.size(); ++i)
    grayOpacity[i] = targetTemperature[i] + targetDensity / 1000.0;
  return grayOpacity;
}

/*!
 * \brief Opacity accessor that returns a vector of opacities that correspond to the provided vector
 *        of densities and a single temperature value.
 *
 *     Opacity[i] = temperature[i] + density/1000
 */
std::vector<double> DummyGrayOpacity::getOpacity(double targetTemperature,
                                                 const std::vector<double> &targetDensity) const {
  std::vector<double> grayOpacity(targetDensity.size());
  for (size_t i = 0; i < targetDensity.size(); ++i)
    grayOpacity[i] = targetTemperature + targetDensity[i] / 1000.0;
  return grayOpacity;
}

} // namespace rtt_cdi_test

//------------------------------------------------------------------------------------------------//
// end of DummyGrayOpacity.cc
//------------------------------------------------------------------------------------------------//
