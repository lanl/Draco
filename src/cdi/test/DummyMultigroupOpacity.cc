//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi/test/DummyMultigroupOpacity.cc
 * \author Kelly Thompson
 * \date   Mon Jan 8 15:17:16 2001
 * \brief  DummyMultigroupOpacity templated class implementation file.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "DummyMultigroupOpacity.hh"
#include <cmath> // pow(x,n)

namespace rtt_cdi_test {

// ------------ //
// Constructors //
// ------------ //

/*!
 * \brief Constructor for DummyMultigroupOpacity object.
 *
 * \sa The constructor assigns fixed values for all of the member data.  Every instance of this
 *     object has the same member data.
 *
 *     Temperatures     = { 1.0, 2.0, 3.0 }
 *     Densities        = { 0.1, 0.2 }
 *     EnergyBoundaries = { 0.05, 0.5, 5.0, 50.0 }
 */
DummyMultigroupOpacity::DummyMultigroupOpacity(rtt_cdi::Reaction reaction, rtt_cdi::Model model)
    : dataFilename("none"), dataDescriptor("DummyMultigroupOpacity"),
      energyPolicyDescriptor("Multigroup"), groupBoundaries(), temperatureGrid(), densityGrid(),
      reaction_type(reaction), model_type(model) {
  temperatureGrid.resize(numTemperatures);
  densityGrid.resize(numDensities);
  groupBoundaries.resize(numGroupBoundaries);

  for (size_t i = 0; i < numTemperatures; ++i)
    temperatureGrid[i] = static_cast<double>(i + 1);
  for (size_t i = 0; i < numDensities; ++i)
    densityGrid[i] = static_cast<double>(i + 1) * 0.1;
  for (size_t i = 0; i < numGroupBoundaries; ++i)
    groupBoundaries[i] = 5.0 * std::pow(10.0, static_cast<double>(i) - 2.0);
}

// Constructor for entering a different group boundary structure than the default
DummyMultigroupOpacity::DummyMultigroupOpacity(rtt_cdi::Reaction reaction, rtt_cdi::Model model,
                                               size_t num_boundaries)
    : dataFilename("none"), dataDescriptor("DummyMultigroupOpacity"),
      energyPolicyDescriptor("Multigroup"), numGroupBoundaries(num_boundaries), groupBoundaries(),
      temperatureGrid(), densityGrid(), reaction_type(reaction), model_type(model) {
  temperatureGrid.resize(numTemperatures);
  densityGrid.resize(numDensities);
  groupBoundaries.resize(numGroupBoundaries);

  for (size_t i = 0; i < numTemperatures; ++i)
    temperatureGrid[i] = static_cast<double>(i + 1);
  for (size_t i = 0; i < numDensities; ++i)
    densityGrid[i] = static_cast<double>(i + 1) * 0.1;
  for (size_t i = 0; i < numGroupBoundaries; ++i)
    groupBoundaries[i] = 5.0 * std::pow(10.0, static_cast<double>(i) - 2.0);
}

// --------- //
// Accessors //
// --------- //

/*!
 * \brief Opacity accessor that returns a vector of opacities (one for each group) that corresponds
 *     to the provided temperature and density.
 *
 *     Opacity = 2 * ( temperature + density/1000 ) / ( E_high + E_low )
 */
std::vector<double> DummyMultigroupOpacity::getOpacity(double targetTemperature,
                                                       double targetDensity) const {
  std::vector<double> opacity(numGroupBoundaries - 1);
  for (size_t ig = 0; ig < numGroupBoundaries - 1; ++ig)
    opacity[ig] = 2.0 * (targetTemperature + targetDensity / 1000.0) /
                  (groupBoundaries[ig] + groupBoundaries[ig + 1]);
  return opacity;
}

/*!
 * \brief Opacity accessor that returns a vector of multigroup opacities corresponding to the
 *     provided vector of temperatures and a single density.  Each multigroup opacity is in itself a
 *     vector of numGroups opacities.
 *
 *     Opacity = 2 * ( temperature + density/1000 ) / ( E_high + E_low )
 */
std::vector<std::vector<double>>
DummyMultigroupOpacity::getOpacity(const std::vector<double> &targetTemperature,
                                   double targetDensity) const {
  std::vector<std::vector<double>> opacity(targetTemperature.size());

  for (size_t it = 0; it < targetTemperature.size(); ++it) {
    opacity[it].resize(numGroupBoundaries - 1);

    for (size_t ig = 0; ig < numGroupBoundaries - 1; ++ig)
      opacity[it][ig] = 2.0 * (targetTemperature[it] + targetDensity / 1000.0) /
                        (groupBoundaries[ig] + groupBoundaries[ig + 1]);
  }

  return opacity;
}

/*!
 * \brief Opacity accessor that returns a vector of multigroup opacities corresponding to the
 *     provided vector of densities and a single temperature.  Each multigroup opacity is in itself
 *     a vector of numGroups opacities.
 *
 *     Opacity = 2 * ( temperature + density/1000 ) / ( E_high + E_low )
 */
std::vector<std::vector<double>>
DummyMultigroupOpacity::getOpacity(double targetTemperature,
                                   const std::vector<double> &targetDensity) const {
  std::vector<std::vector<double>> opacity(targetDensity.size());

  for (size_t id = 0; id < targetDensity.size(); ++id) {
    opacity[id].resize(numGroupBoundaries - 1);

    for (size_t ig = 0; ig < numGroupBoundaries - 1; ++ig)
      opacity[id][ig] = 2.0 * (targetTemperature + targetDensity[id] / 1000.0) /
                        (groupBoundaries[ig] + groupBoundaries[ig + 1]);
  }

  return opacity;
}

} // end namespace rtt_cdi_test

//------------------------------------------------------------------------------------------------//
// end of DummyMultigroupOpacity.cc
//------------------------------------------------------------------------------------------------//
