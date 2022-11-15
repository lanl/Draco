//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi/test/DummyEoS.cc
 * \author Kelly Thompson
 * \date   Mon Jan 8 16:25:09 2001
 * \brief  DummyEoS class header file (derived from ../EoS)
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "DummyEoS.hh"
#include <cmath>

namespace rtt_cdi_test {

// --------- //
// Accessors //
// --------- //

//------------------------------------------------------------------------------------------------//
/*!
 * \brief EoS accessor that returns a single specific electron internal energy that corresponds to
 *        the provided temperature and density.
 *
 *    internal energy = temperature + 1000*density
 *
 * \param temperature The temperature value for which an opacity value is being requested (Kelvin).
 * \param density The density value for which an opacity value is being requested (g/cm^3).
 * \return A specific electron internal energy (kJ/g).
 */
double DummyEoS::getSpecificElectronInternalEnergy(double temperature, double density) const {
  return temperature + 1000.0 * density;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief EoS accessor that returns a vector of specific electron internal energies that correspond
 *     to the provided vectors of temperatures and densities.
 *
 *    internal energy[i] = temperature[i] + 1000*density[i]
 *
 * \param vtemperature A vector of temperature values for which the EoS values are being requested
 *           (Kelvin).
 * \param vdensity A vector of density values for which the EoS values are being requested (g/cm^3).
 * \return A vector of specific electron internal energies (kJ/g).
 */
std::vector<double>
DummyEoS::getSpecificElectronInternalEnergy(const std::vector<double> &vtemperature,
                                            const std::vector<double> &vdensity) const {
  std::vector<double> seie(vdensity.size());
  for (unsigned i = 0; i < vdensity.size(); ++i)
    seie[i] = vtemperature[i] + 1000.0 * vdensity[i];
  return seie;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Retrieve the electron based heat capacity for this material at the provided density and
 *        temperature.
 *
 *     heat capacity = temperature + density/1000
 *
 * \param density Density of the material in g/cm^3
 * \param temperature Temperature of the material in Kelvin.
 * \return The electron based heat capacity in kJ/g/K.
 */
double DummyEoS::getElectronHeatCapacity(double temperature, double density) const {
  return temperature + density / 1000.0;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Retrieve a set of electron based heat capacities for this material that correspond to the
 *        tuple list of provided densities and temperatures.
 *
 *     heat capacity = temperature + density/1000
 *
 * \param vdensity Density of the material in g/cm^3
 * \param vtemperature Temperature of the material in Kelvin.
 * \return The electron based heat capacity in kJ/g/K.
 */
std::vector<double> DummyEoS::getElectronHeatCapacity(const std::vector<double> &vtemperature,
                                                      const std::vector<double> &vdensity) const {
  std::vector<double> ehc(vdensity.size());
  for (unsigned i = 0; i < vdensity.size(); ++i)
    ehc[i] = vtemperature[i] + vdensity[i] / 1000.0;
  return ehc;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Retrieve the specific ion internal energy for this material at the provided density and
 *        temperature.
 *
 *     internal energy = temperature + 1000*density
 *
 * \param density Density of the material in g/cm^3
 * \param temperature Temperature of the material in Kelvin.
 * \return The specific ion internal energy in kJ/g.
 */
double DummyEoS::getSpecificIonInternalEnergy(double temperature, double density) const {
  return getSpecificElectronInternalEnergy(density, temperature);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Retrieve a set of specific ion internal energies for this material that correspond to the
 *        tuple list of provided densities and temperatures.
 *
 *     internal energy = temperature + 1000*density
 *
 * \param vdensity Density of the material in g/cm^3
 * \param vtemperature Temperature of the material in Kelvin.
 * \return A vector of specific ion internal energies in kJ/g.
 */
std::vector<double>
DummyEoS::getSpecificIonInternalEnergy(const std::vector<double> &vtemperature,
                                       const std::vector<double> &vdensity) const {
  return getSpecificElectronInternalEnergy(vdensity, vtemperature);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Retrieve the ion based heat capacity for this material at the provided density and
 *        temperature.
 *
 *     heat capacity   = temperature + density/1000
 *
 * \param density Density of the material in g/cm^3
 * \param temperature Temperature of the material in Kelvin.
 * \return The ion based heat capacity in kJ/g/K.
 */
double DummyEoS::getIonHeatCapacity(double temperature, double density) const {
  return getElectronHeatCapacity(density, temperature);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Retrieve a set of ion based heat capacities for this material that correspond to the tuple
 *        list of provided densities and temperatures.
 *
 *     heat capacity   = temperature + density/1000
 *
 * \param vdensity Density of the material in g/cm^3
 * \param vtemperature Temperature of the material in Kelvin.
 * \return A vector of ion based heat capacities in kJ/g/K.
 */
std::vector<double> DummyEoS::getIonHeatCapacity(const std::vector<double> &vtemperature,
                                                 const std::vector<double> &vdensity) const {
  return getElectronHeatCapacity(vdensity, vtemperature);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Retrieve the number of free electrons per ion for this material at the provided density
 *        and temperature.
 *
 *     num free electrons = temperature/100
 *
 * \param density Density of the material in g/cm^3
 * \param temperature Temperature of the material in Kelvin.
 * \return The number of free electrons per ion.
 */
double DummyEoS::getNumFreeElectronsPerIon(double temperature, double /*density*/) const {
  return temperature / 100.0;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Retrieve a set of free electrons per ion averages for this material that correspond to the
 *        tuple list of provided densities and temperatures.
 *
 * \param vdensity Density of the material in g/cm^3
 * \param vtemperature Temperature of the material in Kelvin.
 * \return A vector of the number of free electrons per ion.
 */
std::vector<double> DummyEoS::getNumFreeElectronsPerIon(const std::vector<double> &vtemperature,
                                                        const std::vector<double> &vdensity) const {
  std::vector<double> nfepi(vdensity.size());
  for (unsigned i = 0; i < vdensity.size(); ++i)
    nfepi[i] = vtemperature[i] / 100.0;
  return nfepi;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Retrieve the electron based thermal conductivity for this material at the provided density
 *        and temperature.
 *
 *     thermal conductivity = 1000*temperature + density
 *
 * \param density Density of the material in g/cm^3
 * \param temperature Temperature of the material in Kelvin.
 * \return The electron based thermal conductivity in 1/s/cm.
 */
double DummyEoS::getElectronThermalConductivity(double temperature, double density) const {
  return 1000.0 * temperature + density;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Retrieve a set of electron based thermal conductivities for this material that correspond
 *        to the tuple list of provided densities and temperatures.
 *
 *     thermal conductivity = 1000*temperature + density
 *
 * \param vdensity Density of the material in g/cm^3
 * \param vtemperature Temperature of the material in Kelvin.
 * \return A vector of electron based thermal conductivities in 1/s/cm.
 */
std::vector<double>
DummyEoS::getElectronThermalConductivity(const std::vector<double> &vtemperature,
                                         const std::vector<double> &vdensity) const {
  std::vector<double> ebtc(vdensity.size());
  for (unsigned i = 0; i < vdensity.size(); ++i)
    ebtc[i] = 1000.0 * vtemperature[i] + vdensity[i];
  return ebtc;
}

//------------------------------------------------------------------------------------------------//
//! Return Electron Temperature (keV) given Specific Electron Internal Energy (kJ/g).
double DummyEoS::getElectronTemperature(double /*density*/, double SpecificElectronInternalEnergy,
                                        double /*Tguess*/) const {
  return 10.0 * SpecificElectronInternalEnergy;
}

//! Return Ion Temperature (keV) given Specific Ion Internal Energy (kJ/g).
double DummyEoS::getIonTemperature(double /*density*/, double SpecificIonInternalEnergy,
                                   double /*Tguess*/) const {
  return std::sqrt(SpecificIonInternalEnergy);
}

} // end namespace rtt_cdi_test

//------------------------------------------------------------------------------------------------//
// end of cdi/test/DummyEoS.hh
//------------------------------------------------------------------------------------------------//
