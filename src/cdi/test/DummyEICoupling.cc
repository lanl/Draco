//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi/test/DummyEICoupling.cc
 * \author Mathew Cleveland
 * \date   March 2019
 * \brief  DummyEICoupling class header file (derived from ../EICoupling)
 * \note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "DummyEICoupling.hh"
#include <cmath>

namespace rtt_cdi_test {

// --------- //
// Accessors //
// --------- //

//------------------------------------------------------------------------------------------------//
/*!
 * \brief EICoupling accessor that returns a single electron-ion coupling given an electron and ion
 *        temperature, the material density, and the electron and ion screening coeffiecients..
 *
 *    dummy_ei_coupling = etemperature + 10*itemperature + 100*density + 1000*w_e + 10000*w_i
 *
 * \param[in] etemperature The electron temperature value for which an opacity value is being
 *        requested (Kelvin).
 * \param[in] itemperature The ion temperature value for which an opacity value is being requested
 *        (Kelvin).
 * \param[in] density The density value for which an opacity value is being requested (g/cm^3).
 * \param[in] w_e the electron screening coeffiecent [1/s]
 * \param[in] w_i the ion screening coeffiecent [1/s]
 * \return An electron-ion coupling coeffient [kJ/cc/keV/s].
 */
double DummyEICoupling::getElectronIonCoupling(const double etemperature, const double itemperature,
                                               const double density, const double w_e,
                                               const double w_i) const {
  return etemperature + 10.0 * itemperature + 100.0 * density + 1000.0 * w_e + 10000.0 * w_i;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief EICoupling accessor that returns a vector of electron-ion coupling given an electron and
 *        ion temperature, the material density, and the electron and ion screening coeffiecients..
 *
 *    dummy_ei_coupling = etemperature + 10*itemperature + 100*density + 1000*w_e + 10000*w_i
 *
 * \param[in] etemperature The electron temperature vector for which an opacity value is being
 *        requested (Kelvin).
 * \param[in] itemperature The ion temperature vector for which an opacity value is being requested
 *        (Kelvin).
 * \param[in] density The density vector for which an opacity value is being requested (g/cm^3).
 * \param[in] w_e the electron screening coeffiecent vector [1/s]
 * \param[in] w_i the ion screening coeffiecent vector [1/s]
 * \return An electron-ion coupling coeffient vector [kJ/cc/keV/s].
 */
std::vector<double> DummyEICoupling::getElectronIonCoupling(const std::vector<double> &etemperature,
                                                            const std::vector<double> &itemperature,
                                                            const std::vector<double> &density,
                                                            const std::vector<double> &w_e,
                                                            const std::vector<double> &w_i) const {
  std::vector<double> ei_coupling(density.size());
  for (unsigned i = 0; i < density.size(); ++i)
    ei_coupling[i] = etemperature[i] + 10.0 * itemperature[i] + 100.0 * density[i] +
                     1000.0 * w_e[i] + 10000.0 * w_i[i];
  return ei_coupling;
}

} // end namespace rtt_cdi_test

//------------------------------------------------------------------------------------------------//
// end of cdi/test/DummyEICoupling.hh
//------------------------------------------------------------------------------------------------//
