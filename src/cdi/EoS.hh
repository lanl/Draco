//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi/EoS.hh
 * \author Kelly Thompson
 * \date   Fri Apr 13 16:15:59 2001
 * \brief  EoS class header file (an abstract class)
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_EoS_hh
#define rtt_cdi_EoS_hh

#include "ds++/config.h"
#include <vector>

namespace rtt_cdi {

//================================================================================================//
/*!
 * \class EoS
 *
 * \brief This is a pure virtual class that defines a standard interface for all derived EoS
 *        objects.
 *
 * Any derived EoS object must provide as a minumum the functionality outlined in this routine.
 * This functionality includes access to the data grid and the ability to return interpolated
 * opacity values.
 *
 * \example cdi/test/tDummyEoS.cc
 * \sa cdi/test/tCDI.cc
 */
//================================================================================================//

class EoS {
  // DATA

  // There is no data for a pure virtual object.  This class provides an interface and does not
  // preserve state.

public:
  // ---------- //
  // Destructor //
  // ---------- //

  /*!
   * \brief Default EoS() destructor.
   *
   * This is required to correctly release memory when any object derived from EoS is destroyed.
   */
  virtual ~EoS() = default;
  //! Default constructor
  EoS() = default;
  //! Disable copy construction
  EoS(EoS const &rhs) = delete;
  //! Disable move construction
  EoS(EoS &&rhs) noexcept = delete;
  //! Disable assignment
  EoS &operator=(EoS const &rhs) = delete;
  //! Disable move-assignment
  EoS &operator=(EoS &&rhs) noexcept = delete;

  // --------- //
  // Accessors //
  // --------- //

  /*!
   * \brief EoS accessor that returns a single specific electron internal energy that corresponds to
   *        the provided temperature and density.
   *
   * \param temperature The temperature value for which an opacity value is being requested (keV).
   * \param density The density value for which an opacity value is being requested (g/cm^3).
   * \return A specific electron internal energy (kJ/g).
   */
  virtual double getSpecificElectronInternalEnergy(double temperature, double density) const = 0;

  /*!
   * \brief EoS accessor that returns a vector of specific electron internal energies that
   *        correspond to the provided vectors of temperatures and densities.
   *
   * \param vtemperature A vector of temperature values for which the EoS values are being requested
   *           (keV).
   * \param vdensity A vector of density values for which the EoS values are being requested
   *           (g/cm^3).
   * \return A vector of specific electron internal energies (kJ/g).
   */
  virtual std::vector<double>
  getSpecificElectronInternalEnergy(const std::vector<double> &vtemperature,
                                    const std::vector<double> &vdensity) const = 0;

  /*!
   * \brief Retrieve the electron based heat capacity for this material at the provided density and
   *        temperature.
   *
   * \param temperature Temperature of the material in keV.
   * \param density Density of the material in g/cm^3
   * \return The electron based heat capacity in kJ/g/keV.
   */
  virtual double getElectronHeatCapacity(double temperature, double density) const = 0;

  /*!
   * \brief Retrieve a set of electron based heat capacities for this material that correspond to
   *        the tuple list of provided densities and temperatures.
   *
   * \param vtemperature Temperature of the material in keV.
   * \param vdensity Density of the material in g/cm^3
   * \return The electron based heat capacity in kJ/g/keV.
   */
  virtual std::vector<double>
  getElectronHeatCapacity(const std::vector<double> &vtemperature,
                          const std::vector<double> &vdensity) const = 0;

  /*!
   * \brief Retrieve the specific ion internal energy for this material at the provided density and
   *        temperature.
   *
   * \param density Density of the material in g/cm^3
   * \param temperature Temperature of the material in keV.
   * \return The specific ion internal energy in kJ/g.
   */
  virtual double getSpecificIonInternalEnergy(double temperature, double density) const = 0;

  /*!
   * \brief Retrieve a set of specific ion internal energies for this material that correspond to
   *        the tuple list of provided densities and temperatures.
   *
   * \param vdensity Density of the material in g/cm^3
   * \param vtemperature Temperature of the material in keV.
   * \return A vector of specific ion internal energies in kJ/g.
   */
  virtual std::vector<double>
  getSpecificIonInternalEnergy(const std::vector<double> &vtemperature,
                               const std::vector<double> &vdensity) const = 0;

  /*!
   * \brief Retrieve the ion based heat capacity for this material at the provided density and
   *        temperature.
   *
   * \param density Density of the material in g/cm^3
   * \param temperature Temperature of the material in keV.
   * \return The ion based heat capacity in kJ/g/keV.
   */
  virtual double getIonHeatCapacity(double temperature, double density) const = 0;

  /*!
   * \brief Retrieve a set of ion based heat capacities for this material that correspond to the
   *        tuple list of provided densities and temperatures.
   *
   * \param vdensity Density of the material in g/cm^3
   * \param vtemperature Temperature of the material in keV.
   * \return A vector of ion based heat capacities in kJ/g/keV.
   */
  virtual std::vector<double> getIonHeatCapacity(const std::vector<double> &vtemperature,
                                                 const std::vector<double> &vdensity) const = 0;

  /*!
   * \brief Retrieve the number of free electrons per ion for this material at the provided density
   *        and temperature.
   *
   * \param density Density of the material in g/cm^3
   * \param temperature Temperature of the material in keV.
   * \return The number of free electrons per ion.
   */
  virtual double getNumFreeElectronsPerIon(double temperature, double density) const = 0;

  /*!
   * \brief Retrieve a set of free electrons per ion averages for this material that correspond to
   *        the tuple list of provided densities and temperatures.
   *
   * \param vdensity Density of the material in g/cm^3
   * \param vtemperature Temperature of the material in keV.
   * \return A vector of the number of free electrons per ion.
   */
  virtual std::vector<double>
  getNumFreeElectronsPerIon(const std::vector<double> &vtemperature,
                            const std::vector<double> &vdensity) const = 0;

  /*!
   * \brief Retrieve the electron based thermal conductivity for this
   *        material at the provided density and temperature.
   *
   * \param density Density of the material in g/cm^3
   * \param temperature Temperature of the material in keV.
   * \return The electron based thermal conductivity in 1/s/cm.
   */
  virtual double getElectronThermalConductivity(double temperature, double density) const = 0;

  /*!
   * \brief Retrieve a set of electron based thermal conductivities for this material that
   *        correspond to the tuple list of provided densities and temperatures.
   *
   * \param vdensity Density of the material in g/cm^3
   * \param vtemperature Temperature of the material in keV.
   * \return A vector of electron based thermal conductivities in 1/s/cm.
   */
  virtual std::vector<double>
  getElectronThermalConductivity(const std::vector<double> &vtemperature,
                                 const std::vector<double> &vdensity) const = 0;

  /*!
   * \brief Retrieve an electron temperature associated with the provided specific electron internal
   *        energy (kJ/g) and density (g/cm^3).
   *
   * \param density Density of the material in g/cm^3
   * \param SpecificElectronInternalEnergy in kJ/g
   * \param Tguess  A guess for the resulting electron temperature to aid the root finder.
   * \return An electron (material) temperature in keV.
   */
  virtual double getElectronTemperature(double density, double SpecificElectronInternalEnergy,
                                        double Tguess) const = 0;

  /*!
   * \brief Retrieve an ion temperature associated with the provided specific ion internal energy
   *        (kJ/g) and density (g/cm^3).
   *
   * \param density Density of the material in g/cm^3
   * \param SpecificIonInternalEnergy in kJ/g
   * \param Tguess  A guess for the resulting ion temperature to aid the root finder.
   * \return Ionic temperature in keV.
   */
  virtual double getIonTemperature(double density, double SpecificIonInternalEnergy,
                                   double Tguess) const = 0;

  /*!
   * \brief Interface for packing a derived EoS object.
   *
   * Note, the user hands the return value from this function to a derived EoS constructor.  Thus,
   * even though one can pack a EoS through a base class pointer, the client must know the derived
   * type when unpacking.
   */
  virtual std::vector<char> pack() const = 0;

}; // end of class EoS

} // end namespace rtt_cdi

#endif // rtt_cdi_EoS_hh

//------------------------------------------------------------------------------------------------//
// end of cdi/EoS.hh
//------------------------------------------------------------------------------------------------//
