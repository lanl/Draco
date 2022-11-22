//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_analytic/Analytic_EoS.hh
 * \author Thomas M. Evans
 * \date   Tue Oct  2 16:22:32 2001
 * \brief  Analytic_EoS class definition.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_analytic_Analytic_EoS_hh
#define rtt_cdi_analytic_Analytic_EoS_hh

#include "Analytic_Models.hh"
#include "cdi/EoS.hh"
#include <memory>

namespace rtt_cdi_analytic {

//================================================================================================//
/*!
 * \class Analytic_EoS
 *
 * \brief Derived rtt_cdi::EoS class for analytic Equation of State data.
 *
 * The Analytic_EoS class is a derived rtt_cdi::EoS class. It provides analytic Equation of State
 * (EoS) data.  The specific analytic EoS model is derived from the
 * rtt_cdi_analytic::Analytic_EoS_Model base class.  Several pre-built derived classes are provided
 * in Analytic_Models.hh.
 *
 * Clients of this class can provide any analytic model class as long as it conforms to the
 * rtt_cdi_analytic::Analytic_EoS_Model interface.
 *
 * See the member functions for details about the data types and units.
 *
 * \example cdi_analytic/test/tstAnalytic_EoS.cc
 *
 * Example usage of Analytic_EoS, Analytic_EoS_Model, and their incorporation into rtt_cdi::CDI.
 */
//================================================================================================//

class Analytic_EoS : public rtt_cdi::EoS {
public:
  // Useful typedefs.
  using SP_Analytic_Model = std::shared_ptr<Analytic_EoS_Model>;
  using const_SP_Model = std::shared_ptr<const Analytic_EoS_Model>;
  using sf_double = std::vector<double>;
  using sf_char = std::vector<char>;

private:
  // Analytic EoS model.
  SP_Analytic_Model analytic_model;

public:
  // Constructor.
  explicit Analytic_EoS(SP_Analytic_Model model_in);

  // Unpacking constructor.
  explicit Analytic_EoS(const sf_char &packed);

  // >>> ACCESSORS
  const_SP_Model get_Analytic_Model() const { return analytic_model; }

  // >>> INTERFACE SPECIFIED BY rtt_cdi::EoS

  // Get electron internal energy.
  double getSpecificElectronInternalEnergy(double T, double rho) const override;
  sf_double getSpecificElectronInternalEnergy(const sf_double &T,
                                              const sf_double &rho) const override;

  // Get ion internal energy.
  double getSpecificIonInternalEnergy(double T, double rho) const override;
  sf_double getSpecificIonInternalEnergy(const sf_double &T, const sf_double &rho) const override;

  // Get electron heat capacity.
  double getElectronHeatCapacity(double T, double rho) const override;
  sf_double getElectronHeatCapacity(const sf_double &T, const sf_double &rho) const override;

  // Get ion heat capacity.
  double getIonHeatCapacity(double T, double rho) const override;
  sf_double getIonHeatCapacity(const sf_double &T, const sf_double &rho) const override;

  // Get the number of free electrons per ion.
  double getNumFreeElectronsPerIon(double T, double rho) const override;
  sf_double getNumFreeElectronsPerIon(const sf_double &T, const sf_double &rho) const override;

  // Get the electron thermal conductivity.
  double getElectronThermalConductivity(double T, double rho) const override;
  sf_double getElectronThermalConductivity(const sf_double &T, const sf_double &rho) const override;

  // Get the new Te, given delta Ue, Te0.
  double getElectronTemperature(double /*rho*/, double Ue, double Tguess = 1.0) const override;

  // Get the new Ti, given delta Uic, Ti0.
  double getIonTemperature(double /*rho*/, double Uic, double Tguess = 1.0) const override;

  // Pack the Analytic_EoS into a character string.
  sf_char pack() const override;
};

} // end namespace rtt_cdi_analytic

#endif // rtt_cdi_analytic_Analytic_EoS_hh

//------------------------------------------------------------------------------------------------//
// end of cdi_analytic/Analytic_EoS.hh
//------------------------------------------------------------------------------------------------//
