//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_analytic/Analytic_EoS.cc
 * \author Thomas M. Evans
 * \date   Tue Oct  2 16:22:32 2001
 * \brief  Analytic_EoS member definitions.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Analytic_EoS.hh"
#include "ds++/Packing_Utils.hh"

namespace rtt_cdi_analytic {

//------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for an analytic EoS model.
 *
 * This constructor builds an analytic EoS model defined by the
 * rtt_cdi_analytic::Analytic_EoS_Model derived class argument.
 *
 * \param model_in shared_ptr to a derived
 *        rtt_cdi_analytic::Analytic_EoS_Model object
 */
Analytic_EoS::Analytic_EoS(SP_Analytic_Model model_in) : analytic_model(std::move(model_in)) {
  Ensure(analytic_model);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Unpacking constructor.
 *
 * This constructor rebuilds and Analytic_EoS from a vector<char> that was
 * created by a call to pack().  It can only rebuild Analytic_Model types
 * that have been registered in the rtt_cdi_analytic::EoS_Models enumeration.
 */
Analytic_EoS::Analytic_EoS(const sf_char &packed) : analytic_model() {
  // the packed size must be at least 2 integers (size,
  // analytic model indicator)
  Require(packed.size() >= 2 * sizeof(int));

  // make an unpacker
  rtt_dsxx::Unpacker unpacker;

  // set the buffer
  unpacker.set_buffer(packed.size(), &packed[0]);

  // unpack the size of the analytic model
  int size_analytic;
  unpacker >> size_analytic;
  Check(static_cast<size_t>(size_analytic) >= sizeof(int));

  // unpack the packed analytic model
  std::vector<char> packed_analytic(size_analytic);
  for (size_t i = 0; i < static_cast<size_t>(size_analytic); ++i)
    unpacker >> packed_analytic[i];

  Check(unpacker.get_ptr() == &packed[0] + packed.size());

  // now reset the buffer so that we can determine the analytic model
  // indicator
  unpacker.set_buffer(size_analytic, &packed_analytic[0]);

  // unpack the indicator
  int indicator;
  unpacker >> indicator;

  // now determine which analytic model we need to build
  if (indicator == POLYNOMIAL_SPECIFIC_HEAT_ANALYTIC_EOS_MODEL) {
    analytic_model = std::make_shared<Polynomial_Specific_Heat_Analytic_EoS_Model>(packed_analytic);
  } else {
    Insist(0, "Unregistered analytic EoS model!");
  }

  Ensure(analytic_model);
}

//------------------------------------------------------------------------------------------------//
// EoS INTERFACE FUNCTIONS
//------------------------------------------------------------------------------------------------//
/*!
 *
 * \brief Return a scalar electron internal energy given a scalar temperature
 * and density.
 *
 * Given a scalar temperature and density, return the electron internal
 * energy.  The electron internal energy is defined by the analytic model
 * given to the constructor (Analytic_EoS()).
 *
 * \param T material temperature in keV
 * \param rho density in g/cm^3
 * \return electron internal energy in kJ/g
 */
double Analytic_EoS::getSpecificElectronInternalEnergy(double T, double rho) const {
  Require(T >= 0.0);
  Require(rho >= 0.0);

  double internal_energy = analytic_model->calculate_electron_internal_energy(T, rho);

  Ensure(internal_energy >= 0.0);
  return internal_energy;
}

//------------------------------------------------------------------------------------------------//
/*!
 *
 * \brief Return a field of electron internal energies given fields of
 * temperature and density.
 *
 * Given temperature and density fields, return an electron internal energy
 * field.  The electron internal energy field is the same size as the input
 * temperature and density fields. The electron internal energy is defined by
 * the analytic model given to the constructor (Analytic_EoS()).
 *
 * \param T material temperature field in keV
 * \param rho density field in g/cm^3
 * \return electron internal energy field in kJ/g
 */
Analytic_EoS::sf_double
Analytic_EoS::getSpecificElectronInternalEnergy(const sf_double &T, const sf_double &rho) const {
  Require(T.size() == rho.size());

  // define the return electron internal energy field
  sf_double internal_energy(T.size(), 0.0);

  // loop through T/rho field and solve for internal energy
  for (size_t i = 0; i < T.size(); i++) {
    Check(T[i] >= 0.0);
    Check(rho[i] >= 0.0);

    internal_energy[i] = analytic_model->calculate_electron_internal_energy(T[i], rho[i]);

    Check(internal_energy[i] >= 0.0);
  }

  return internal_energy;
}

//------------------------------------------------------------------------------------------------//
/*!
 *
 * \brief Return a scalar ion internal energy given a scalar temperature and
 * density.
 *
 * Given a scalar temperature and density, return the ion internal energy.
 * The ion internal energy is defined by the analytic model given to the
 * constructor (Analytic_EoS()).
 *
 * \param T material temperature in keV
 * \param rho density in g/cm^3
 * \return ion internal energy in kJ/g
 */
double Analytic_EoS::getSpecificIonInternalEnergy(double T, double rho) const {
  Require(T >= 0.0);
  Require(rho >= 0.0);

  double internal_energy = analytic_model->calculate_ion_internal_energy(T, rho);

  Ensure(internal_energy >= 0.0);
  return internal_energy;
}

//------------------------------------------------------------------------------------------------//
/*!
 *
 * \brief Return a field of ion internal energies given fields of temperature
 * and density.
 *
 * Given temperature and density fields, return an ion internal energy field.
 * The ion internal energy field is the same size as the input temperature
 * and density fields. The ion internal energy is defined by the analytic
 * model given to the constructor (Analytic_EoS()).
 *
 * \param T material temperature field in keV
 * \param rho density field in g/cm^3
 * \return ion internal energy field in kJ/g
 */
Analytic_EoS::sf_double Analytic_EoS::getSpecificIonInternalEnergy(const sf_double &T,
                                                                   const sf_double &rho) const {
  Require(T.size() == rho.size());

  // define the return ion internal energy field
  sf_double internal_energy(T.size(), 0.0);

  // loop through T/rho field and solve for internal energy
  for (size_t i = 0; i < T.size(); i++) {
    Check(T[i] >= 0.0);
    Check(rho[i] >= 0.0);

    internal_energy[i] = analytic_model->calculate_ion_internal_energy(T[i], rho[i]);

    Check(internal_energy[i] >= 0.0);
  }

  return internal_energy;
}

//------------------------------------------------------------------------------------------------//
/*!
 *
 * \brief Return a scalar electron heat capacity given a scalar temperature
 * and density.
 *
 * Given a scalar temperature and density, return the electron heat capacity.
 * The electron heat capacity is defined by the analytic model given to the
 * constructor (Analytic_EoS()).
 *
 * \param T material temperature in keV
 * \param rho density in g/cm^3
 * \return electron heat capacity in kJ/g/keV
 */
double Analytic_EoS::getElectronHeatCapacity(double T, double rho) const {
  Require(T >= 0.0);
  Require(rho >= 0.0);

  double heat_capacity = analytic_model->calculate_electron_heat_capacity(T, rho);

  Ensure(heat_capacity >= 0.0);
  return heat_capacity;
}

//------------------------------------------------------------------------------------------------//
/*!
 *
 * \brief Return a field of electron heat capacities given fields of
 * temperature and density.
 *
 * Given temperature and density fields, return an electron heat capacity
 * field.  The electron heat capacity field is the same size as the input
 * temperature and density fields. The electron heat capacity is defined by
 * the analytic model given to the constructor (Analytic_EoS()).
 *
 * \param T material temperature field in keV
 * \param rho density field in g/cm^3
 * \return electron heat capacity field in kJ/g/keV
 */
Analytic_EoS::sf_double Analytic_EoS::getElectronHeatCapacity(const sf_double &T,
                                                              const sf_double &rho) const {
  Require(T.size() == rho.size());

  // define the return electron heat capacity field
  sf_double heat_capacity(T.size(), 0.0);

  // loop through T/rho field and solve for the heat capacity
  for (size_t i = 0; i < T.size(); i++) {
    Check(T[i] >= 0.0);
    Check(rho[i] >= 0.0);

    heat_capacity[i] = analytic_model->calculate_electron_heat_capacity(T[i], rho[i]);

    Check(heat_capacity[i] >= 0.0);
  }

  return heat_capacity;
}

//------------------------------------------------------------------------------------------------//
/*!
 *
 * \brief Return a scalar ion heat capacity given a scalar temperature and
 * density.
 *
 * Given a scalar temperature and density, return the ion heat capacity.  The
 * ion heat capacity is defined by the analytic model given to the
 * constructor (Analytic_EoS()).
 *
 * \param T material temperature in keV
 * \param rho density in g/cm^3
 * \return ion heat capacity in kJ/g/keV
 */
double Analytic_EoS::getIonHeatCapacity(double T, double rho) const {
  Require(T >= 0.0);
  Require(rho >= 0.0);

  double heat_capacity = analytic_model->calculate_ion_heat_capacity(T, rho);

  Ensure(heat_capacity >= 0.0);
  return heat_capacity;
}

//------------------------------------------------------------------------------------------------//
/*!
 *
 * \brief Return a field of ion heat capacities given fields of temperature
 * and density.
 *
 * Given temperature and density fields, return an ion heat capacity field.
 * The ion heat capacity field is the same size as the input temperature and
 * density fields. The ion heat capacity is defined by the analytic model
 * given to the constructor (Analytic_EoS()).
 *
 * \param T material temperature field in keV
 * \param rho density field in g/cm^3
 * \return ion heat capacity field in kJ/g/keV
 */
Analytic_EoS::sf_double Analytic_EoS::getIonHeatCapacity(const sf_double &T,
                                                         const sf_double &rho) const {
  Require(T.size() == rho.size());

  // define the return ion heat capacity field
  sf_double heat_capacity(T.size(), 0.0);

  // loop through T/rho field and solve for the heat capacity
  for (size_t i = 0; i < T.size(); i++) {
    Check(T[i] >= 0.0);
    Check(rho[i] >= 0.0);

    heat_capacity[i] = analytic_model->calculate_ion_heat_capacity(T[i], rho[i]);

    Check(heat_capacity[i] >= 0.0);
  }

  return heat_capacity;
}

//------------------------------------------------------------------------------------------------//
/*!
 *
 * \brief Return a scalar number of free electrons per ion given a scalar
 * temperature and density.
 *
 * Given a scalar temperature and density, return the number of free
 * electrons per ion.  The number of free electrons per ion is defined by the
 * analytic model given to the constructor (Analytic_EoS()).
 *
 * \param T material temperature in keV
 * \param rho density in g/cm^3
 * \return number of free electrons per ion
 */
double Analytic_EoS::getNumFreeElectronsPerIon(double T, double rho) const {
  Require(T >= 0.0);
  Require(rho >= 0.0);

  double number = analytic_model->calculate_num_free_elec_per_ion(T, rho);

  Ensure(number >= 0.0);
  return number;
}

//------------------------------------------------------------------------------------------------//
/*!
 *
 * \brief Return a field of numbers of free electrons per ion given fields of
 * temperature and density.
 *
 * Given temperature and density fields, return a number of free electrons
 * per ion field.  The number of free electrons per ion field is the same
 * size as the input temperature and density fields. The number of free
 * electrons per ion is defined by the analytic model given to the
 * constructor (Analytic_EoS()).
 *
 * \param T material temperature field in keV
 * \param rho density field in g/cm^3
 * \return number of free electrons per ion
 */
Analytic_EoS::sf_double Analytic_EoS::getNumFreeElectronsPerIon(const sf_double &T,
                                                                const sf_double &rho) const {
  Require(T.size() == rho.size());

  // define the return ion heat capacity field
  sf_double number(T.size(), 0.0);

  // loop through T/rho field and solve for the number of free electrons
  // per ion
  for (size_t i = 0; i < T.size(); i++) {
    Check(T[i] >= 0.0);
    Check(rho[i] >= 0.0);

    number[i] = analytic_model->calculate_num_free_elec_per_ion(T[i], rho[i]);

    Check(number[i] >= 0.0);
  }

  return number;
}

//------------------------------------------------------------------------------------------------//
/*!
 *
 * \brief Return a scalar electron thermal conductivity given a scalar
 * temperature and density.
 *
 * Given a scalar temperature and density, return the electron themal
 * conductivity.  The electron thermal conductivity is defined by the
 * analytic model given to the constructor (Analytic_EoS()).
 *
 * \param T material temperature in keV
 * \param rho density in g/cm^3
 * \return electron thermal conductivity in /s/cm
 */
double Analytic_EoS::getElectronThermalConductivity(double T, double rho) const {
  Require(T >= 0.0);
  Require(rho >= 0.0);

  double thermal_conductivity = analytic_model->calculate_elec_thermal_conductivity(T, rho);

  Ensure(thermal_conductivity >= 0.0);
  return thermal_conductivity;
}

//------------------------------------------------------------------------------------------------//
/*!
 *
 * \brief Return a field of electron thermal conductivities given fields of
 * temperature and density.
 *
 * Given temperature and density fields, return an electron thermal
 * conductivity field.  The electron thermal conductivity field is the same
 * size as the input temperature and density fields. The electron thermal
 * conductivity is defined by the analytic model given to the constructor
 * (Analytic_EoS()).
 *
 * \param T material temperature field in keV
 * \param rho density field in g/cm^3
 * \return electron heat capacity field in /s/cm
 */
Analytic_EoS::sf_double Analytic_EoS::getElectronThermalConductivity(const sf_double &T,
                                                                     const sf_double &rho) const {
  Require(T.size() == rho.size());

  // define the return electron thermal conductivity field
  sf_double thermal_conductivity(T.size(), 0.0);

  // loop through T/rho field and solve for the thermal conductivity
  for (size_t i = 0; i < T.size(); i++) {
    Check(T[i] >= 0.0);
    Check(rho[i] >= 0.0);

    thermal_conductivity[i] = analytic_model->calculate_elec_thermal_conductivity(T[i], rho[i]);

    Check(thermal_conductivity[i] >= 0.0);
  }

  return thermal_conductivity;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the electron temperature given density and the specific
 *        electron energy.
 *
 * Given density and specific electron energy fields, return an electron
 * temperature.  The electron thermal conductivity is defined by the analytic
 * model given to the constructor (Analytic_EoS()).
 *
 * \param rho density field in g/cm^3 (unused by the current analytic models)
 * \param Ue specific electron energy (kJ/g)
 * \param Tguess An electron temperature (keV) guess that can be used to help
 *        the root finder.
 * \return electron temperature in keV.
 */
double Analytic_EoS::getElectronTemperature(double rho, double Ue, double Tguess) const {
  Require(Tguess >= 0.0);
  Require(Ue >= 0.0);
  Require(rho > 0.0);

  double T_new = analytic_model->calculate_elec_temperature(rho, Ue, Tguess);

  Ensure(T_new >= 0.0);
  return T_new;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the ion temperature given density and the specific ion
 *        energy.
 *
 * Given density and specific ion energy fields, return an ion temperature.
 * The ion thermal conductivity is defined by the analytic model given to
 * the constructor (Analytic_EoS()).
 *
 * \param rho density field in g/cm^3 (unused by the current analytic models)
 * \param Uic specific ion energy plus cold curve energy density (kJ/g)
 * \param Tguess An ion temperature (keV) guess that can be used to help the
 *        root finder.
 * \return ion temperature in keV.
 */
double Analytic_EoS::getIonTemperature(double rho, double Uic, double Tguess) const {
  Require(Tguess >= 0.0);
  Require(Uic >= 0.0);
  Require(rho > 0.0);

  double T_new = analytic_model->calculate_ion_temperature(rho, Uic, Tguess);

  Ensure(T_new >= 0.0);
  return T_new;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Pack an analytic EoS opacity.
 *
 * This function will pack up the Analytic_EoS into a char
 * array (represented by a vector<char>).  The Analytic_EoS_Model derived
 * class must have a pack function; this is enforced by the virtual
 * Analytic_EoS_Model base class.
 */
Analytic_EoS::sf_char Analytic_EoS::pack() const {
  Require(analytic_model);

  // make a packer
  rtt_dsxx::Packer packer;

  // first pack up the analytic model
  sf_char anal_model = analytic_model->pack();

  // now add up the total size (in bytes): size of analytic model + 1
  // int for size of analytic model
  Check(anal_model.size() + 1 * sizeof(int) < INT_MAX);
  int size = static_cast<int>(anal_model.size() + 1 * sizeof(int));

  // make a char array
  sf_char packed(size);

  // set the buffer
  packer.set_buffer(size, &packed[0]);

  // pack the anal_model size
  packer << static_cast<int>(anal_model.size());

  // now pack the anal model
  for (auto &am : anal_model)
    packer << am;

  Ensure(packer.get_ptr() == &packed[0] + size);
  return packed;
}

} // end namespace rtt_cdi_analytic

//------------------------------------------------------------------------------------------------//
// end of Analytic_EoS.cc
//------------------------------------------------------------------------------------------------//
