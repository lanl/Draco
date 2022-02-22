//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_analytic/Analytic_EICoupling.cc
 * \author Mathew Cleveland
 * \date   March 2019
 * \brief  Analytic_EICoupling member definitions.
 * \note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Analytic_EICoupling.hh"
#include "ds++/Packing_Utils.hh"

namespace rtt_cdi_analytic {

//------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for an analytic EICoupling model.
 *
 * This constructor builds an analytic EICoupling model defined by the
 * rtt_cdi_analytic::Analytic_EICoupling_Model derived class argument.
 *
 * \param model_in shared_ptr to a derived
 *        rtt_cdi_analytic::Analytic_EICoupling_Model object
 */
Analytic_EICoupling::Analytic_EICoupling(SP_Analytic_Model model_in)
    : analytic_model(std::move(model_in)) {
  Ensure(analytic_model);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Unpacking constructor.
 *
 * This constructor rebuilds and Analytic_EICoupling from a vector<char> that
 * was created by a call to pack().  It can only rebuild Analytic_Model types
 * that have been registered in the rtt_cdi_analytic::EICoupling_Models
 * enumeration.
 */
Analytic_EICoupling::Analytic_EICoupling(const sf_char &packed) : analytic_model() {
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
  if (indicator == CONSTANT_ANALYTIC_EICOUPLING_MODEL) {
    analytic_model = std::make_shared<Constant_Analytic_EICoupling_Model>(packed_analytic);
  } else {
    Insist(0, "Unregistered analytic EICoupling model!");
  }

  Ensure(analytic_model);
}

//------------------------------------------------------------------------------------------------//
// EICoupling INTERFACE FUNCTIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief EICoupling accessor that returns a single electron-ion coupling
 *        coefficient.
 *
 * \param[in] eTemperature The electron temperature value for which an
 *     opacity value is being requested (keV).
 * \param[in] iTemperature The electron temperature value for which an
 *     opacity value is being requested (keV).
 * \param[in] density The density value for which an opacity
 *     value is being requested (g/cm^3).
 * \param[in] w_e is the plasma electron frequency (as defined by Eq. 3.41 in
 *     Brown, Preston, and Singleton, 'Physics Reports', V410, Issue 4, 2005)
 * \param[in] w_i is the average plasma ion frequency (as defined by Eq. 3.61 in
 *     Brown, Preston, and Singleton, 'Physics Reports', V410, Issue 4, 2005)
 * \return A electron-ion coupling coeffiecent (1/s).
 */
double Analytic_EICoupling::getElectronIonCoupling(const double eTemperature,
                                                   const double iTemperature, const double density,
                                                   const double w_e, const double w_i) const {
  Require(eTemperature >= 0.0);
  Require(iTemperature >= 0.0);
  Require(density >= 0.0);
  Require(w_e >= 0.0);
  Require(w_i >= 0.0);

  double ei_coupling =
      analytic_model->calculate_ei_coupling(eTemperature, iTemperature, density, w_e, w_i);

  Ensure(ei_coupling >= 0.0);
  return ei_coupling;
}

/*!
 * \brief EICoupling accessor that returns a vector electron-ion coupling
 * coefficients.
 *
 * \param[in] vetemperature The electron temperature value for which an
 *     opacity value is being requested (keV).
 * \param[in] vitemperature The electron temperature value for which an
 *     opacity value is being requested (keV).
 * \param[in] vdensity The density value for which an opacity
 *     value is being requested (g/cm^3).
 * \param[in] vw_e is the plasma electron frequency (as defined by Eq. 3.41 in
 *     Brown, Preston, and Singleton, 'Physics Reports', V410, Issue 4, 2005)
 * \param[in] vw_i is the average plasma ion frequency (as defined by Eq. 3.61 in
 *     Brown, Preston, and Singleton, 'Physics Reports', V410, Issue 4, 2005)
 * \return A vector of electron-ion coupling coeffiecent (kJ/cc/keV/s).

 */
Analytic_EICoupling::sf_double Analytic_EICoupling::getElectronIonCoupling(
    const std::vector<double> &vetemperature, const std::vector<double> &vitemperature,
    const std::vector<double> &vdensity, const std::vector<double> &vw_e,
    const std::vector<double> &vw_i) const {
  Require(vetemperature.size() == vdensity.size());
  Require(vitemperature.size() == vdensity.size());
  Require(vw_e.size() == vdensity.size());
  Require(vw_i.size() == vdensity.size());

  // define the return electron-ion coupling field
  sf_double ei_coupling(vetemperature.size(), 0.0);

  // loop through
  for (size_t i = 0; i < vetemperature.size(); i++) {
    Check(vetemperature[i] >= 0.0);
    Check(vitemperature[i] >= 0.0);
    Check(vdensity[i] >= 0.0);
    Check(vw_e[i] >= 0.0);
    Check(vw_i[i] >= 0.0);

    ei_coupling[i] = analytic_model->calculate_ei_coupling(vetemperature[i], vitemperature[i],
                                                           vdensity[i], vw_e[i], vw_i[i]);

    Check(ei_coupling[i] >= 0.0);
  }

  return ei_coupling;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Pack an analytic EICoupling.
 *
 * This function will pack up the Analytic_EICoupling into a char
 * array (represented by a vector<char>).  The Analytic_EICoupling_Model derived
 * class must have a pack function; this is enforced by the virtual
 * Analytic_EICoupling_Model base class.
 */
Analytic_EICoupling::sf_char Analytic_EICoupling::pack() const {
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
// end of Analytic_EICoupling.cc
//------------------------------------------------------------------------------------------------//
