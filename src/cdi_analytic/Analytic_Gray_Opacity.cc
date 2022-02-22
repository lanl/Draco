//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_analytic/Analytic_Gray_Opacity.cc
 * \author Thomas M. Evans
 * \date   Fri Aug 24 13:13:46 2001
 * \brief  Analytic_Gray_Opacity member definitions.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Analytic_Gray_Opacity.hh"
#include "ds++/Packing_Utils.hh"

namespace rtt_cdi_analytic {

//------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for an analytic gray opacity model.
 *
 * This constructor builds an opacity model defined by the
 * rtt_cdi_analytic::Analytic_Opacity_Model derived class argument.
 *
 * The reaction type for this instance of the class is determined by the
 * rtt_cdi::Reaction argument.
 *
 * \param model_in shared_ptr to a derived
 *                 rtt_cdi_analytic::Analytic_Opacity_Model object
 * \param reaction_in rtt_cdi::Reaction type (enumeration)
 * \param cdi_model_in CDI model type
 */
Analytic_Gray_Opacity::Analytic_Gray_Opacity(SP_Analytic_Model model_in,
                                             rtt_cdi::Reaction reaction_in,
                                             rtt_cdi::Model cdi_model_in)
    : analytic_model(std::move(model_in)), reaction(reaction_in), model(cdi_model_in) {
  Require(reaction == rtt_cdi::TOTAL || reaction == rtt_cdi::ABSORPTION ||
          reaction == rtt_cdi::SCATTERING);

  Ensure(analytic_model);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Unpacking constructor.
 *
 * This constructor rebuilds and Analytic_Gray_Opacity from a vector<char>
 * that was created by a call to pack().  It can only rebuild Analytic_Model
 * types that have been registered in the rtt_cdi_analytic::Opacity_Models
 * enumeration.
 */
Analytic_Gray_Opacity::Analytic_Gray_Opacity(const sf_char &packed)
    : analytic_model(), reaction(), model() {
  // the packed size must be at least 4 integers (size, reaction type,
  // model type, analytic model indicator)
  Require(packed.size() >= 4 * sizeof(int));

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
  for (int i = 0; i < size_analytic; i++)
    unpacker >> packed_analytic[i];

  // unpack the reaction and model type
  int react_int, model_int;
  unpacker >> react_int >> model_int;
  Check(unpacker.get_ptr() == &packed[0] + packed.size());

  // assign the reaction and type
  reaction = static_cast<rtt_cdi::Reaction>(react_int);
  model = static_cast<rtt_cdi::Model>(model_int);
  Check(reaction == rtt_cdi::ABSORPTION || reaction == rtt_cdi::SCATTERING ||
        reaction == rtt_cdi::TOTAL);

  // now reset the buffer so that we can determine the analytic model
  // indicator
  unpacker.set_buffer(size_analytic, &packed_analytic[0]);

  // unpack the indicator
  int indicator;
  unpacker >> indicator;

  // now determine which analytic model we need to build
  if (indicator == CONSTANT_ANALYTIC_OPACITY_MODEL) {
    analytic_model = std::make_shared<Constant_Analytic_Opacity_Model>(packed_analytic);
  } else if (indicator == POLYNOMIAL_ANALYTIC_OPACITY_MODEL) {
    analytic_model = std::make_shared<Polynomial_Analytic_Opacity_Model>(packed_analytic);
  } else {
    Insist(0, "Unregistered analytic opacity model!");
  }

  Ensure(analytic_model);
}

//------------------------------------------------------------------------------------------------//
// OPACITY INTERFACE FUNCTIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return a scalar opacity given a scalar temperature and density.
 *
 * Given a scalar temperature and density, return an opacity for the reaction
 * type specified by the constructor.  The analytic opacity model is
 * specified in the constructor (Analytic_Gray_Opacity()).
 *
 * \param temperature material temperature in keV
 * \param density material density in g/cm^3
 * \return opacity (coefficient) in cm^2/g
 *
 */
double Analytic_Gray_Opacity::getOpacity(double temperature, double density) const {
  Require(temperature >= 0.0);
  Require(density >= 0.0);

  // define return opacity
  double opacity = analytic_model->calculate_opacity(temperature, density);

  Ensure(opacity >= 0.0);
  return opacity;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return a field of opacities given a field of temperatures and a
 * scalar density.
 *
 * Given a field of temperatures and a scalar density, return an opacity
 * field for the reaction type specified by the constructor.  The analytic
 * opacity model is specified in the constructor
 * (Analytic_Gray_Opacity()).  The returned opacity field has the same number
 * of elements as the temperature field.
 *
 * The field type for temperatures is an std::vector.
 *
 * \param temperature std::vector of material temperatures in keV
 * \param density material density in g/cm^3
 * \return std::vector of opacities (coefficients) in cm^2/g
 */
Analytic_Gray_Opacity::sf_double Analytic_Gray_Opacity::getOpacity(const sf_double &temperature,
                                                                   double density) const {
  Require(density >= 0.0);

  // define the return opacity field (same size as temperature field)
  sf_double opacity(temperature.size(), 0.0);

  // define an opacity iterator
  auto sig = opacity.begin();

  // loop through temperatures and solve for opacity
  for (auto T = temperature.begin(); T != temperature.end(); T++, sig++) {
    Check(*T >= 0.0);

    // define opacity
    *sig = analytic_model->calculate_opacity(*T, density);

    Check(*sig >= 0.0);
  }

  return opacity;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return a field of opacities given a field of densities and a scalar
 * temperature.
 *
 * Given a field of densities and a scalar temperature, return an opacity
 * field for the reaction type specified by the constructor.  The analytic
 * opacity model is specified in the constructor
 * (Analytic_Gray_Opacity()).  The returned opacity field has the same number
 * of elements as the density field.
 *
 * The field type for densities is an std::vector.
 *
 * \param temperature material temperature in keV
 * \param density std::vector of material densities in g/cc
 * \return std::vector of opacities (coefficients) in cm^2/g
 */
Analytic_Gray_Opacity::sf_double Analytic_Gray_Opacity::getOpacity(double temperature,
                                                                   const sf_double &density) const {
  Require(temperature >= 0.0);

  // define the return opacity field (same size as density field)
  sf_double opacity(density.size(), 0.0);

  // define an opacity iterator
  auto sig = opacity.begin();

  // loop through densities and solve for opacity
  for (auto rho = density.begin(); rho != density.end(); rho++, sig++) {
    Check(*rho >= 0.0);

    // define opacity
    *sig = analytic_model->calculate_opacity(temperature, *rho);

    Check(*sig >= 0.0);
  }

  return opacity;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Pack an analytic gray opacity.
 *
 * This function will pack up the Analytic_Gray_Opacity into a char array
 * (represented by a vector<char>).  The Analytic_Opacity_Model derived class
 * must have a pack function; this is enforced by the virtual
 * Analytic_Opacity_Model base class.
 */
Analytic_Gray_Opacity::sf_char Analytic_Gray_Opacity::pack() const {
  Require(analytic_model);

  // make a packer
  rtt_dsxx::Packer packer;

  // first pack up the analytic model
  sf_char anal_model = analytic_model->pack();

  // now add up the total size (in bytes): size of analytic model + 3
  // int--one for reaction type, one for model type, and one for size of
  // analytic model
  Check(anal_model.size() + 3 * sizeof(int) < INT_MAX);
  int size = static_cast<int>(anal_model.size() + 3 * sizeof(int));

  // make a char array
  sf_char packed(size);

  // set the buffer
  packer.set_buffer(size, &packed[0]);

  // pack the anal_model size
  packer << static_cast<int>(anal_model.size());

  // now pack the anal model
  for (auto &am : anal_model)
    packer << am;

  // pack the reaction and model type
  packer << static_cast<int>(reaction) << static_cast<int>(model);

  Ensure(packer.get_ptr() == &packed[0] + size);
  return packed;
}

} // end namespace rtt_cdi_analytic

//------------------------------------------------------------------------------------------------//
// end of Analytic_Gray_Opacity.cc
//------------------------------------------------------------------------------------------------//
