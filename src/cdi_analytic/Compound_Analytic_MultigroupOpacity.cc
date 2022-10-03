//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_analytic/Compound_Analytic_MultigroupOpacity.cc
 * \brief  Compound_Analytic_MultigroupOpacity class member definitions.
 * \note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Compound_Analytic_MultigroupOpacity.hh"
#include "ds++/Packing_Utils.hh"
#include "ds++/dbc.hh"

namespace rtt_cdi_analytic {

//------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for an analytic multigroup opacity model.
 *
 * This constructor builds an opacity model defined by the
 * rtt_cdi_analytic::Analytic_Opacity_Model derived class argument.
 *
 * The reaction type for this instance of the class is determined by the
 * rtt_cdi::Reaction argument.
 *
 * The group structure (in keV) must be provided by the groups argument.  The
 * number of Analytic_Opacity_Model objects given in the models argument must
 * be equal to the number of groups.
 *
 * \param groups vector containing the group boundaries in keV from lowest to
 *        highest
 * \param models vector containing shared_ptrs to Analytic_Model derived types
 *        for each group, the size should be groups.size() - 1
 * \param reaction_in rtt_cdi::Reaction type (enumeration)
 * \param model_in Enum specifying CDI model.
 */
Compound_Analytic_MultigroupOpacity::Compound_Analytic_MultigroupOpacity(
    const sf_double &groups, const sf_Analytic_Model &models, rtt_cdi::Reaction reaction_in,
    rtt_cdi::Model model_in)
    : Analytic_MultigroupOpacity(groups, reaction_in, model_in), group_models(models) {
  Require(groups.size() - 1 == models.size());
  Require(rtt_dsxx::is_strict_monotonic_increasing(groups.begin(), groups.end()));

  Ensure(check_class_invariant());
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Unpacking constructor.
 *
 * This constructor rebuilds and Compound_Analytic_MultigroupOpacity from a
 * vector<char> that was created by a call to pack().  It can only rebuild
 * Analytic_Model types that have been registered in the
 * rtt_cdi_analytic::Opacity_Models enumeration.
 */
Compound_Analytic_MultigroupOpacity::Compound_Analytic_MultigroupOpacity(const sf_char &packed)
    : Analytic_MultigroupOpacity(packed), group_models() {
  // get the number of group boundaries
  size_t const num_groups = getGroupBoundaries().size() - 1;
  unsigned const base_size = Analytic_MultigroupOpacity::packed_size();

  // make an unpacker
  rtt_dsxx::Unpacker unpacker;

  // register the unpacker
  unpacker.set_buffer(packed.size() - base_size, &packed[base_size]);

  // make the model vectors
  group_models.resize(num_groups);

  // unpack the models
  std::vector<sf_char> models(num_groups);
  int model_size = 0;
  for (auto &model_1dvec : models) {
    // unpack the size of the analytic model
    unpacker >> model_size;
    Check(static_cast<size_t>(model_size) >= sizeof(int));

    model_1dvec.resize(model_size);

    // unpack the model
    for (auto &elem : model_1dvec)
      unpacker >> elem;
  }

  // now rebuild the analytic models
  int indicator = 0;
  for (size_t i = 0; i < models.size(); ++i) {
    // reset the buffer
    unpacker.set_buffer(models[i].size(), &models[i][0]);

    // get the indicator for this model (first packed datum)
    unpacker >> indicator;

    // now determine which analytic model we need to build
    if (indicator == CONSTANT_ANALYTIC_OPACITY_MODEL) {
      group_models[i] = std::make_shared<Constant_Analytic_Opacity_Model>(models[i]);
    } else if (indicator == POLYNOMIAL_ANALYTIC_OPACITY_MODEL) {
      group_models[i] = std::make_shared<Polynomial_Analytic_Opacity_Model>(models[i]);
    } else {
      Insist(0, "Unregistered analytic opacity model!");
    }

    Ensure(group_models[i]);
  }

  Ensure(check_class_invariant());
}
//------------------------------------------------------------------------------------------------//
bool Compound_Analytic_MultigroupOpacity::check_class_invariant() const {
  return group_models.size() + 1 == getGroupBoundaries().size();
}

//------------------------------------------------------------------------------------------------//
// OPACITY INTERFACE FUNCTIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return the group opacities given a scalar temperature and density.
 *
 * Given a scalar temperature and density, return the group opacities
 * (vector<double>) for the reaction type specified by the constructor.  The
 * analytic opacity model is specified in the constructor
 * (Compound_Analytic_MultigroupOpacity()).
 *
 * \param temperature material temperature in keV
 * \param density material density in g/cm^3
 * \return group opacities (coefficients) in cm^2/g
 *
 */
Compound_Analytic_MultigroupOpacity::sf_double
Compound_Analytic_MultigroupOpacity::getOpacity(double temperature, double density) const {
  Require(temperature >= 0.0);
  Require(density >= 0.0);

  // return opacities
  sf_double opacities(group_models.size(), 0.0);

  // loop through groups and get opacities
  auto nu = getGroupBoundaries();
  Check(nu.size() == group_models.size() + 1);

  for (size_t i = 0; i < opacities.size(); ++i) {
    Check(group_models[i]);

    // assign the opacity based on the group model
    opacities[i] = group_models[i]->calculate_opacity(temperature, density, nu[i], nu[i + 1]);

    Check(opacities[i] >= 0.0);
  }

  return opacities;
}

//------------------------------------------------------------------------------------------------//
/*!
 *
 * \brief Return a vector of multigroup opacities given a vector of
 * temperatures and a scalar density.
 *
 * Given a field of temperatures and a scalar density, return a vector of
 * multigroup opacities (vector<vector<double>>) for the reaction type
 * specified by the constructor.  The analytic opacity model is specified in
 * the constructor (Analytic_Gray_Opacity()).  The returned opacity field is
 * indexed [num_temperatures][num_groups].
 *
 * The field type for temperatures is an std::vector.
 *
 * \param temperature std::vector of material temperatures in keV
 *
 * \param density material density in g/cm^3
 *
 * \return std::vector<std::vector> of multigroup opacities (coefficients) in
 * cm^2/g indexed [temperature][group]
 */
Compound_Analytic_MultigroupOpacity::vf_double
Compound_Analytic_MultigroupOpacity::getOpacity(const sf_double &temperature,
                                                double density) const {
  Require(density >= 0.0);

  // define the return opacity field (same size as temperature field); each
  // entry has a vector sized by the number of groups
  vf_double opacities(temperature.size(), sf_double(group_models.size(), 0.0));

  // loop through temperatures and solve for opacity
  for (size_t i = 0; i < opacities.size(); ++i) {
    Check(temperature[i] >= 0.0);

    // loop through groups
    for (size_t j = 0; j < opacities[i].size(); ++j) {
      Check(group_models[j]);

      // assign the opacity based on the group model
      opacities[i][j] = group_models[j]->calculate_opacity(temperature[i], density);

      Check(opacities[i][j] >= 0.0);
    }
  }

  return opacities;
}

//------------------------------------------------------------------------------------------------//
/*!
 *
 * \brief Return a vector of multigroup opacities given a vector of
 * density and a scalar temperature.
 *
 * Given a field of densities and a scalar temperature, return a vector of
 * multigroup opacities (vector<vector<double>>) for the reaction type
 * specified by the constructor.  The analytic opacity model is specified in
 * the constructor (Analytic_Gray_Opacity()).  The returned opacity field is
 * indexed [num_density][num_groups].
 *
 * The field type for density is an std::vector.
 *
 * \param temperature in keV
 *
 * \param density std::vector of material densities in g/cm^3
 *
 * \return std::vector<std::vector> of multigroup opacities (coefficients) in
 * cm^2/g indexed [density][group]
 */
Compound_Analytic_MultigroupOpacity::vf_double
Compound_Analytic_MultigroupOpacity::getOpacity(double temperature,
                                                const sf_double &density) const {
  Require(temperature >= 0.0);

  // define the return opacity field (same size as density field); each
  // entry has a vector sized by the number of groups
  vf_double opacities(density.size(), sf_double(group_models.size(), 0.0));

  // loop through densities and solve for opacity
  for (size_t i = 0; i < opacities.size(); ++i) {
    Check(density[i] >= 0.0);

    // loop through groups
    for (size_t j = 0; j < opacities[i].size(); ++j) {
      Check(group_models[j]);

      // assign the opacity based on the group model
      opacities[i][j] = group_models[j]->calculate_opacity(temperature, density[i]);

      Check(opacities[i][j] >= 0.0);
    }
  }

  return opacities;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Pack an analytic multigroup opacity.
 *
 * This function will pack up the Analytic_Mulitgroup_Opacity into a char
 * array (represented by a vector<char>).  The Analytic_Opacity_Model derived
 * class must have a pack function; this is enforced by the virtual
 * Analytic_Opacity_Model base class.
 */
Compound_Analytic_MultigroupOpacity::sf_char Compound_Analytic_MultigroupOpacity::pack() const {
  // make a packer
  rtt_dsxx::Packer packer;

  // make a char array
  sf_char packed = Analytic_MultigroupOpacity::pack();

  // first pack up models
  std::vector<sf_char> models(group_models.size());
  int num_bytes_models = 0;

  // loop through and pack up the models
  for (size_t i = 0; i < models.size(); ++i) {
    Check(group_models[i]);

    // (KPL) The original logic in this loop was simply:
    // models[i] = group_models[i]->pack();
    // However, that code causes a link-time failure for xl-16.1.1+gcc-9.3 DEBUG builds.
    // The following logic (create temporary, resize destination char vector, assign destination
    // to temporary) appears to work for all supported systems+compilers as of October 2022.
    sf_char this_pack = group_models[i]->pack();
    models[i].resize(this_pack.size());
    models[i] = this_pack;
    Check(num_bytes_models + models[i].size() < INT_MAX);
    num_bytes_models += static_cast<int>(models[i].size());
  }

  // now add up the total size; number of groups + 1 int for number of
  // groups, number of models + size in each model + models, 1 int for
  // reaction type, 1 int for model type
  size_t const base_size = packed.size();
  size_t const size = models.size() * sizeof(int) + num_bytes_models;

  // extend the char array
  packed.resize(size + base_size);

  // set the buffer
  packer.set_buffer(size, &packed[base_size]);

  // pack each models size and data
  for (auto &model_1dvec : models) {
    // pack the size of this model
    packer << static_cast<int>(model_1dvec.size());

    // now pack the model data
    for (auto &elem : model_1dvec)
      packer << elem;
  }

  Ensure(packer.get_ptr() == &packed[base_size] + size);
  return packed;
}

} // end namespace rtt_cdi_analytic

//------------------------------------------------------------------------------------------------//
// end of Compound_Analytic_MultigroupOpacity.cc
//------------------------------------------------------------------------------------------------//
