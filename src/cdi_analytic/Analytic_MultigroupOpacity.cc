//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_analytic/Analytic_MultigroupOpacity.cc
 * \author Thomas M. Evans
 * \date   Tue Nov 13 11:19:59 2001
 * \brief  Analytic_MultigroupOpacity class member definitions.
 * \note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Analytic_MultigroupOpacity.hh"
#include "ds++/Packing_Utils.hh"

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
 * \param[in] groups vector containing the group boundaries in keV from lowest
 *               to highest
 * \param[in] reaction_in rtt_cdi::Reaction type (enumeration)
 * \param[in] model_in CDI type
 */
Analytic_MultigroupOpacity::Analytic_MultigroupOpacity(sf_double groups,
                                                       rtt_cdi::Reaction reaction_in,
                                                       rtt_cdi::Model model_in)
    : group_boundaries(std::move(groups)), reaction(reaction_in), model(model_in) {
  Require(reaction == rtt_cdi::TOTAL || reaction == rtt_cdi::ABSORPTION ||
          reaction == rtt_cdi::SCATTERING);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Unpacking constructor.
 *
 * This constructor rebuilds and Analytic_MultigroupOpacity from a vector<char>
 * that was created by a call to pack().  It can only rebuild Analytic_Model
 * types that have been registered in the rtt_cdi_analytic::Opacity_Models
 * enumeration.
 */
Analytic_MultigroupOpacity::Analytic_MultigroupOpacity(const sf_char &packed)
    : group_boundaries(), reaction(), model() {
  // the packed size must be at least 4 integers (number of groups, reaction
  // type, model type, analytic model indicator)
  Require(packed.size() >= 4 * sizeof(int));

  // make an unpacker
  rtt_dsxx::Unpacker unpacker;

  // register the unpacker
  unpacker.set_buffer(packed.size(), &packed[0]);

  // unpack the number of group boundaries
  int ngrp_bounds = 0;
  unpacker >> ngrp_bounds;

  // make the group boundaries and model vectors
  group_boundaries.resize(ngrp_bounds);

  // unpack the group boundaries
  for (int i = 0; i < ngrp_bounds; i++)
    unpacker >> group_boundaries[i];

  // unpack the reaction and model type
  int reaction_int, model_int;
  unpacker >> reaction_int >> model_int;

  // assign the reaction and model type
  reaction = static_cast<rtt_cdi::Reaction>(reaction_int);
  model = static_cast<rtt_cdi::Model>(model_int);
  Check(reaction == rtt_cdi::ABSORPTION || reaction == rtt_cdi::SCATTERING ||
        reaction == rtt_cdi::TOTAL);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Pack an analytic multigroup opacity.
 *
 * This function will pack up the Analytic_Mulitgroup_Opacity into a char array
 * (represented by a vector<char>).  The Analytic_Opacity_Model derived class
 * must have a pack function; this is enforced by the virtual
 * Analytic_Opacity_Model base class.
 */
Analytic_MultigroupOpacity::sf_char Analytic_MultigroupOpacity::pack() const {
  // make a packer
  rtt_dsxx::Packer packer;

  // now add up the total size; number of groups + 1 int for number of groups,
  // number of models + size in each model + models, 1 int for reaction type,
  // 1 int for model type
  Check(3 * sizeof(int) + group_boundaries.size() * sizeof(double) < INT_MAX);
  int size = static_cast<int>(3 * sizeof(int) + group_boundaries.size() * sizeof(double));

  // make a char array
  sf_char packed(size);

  // set the buffer
  packer.set_buffer(size, &packed[0]);

  // pack the number of groups and group boundaries
  packer << static_cast<int>(group_boundaries.size());
  for (auto &gb : group_boundaries)
    packer << gb;

  // now pack the reaction and model type
  packer << static_cast<int>(reaction) << static_cast<int>(model);

  Ensure(packer.get_ptr() == &packed[0] + size);
  return packed;
}

//------------------------------------------------------------------------------------------------//
unsigned Analytic_MultigroupOpacity::packed_size() const {
  // This must match the size calculated in the previous function
  Check(3 * sizeof(int) + group_boundaries.size() * sizeof(double) < UINT_MAX);
  return static_cast<unsigned>(3 * sizeof(int) + group_boundaries.size() * sizeof(double));
}
} // end namespace rtt_cdi_analytic

//------------------------------------------------------------------------------------------------//
// end of Analytic_MultigroupOpacity.cc
//------------------------------------------------------------------------------------------------//
