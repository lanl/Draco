//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_analytic/Compound_Analytic_MultigroupOpacity.hh
 * \author Thomas M. Evans
 * \date   Tue Nov 13 11:19:59 2001
 * \brief  Compound_Analytic_MultigroupOpacity class definition.
 * \note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_analytic_Compound_Analytic_MultigroupOpacity_hh
#define rtt_cdi_analytic_Compound_Analytic_MultigroupOpacity_hh

#include "Analytic_MultigroupOpacity.hh"

namespace rtt_cdi_analytic {

//================================================================================================//
/*!
 * \class Compound_Analytic_MultigroupOpacity
 *
 * \brief Derived Analytic_MultigroupOpacity class for analytic opacities.
 *
 * The Compound_Analytic_MultigroupOpacity class is a derived
 * Analytic_MultigroupOpacity class.  It provides analytic opacity data. The
 * specific analytic opacity model is derived from the
 * rtt_cdi_analytic::Analytic_Opacity_Model base class.  Several pre-built
 * derived classes are provided in Analytic_Models.hh.
 *
 * Clients of this class can provide any analytic model class as long as it
 * conforms to the rtt_cdi_analytic::Analytic_Opacity_Model interface.  This
 * interface consists of a single function,
 * Analytic_Opacity_Model::calculate_opacity().
 *
 * An rtt_cdi_analytic::Analytic_Opacity_Model is provided for each analytic
 * energy group.  In other words, the user specifies the group structure for
 * the multigroup data and applies a Analytic_Model for each group.  This way
 * a different rtt_cdi_analytic::Analytic_Model can be used in each group.
 * For example, the client could choose to use a
 * rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model for the low frequency
 * groups and a rtt_cdi_analytic::Constant_Analytic_Opacity_Model for the high
 * frequency groups.
 *
 * Note that opacities are returned in units of cm^2/g. Thus the resultant
 * opacity must be multiplied by density to get units of 1/cm.  See the
 * documentation in rtt_cdi_analytic::Analytic_Model for more info.
 *
 * The constructors take a rtt_cdi::Reaction argument to determine the
 * reaction type.  The enumeration rtt_cdi::Reaction can have the value
 * TOTAL, ABSORPTION, or SCATTERING.
 *
 * The default rtt_cdi::Model for an Analytic_Gray_Opacity is
 * rtt_cdi::ANALYTIC.  However, this can be overridden in the constructor.
 *
 * This class conforms to the interface specified by
 * Analytic_MultigroupOpacity and can be used with rtt_cdi::CDI to get
 * analytic opacities.
 *
 * \example cdi_analytic/test/tstCompound_Analytic_MultigroupOpacity.cc
 *
 * Example usage of Compound_Analytic_MultigroupOpacity, Analytic_Opacity_Model,
 * and their incorporation into rtt_cdi::CDI.
 */
//================================================================================================//

class Compound_Analytic_MultigroupOpacity : public Analytic_MultigroupOpacity {
public:
  // Useful typedefs.
  using SP_Analytic_Model = std::shared_ptr<Analytic_Opacity_Model>;
  using const_Model = std::shared_ptr<const Analytic_Opacity_Model>;
  using sf_Analytic_Model = std::vector<SP_Analytic_Model>;
  using sf_double = std::vector<double>;
  using vf_double = std::vector<sf_double>;
  using std_string = std::string;
  using sf_char = std::vector<char>;

private:
  // Analytic models for each group.
  sf_Analytic_Model group_models;

public:
  // Constructor.
  Compound_Analytic_MultigroupOpacity(const sf_double &groups, const sf_Analytic_Model &models,
                                      rtt_cdi::Reaction reaction_in,
                                      rtt_cdi::Model model_in = rtt_cdi::ANALYTIC);

  // Constructor for packed Analytic_Multigroup_Opacities
  explicit Compound_Analytic_MultigroupOpacity(const sf_char &packed);

  // Check that the class invariant is not violated.
  bool check_class_invariant() const;

  // >>> ACCESSORS
  const_Model get_Analytic_Model(size_t g) const { return group_models[g - 1]; }

  // >>> INTERFACE SPECIFIED BY rtt_cdi::MultigroupOpacity

  // Get the group opacities.
  sf_double getOpacity(double temperature, double density) const override;

  // Get the group opacity fields given a field of temperatures.
  vf_double getOpacity(const sf_double &temperature, double density) const override;

  // Get the group opacity fields given a field of densities.
  vf_double getOpacity(double temperature, const sf_double &density) const override;

  // Get the data description of the opacity.
  inline std_string getDataDescriptor() const override;

  // Pack the Compound_Analytic_MultigroupOpacity into a character string.
  sf_char pack() const override;
};

//------------------------------------------------------------------------------------------------//
// INLINE FUNCTIONS
//------------------------------------------------------------------------------------------------//
//! Return a string describing the opacity model.
Compound_Analytic_MultigroupOpacity::std_string
Compound_Analytic_MultigroupOpacity::getDataDescriptor() const {
  std_string descriptor;

  rtt_cdi::Reaction const rxn = getReactionType();

  if (rxn == rtt_cdi::TOTAL)
    descriptor = "Compound Multigroup Total";
  else if (rxn == rtt_cdi::ABSORPTION)
    descriptor = "Compound Multigroup Absorption";
  else if (rxn == rtt_cdi::SCATTERING)
    descriptor = "Compound Multigroup Scattering";
  else {
    Insist(rxn == rtt_cdi::TOTAL || rxn == rtt_cdi::ABSORPTION || rxn == rtt_cdi::SCATTERING,
           "Invalid Compound multigroup model opacity!");
  }

  return descriptor;
}

} // namespace rtt_cdi_analytic

#endif // rtt_cdi_analytic_Compound_Analytic_MultigroupOpacity_hh

//------------------------------------------------------------------------------------------------//
// end of cdi_analytic/Compound_Analytic_MultigroupOpacity.hh
//------------------------------------------------------------------------------------------------//
