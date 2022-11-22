//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_analytic/Analytic_MultigroupOpacity.hh
 * \author Thomas M. Evans
 * \date   Tue Nov 13 11:19:59 2001
 * \brief  Analytic_MultigroupOpacity class definition.
 * \note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_analytic_Analytic_MultigroupOpacity_hh
#define rtt_cdi_analytic_Analytic_MultigroupOpacity_hh

#include "Analytic_Models.hh"
#include "cdi/MultigroupOpacity.hh"
#include <memory>

namespace rtt_cdi_analytic {

//================================================================================================//
/*!
 * \class Analytic_MultigroupOpacity
 *
 * \brief Derived rtt_cdi::MultigroupOpacity class for analytic opacities.
 *
 * The Analytic_MultigroupOpacity class is an abstract class derived from rtt_cdi::MultigroupOpacity
 * class.  It provides analytic opacity data. The specific analytic opacity model is specified by
 * concrete classes derived from Analytic_MultigroupOpacity.
 *
 * Note that opacities are returned in units of cm^2/g. Thus the resultant opacity must be
 * multiplied by density to get units of 1/cm.  See the documentation in
 * rtt_cdi_analytic::Analytic_Model for more info.
 *
 * The constructors take a rtt_cdi::Reaction argument to determine the reaction type.  The
 * enumeration rtt_cdi::Reaction can have the value TOTAL, ABSORPTION, or SCATTERING.
 *
 * The default rtt_cdi::Model for an Analytic_MultigroupOpacity is rtt_cdi::ANALYTIC.  However, this
 * can be overridden in the constructor.
 *
 * \sa cdi_analytic/nGray_Analytic_MultigroupOpacity.hh
 * Example usage of Analytic_MultigroupOpacity.
 */
//================================================================================================//

class Analytic_MultigroupOpacity : public rtt_cdi::MultigroupOpacity {
public:
  // Useful typedefs.
  using sf_double = std::vector<double>;
  using vf_double = std::vector<sf_double>;
  using sf_char = std::vector<char>;

private:
  // Group structure.
  sf_double group_boundaries;

  // Reaction model.
  rtt_cdi::Reaction reaction;

  // CDI model.
  rtt_cdi::Model model;

protected:
  // Constructor.
  Analytic_MultigroupOpacity(sf_double groups, rtt_cdi::Reaction reaction_in,
                             rtt_cdi::Model model_in = rtt_cdi::ANALYTIC);

  // Constructor for packed Analytic_Multigroup_Opacities
  explicit Analytic_MultigroupOpacity(const sf_char &packed);

  // Get the packed size of the object
  unsigned packed_size() const;

public:
  // >>> ACCESSORS

  // >>> INTERFACE SPECIFIED BY rtt_cdi::MultigroupOpacity

  ~Analytic_MultigroupOpacity() override = default;
  Analytic_MultigroupOpacity(Analytic_MultigroupOpacity const &rhs) = delete;
  Analytic_MultigroupOpacity(Analytic_MultigroupOpacity &&rhs) noexcept = delete;
  Analytic_MultigroupOpacity &operator=(Analytic_MultigroupOpacity const &rhs) = delete;
  Analytic_MultigroupOpacity &operator=(Analytic_MultigroupOpacity &&rhs) noexcept = delete;

  // Get the group opacities.
  sf_double getOpacity(double targetTemperature, double targetDensity) const override = 0;

  // Get the group opacity fields given a field of temperatures.
  vf_double getOpacity(const sf_double &targetTemperature, double targetDensity) const override = 0;

  // Get the group opacity fields given a field of densities.
  vf_double getOpacity(double targetTemperature, const sf_double &targetDensity) const override = 0;

  //! Query to see if data is in tabular or functional form (false).
  bool data_in_tabular_form() const override { return false; }

  //! Query to get the reaction type.
  rtt_cdi::Reaction getReactionType() const override { return reaction; }

  //! Query for model type.
  rtt_cdi::Model getModelType() const override { return model; }

  //! Return the energy policy descriptor (mg).
  inline std::string getEnergyPolicyDescriptor() const override { return "mg"; }

  // Get the data description of the opacity.
  std::string getDataDescriptor() const override = 0;

  // Get the name of the associated data file.
  inline std::string getDataFilename() const override { return ""; }

  //! Get the temperature grid (size 0 for function-based analytic data).
  sf_double getTemperatureGrid() const override { return sf_double(0); }

  //! Get the density grid (size 0 for function-based analytic data).
  sf_double getDensityGrid() const override { return sf_double(0); }

  //! Get the group boundaries (keV) of the multigroup set.
  sf_double getGroupBoundaries() const override { return group_boundaries; }

  //! Get the size of the temperature grid (size 0).
  size_t getNumTemperatures() const override { return 0; }

  //! Get the size of the density grid (size 0).
  size_t getNumDensities() const override { return 0; }

  //! Get the number of frequency group boundaries.
  size_t getNumGroupBoundaries() const override { return group_boundaries.size(); }

  //! Get the number of frequency group boundaries.
  size_t getNumGroups() const override { return group_boundaries.size() - 1; }

  // Pack the Analytic_MultigroupOpacity into a character string.
  sf_char pack() const override = 0;

  /*!
   * \brief Returns the general opacity model type, defined in OpacityCommon.hh
   *
   * Since this is an analytic model, return 1 (rtt_cdi::ANALYTIC_TYPE)
   */
  rtt_cdi::OpacityModelType getOpacityModelType() const override { return rtt_cdi::ANALYTIC_TYPE; }
};

} // end namespace rtt_cdi_analytic

#endif // rtt_cdi_analytic_Analytic_MultigroupOpacity_hh

//------------------------------------------------------------------------------------------------//
// end of cdi_analytic/Analytic_MultigroupOpacity.hh
//------------------------------------------------------------------------------------------------//
