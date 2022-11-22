//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_analytic/Analytic_Gray_Opacity.hh
 * \author Thomas M. Evans
 * \date   Fri Aug 24 13:13:46 2001
 * \brief  Analytic_Gray_Opacity class definition.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_analytic_Analytic_Gray_Opacity_hh
#define rtt_cdi_analytic_Analytic_Gray_Opacity_hh

#include "Analytic_Models.hh"
#include "cdi/GrayOpacity.hh"
#include <memory>

namespace rtt_cdi_analytic {

//================================================================================================//
/*!
 * \class Analytic_Gray_Opacity
 *
 * \brief Derived rtt_cdi::GrayOpacity class for analytic opacities.
 *
 * The Analytic_Gray_Opacity class is a derived rtt_cdi::GrayOpacity class.  It provides analytic
 * opacity data. The specific analytic opacity model is derived from the
 * rtt_cdi_analytic::Analytic_Opacity_Model base class.  Several pre-built derived classes are
 * provided in Analytic_Models.hh.
 *
 * Clients of this class can provide any analytic model class as long as it conforms to the
 * rtt_cdi_analytic::Analytic_Opacity_Model interface.  This interface consists of a single
 * function, Analytic_Opacity_Model::calculate_opacity().
 *
 * Note that opacities are returned in units of cm^2/g. Thus the resultant opacity must be
 * multiplied by density to get units of 1/cm.  See the documentation in
 * rtt_cdi_analytic::Analytic_Model for more info.
 *
 * The constructors take a rtt_cdi::Reaction argument to determine the reaction type.  The
 * enumeration rtt_cdi::Reaction can have the value TOTAL, ABSORPTION, or SCATTERING.
 *
 * The default rtt_cdi::Model for an Analytic_Gray_Opacity is rtt_cdi::ANALYTIC.  However, this can
 * be overridden in the constructor.
 *
 * This class conforms to the interface specified by rtt_cdi::GrayOpacity and can be used with
 * rtt_cdi::CDI to get analytic opacities.
 *
 * \example cdi_analytic/test/tstAnalytic_Gray_Opacity.cc
 *
 * Example usage of Analytic_Gray_Opacity, Analytic_Opacity_Model, and their incorporation into
 * rtt_cdi::CDI.
 */
//================================================================================================//

class Analytic_Gray_Opacity : public rtt_cdi::GrayOpacity {
public:
  // Useful typedefs.
  using SP_Analytic_Model = std::shared_ptr<Analytic_Opacity_Model>;
  using const_SP_Model = std::shared_ptr<const Analytic_Opacity_Model>;
  using sf_double = std::vector<double>;
  using std_string = std::string;
  using sf_char = std::vector<char>;

private:
  // Analytic opacity model.
  SP_Analytic_Model analytic_model;

  // Reaction model.
  rtt_cdi::Reaction reaction;

  // CDI model.
  rtt_cdi::Model model;

public:
  // Constructor.
  Analytic_Gray_Opacity(SP_Analytic_Model model_in, rtt_cdi::Reaction reaction_in,
                        rtt_cdi::Model cdi_model_in = rtt_cdi::ANALYTIC);

  // Constructor for packed Analytic_Gray_Opacities.
  explicit Analytic_Gray_Opacity(const sf_char &packed);

  // >>> ACCESSORS
  const_SP_Model get_Analytic_Model() const { return analytic_model; }

  // >>> INTERFACE SPECIFIED BY rtt_cdi::GrayOpacity

  // Get an opacity.
  double getOpacity(double temperature, double density) const override;

  // Get an opacity field given a field of temperatures.
  sf_double getOpacity(const sf_double &temperature, double density) const override;

  // Get an opacity field given a field of densities.
  sf_double getOpacity(double temperature, const sf_double &density) const override;

  //! Query to see if data is in tabular or functional form (false).
  bool data_in_tabular_form() const override { return false; }

  //! Query to get the reaction type.
  rtt_cdi::Reaction getReactionType() const override { return reaction; }

  //! Query for model type.
  rtt_cdi::Model getModelType() const override { return model; }

  // Return the energy policy (gray).
  inline std_string getEnergyPolicyDescriptor() const override;

  // Get the data description of the opacity.
  inline std_string getDataDescriptor() const override;

  // Get the name of the associated data file.
  inline std_string getDataFilename() const override;

  //! Get the temperature grid (size 0 for function-based analytic data).
  sf_double getTemperatureGrid() const override { return sf_double(0); }

  //! Get the density grid (size 0 for function-based analytic data).
  sf_double getDensityGrid() const override { return sf_double(0); }

  //! Get the size of the temperature grid (size 0).
  size_t getNumTemperatures() const override { return 0; }

  //! Get the size of the density grid (size 0).
  size_t getNumDensities() const override { return 0; }

  // Pack the Analytic_Gray_Opacity into a character string.
  sf_char pack() const override;

  /*!
   * \brief Returns the general opacity model type, defined in OpacityCommon.hh
   *
   * Since this is an analytic model, return 1 (rtt_cdi::ANALYTIC_TYPE)
   */
  rtt_cdi::OpacityModelType getOpacityModelType() const override { return rtt_cdi::ANALYTIC_TYPE; }
};

//------------------------------------------------------------------------------------------------//
// INLINE FUNCTIONS
//------------------------------------------------------------------------------------------------//
//! Return the energy policy descriptor (gray for Analytic_Gray_Opacity).
Analytic_Gray_Opacity::std_string Analytic_Gray_Opacity::getEnergyPolicyDescriptor() const {
  return "gray";
}

//------------------------------------------------------------------------------------------------//
//! Return a string describing the opacity model.
Analytic_Gray_Opacity::std_string Analytic_Gray_Opacity::getDataDescriptor() const {
  std_string descriptor;

  if (reaction == rtt_cdi::TOTAL)
    descriptor = "Analytic Gray Total";
  else if (reaction == rtt_cdi::ABSORPTION)
    descriptor = "Analytic Gray Absorption";
  else if (reaction == rtt_cdi::SCATTERING)
    descriptor = "Analytic Gray Scattering";
  else
    Insist(0, "Invalid analytic gray model opacity!");

  return descriptor;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return NULL string for the data filename.
 */
Analytic_Gray_Opacity::std_string Analytic_Gray_Opacity::getDataFilename() const { return ""; }

} // end namespace rtt_cdi_analytic

#endif // rtt_cdi_analytic_Analytic_Gray_Opacity_hh

//------------------------------------------------------------------------------------------------//
// end of cdi_analytic/Analytic_Gray_Opacity.hh
//------------------------------------------------------------------------------------------------//
