//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_analytic/Pseudo_Line_Analytic_MultigroupOpacity.hh
 * \author Kent G. Budge
 * \date   Tue Apr  5 08:36:13 MDT 2011
 * \note   Copyright (C) 2011-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_analytic_Pseudo_Line_Analytic_MultigroupOpacity_hh
#define rtt_cdi_analytic_Pseudo_Line_Analytic_MultigroupOpacity_hh

#include "Analytic_MultigroupOpacity.hh"
#include "Pseudo_Line_Base.hh"

namespace rtt_cdi_analytic {
using rtt_parser::Expression;

//------------------------------------------------------------------------------------------------//
/*!
 * \class Pseudo_Line_Analytic_Opacity_Model
 *
 * \brief Derived Analytic_Opacity_Model class that defines a random line spectrum for the opacity.
 *
 * The opacity function is a continuum on which is superimposed a number of lines of the specified
 * peak and width. The line locations are chosen at random.
 *
 * The mass opacity coefficient is assumed independent of temperature or density, which allows
 * precalculation of the opacity structure, an important time saver.
 */
class Pseudo_Line_Analytic_MultigroupOpacity : public Analytic_MultigroupOpacity,
                                               public Pseudo_Line_Base {
private:
  Averaging averaging_;
  //! value of 0 indicates to use adaptive Romberg integration
  unsigned qpoints_;

  friend class PLR_Functor; //!< used in calculation of Rosseland averages
  friend class PLP_Functor; //!< used in calculation of Planck averages

public:
  Pseudo_Line_Analytic_MultigroupOpacity(sf_double const &group_bounds, rtt_cdi::Reaction reaction,
                                         std::shared_ptr<Expression const> const &cont,
                                         unsigned number_of_lines, double line_peak,
                                         double line_width, unsigned number_of_edges,
                                         double edge_ratio, double Tref, double Tpow, double emin,
                                         double emax, Averaging averaging, unsigned qpoints,
                                         unsigned seed);

  //! Constructor for packed state.
  explicit Pseudo_Line_Analytic_MultigroupOpacity(const sf_char &packed);

  //! Get the group opacities.
  sf_double getOpacity(double targetTemperature, double targetDensity) const override;

  //! Get the group opacity fields given a field of temperatures.
  vf_double getOpacity(const sf_double &targetTemperature, double targetDensity) const override;

  //! Get the group opacity fields given a field of densities.
  vf_double getOpacity(double targetTemperature, const sf_double &targetDensity) const override;

  //! Get the data description of the opacity.
  std::string getDataDescriptor() const override;

  //! Pack up the class for persistence.
  sf_char pack() const override;
};

} // end namespace rtt_cdi_analytic

#endif // rtt_cdi_analytic_Pseudo_Line_Analytic_MultigroupOpacity_hh

//------------------------------------------------------------------------------------------------//
// end of cdi_analytic/Pseudo_Line_Analytic_MultigroupOpacity.hh
//------------------------------------------------------------------------------------------------//
