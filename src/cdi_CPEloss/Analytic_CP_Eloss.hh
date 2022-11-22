//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_CPEloss/Analytic_CP_Eloss.hh
 * \author Kendra P. Long
 * \date   Fri Aug  2 14:28:08 2019
 * \brief  Analytic_CP_Eloss class definition.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_CPEloss_Analytic_CP_Eloss_hh
#define rtt_cdi_CPEloss_Analytic_CP_Eloss_hh

#include "Analytic_Eloss_Model.hh"
#include "cdi/CPCommon.hh"
#include "cdi/CPEloss.hh"
#include <memory>

namespace rtt_cdi_cpeloss {
//================================================================================================//
/*!
 * \class Analytic_CP_Eloss
 *
 * \brief Derived rtt_cdi::CPEloss class for analytic eloss.  This class implements the interface
 *        found in cdi/CPEloss.hh for the case where CP energy loss data is in analytic form.
 */
//================================================================================================//

class Analytic_CP_Eloss : public rtt_cdi::CPEloss {
public:
  // Useful typedefs.
  using const_SP_Model = std::shared_ptr<const Analytic_Eloss_Model>;
  using SP_Analytic_Model = std::shared_ptr<Analytic_Eloss_Model>;
  using sf_double = std::vector<double>;
  using std_string = std::string;
  using CParticle = rtt_cdi::CParticle;
  using CPModelAngleCutoff = rtt_cdi::CPModelAngleCutoff;

private:
  // Analytic eloss model.
  SP_Analytic_Model analytic_model;

  // Data model (e.g. Eloss, large-angle scatter, etc.)  Currently only ELOSS is implemented.
  // rtt_cdi::CPModel cpmodel;

public:
  // Constructor.
  Analytic_CP_Eloss(SP_Analytic_Model model_in, CParticle target_in, CParticle projectile_in,
                    CPModelAngleCutoff model_angle_cutoff_in);

  // >>> ACCESSORS
  const_SP_Model get_Analytic_Model() const { return analytic_model; }

  // >>> INTERFACE SPECIFIED BY rtt_cdi::CPEloss

  // Get an eloss.
  double getEloss(const double temperature, const double density, const double v0) const override;

  //! Query to see if data is in tabular or functional form (false).
  constexpr static bool is_data_in_tabular_form() { return false; }

  //! Get the temperature grid (size 0 for function-based analytic data).
  sf_double getTemperatureGrid() const override { return sf_double(0); }

  //! Get the density grid (size 0 for function-based analytic data).
  sf_double getDensityGrid() const override { return sf_double(0); }

  //! Get the density grid (size 0 for function-based analytic data).
  sf_double getEnergyGrid() const override { return sf_double(0); }

  //! Get the size of the temperature grid (size 0).
  size_t getNumTemperatures() const override { return 0; }

  //! Get the size of the density grid (size 0).
  size_t getNumDensities() const override { return 0; }

  //! Get the size of the density grid (size 0).
  size_t getNumEnergies() const override { return 0; }

  /*!
   * \brief Returns the general eloss model type
   *
   * Since this is an analytic model, return 1 (rtt_cdi::ANALYTIC_ETYPE)
   */
  rtt_cdi::CPModelType getModelType() const { return rtt_cdi::CPModelType::ANALYTIC_ETYPE; }
};

} // namespace rtt_cdi_cpeloss

#endif

//------------------------------------------------------------------------------------------------//
// End cdi_CPEloss/Analytic_CP_Eloss.hh
//------------------------------------------------------------------------------------------------//
