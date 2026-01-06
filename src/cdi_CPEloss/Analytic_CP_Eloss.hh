//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_CPEloss/Analytic_CP_Eloss.hh
 * \author Kendra P. Long
 * \date   Fri Aug  2 14:28:08 2019
 * \brief  Analytic_CP_Eloss class definition.
 * \note   Copyright (C) 2020-2026 Triad National Security, LLC., All rights reserved. */
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
  using CPModelAngleCutoff = rtt_cdi::CPModelAngleCutoff;

private:
  // Analytic eloss model.
  SP_Analytic_Model analytic_model;

  // Data model (e.g. Eloss, large-angle scatter, etc.)  Currently only ELOSS is implemented.
  // rtt_cdi::CPModel cpmodel;

public:
  // Constructor.
  Analytic_CP_Eloss(SP_Analytic_Model model_in, CPModelAngleCutoff model_angle_cutoff_in);

  // >>> ACCESSORS
  [[nodiscard]] const_SP_Model get_Analytic_Model() const { return analytic_model; }

  // >>> INTERFACE SPECIFIED BY rtt_cdi::CPEloss

  // Get an eloss.
  [[nodiscard]] double getEloss(const double temperature, const double ndensity, const double v0,
                                const double debye_length = 0.0) const override;

  /*!
   * \brief Get target zaid
   */
  [[nodiscard]] int targ_zaid() const override { return analytic_model->targ_zaid(); }

  /*!
   * \brief Get target mass
   */
  [[nodiscard]] double targ_mass() const override { return analytic_model->targ_mass(); }

  /*!
   * \brief Get projectile zaid
   */
  [[nodiscard]] int proj_zaid() const override { return analytic_model->proj_zaid(); }

  /*!
   * \brief Get projectile mass
   */
  [[nodiscard]] double proj_mass() const override { return analytic_model->proj_mass(); }
  //! Query to see if data is in tabular or functional form (false).
  constexpr static bool is_data_in_tabular_form() { return false; }

  //! Get the temperature grid (size 0 for function-based analytic data).
  [[nodiscard]] sf_double getTemperatureGrid() const override { return sf_double(0); }

  //! Get the density grid (size 0 for function-based analytic data).
  [[nodiscard]] sf_double getDensityGrid() const override { return sf_double(0); }

  //! Get the density grid (size 0 for function-based analytic data).
  [[nodiscard]] sf_double getEnergyGrid() const override { return sf_double(0); }

  //! Get the size of the temperature grid (size 0).
  [[nodiscard]] size_t getNumTemperatures() const override { return 0; }

  //! Get the size of the density grid (size 0).
  [[nodiscard]] size_t getNumDensities() const override { return 0; }

  //! Get the size of the density grid (size 0).
  [[nodiscard]] size_t getNumEnergies() const override { return 0; }

  /*!
   * \brief Returns the general eloss model type
   *
   * Since this is an analytic model, return 1 (rtt_cdi::ANALYTIC_ETYPE)
   */
  [[nodiscard]] rtt_cdi::CPModelType getModelType() const {
    return rtt_cdi::CPModelType::ANALYTIC_ETYPE;
  }
};

} // namespace rtt_cdi_cpeloss

#endif

//------------------------------------------------------------------------------------------------//
// End cdi_CPEloss/Analytic_CP_Eloss.hh
//------------------------------------------------------------------------------------------------//
