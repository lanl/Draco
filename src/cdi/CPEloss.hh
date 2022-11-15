//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi/CPEloss.hh
 * \author Kendra Long
 * \date   Fri Aug  2 14:28:08 2019
 * \brief  CPEloss class header file (an abstract class)
 * \note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_CPEloss_hh
#define rtt_cdi_CPEloss_hh

#include "CPCommon.hh"
#include "ds++/config.h"
#include <string>
#include <vector>

namespace rtt_cdi {

//================================================================================================//
/*!
 * \class CPEloss
 *
 * \brief This is a pure virtual class that defines a standard interface for all derived CPEloss
 *        objects.
 *
 * Any derived CPEloss object must provide as a minumum the functionality outlined in this routine.
 * This functionality includes access to the data grid and the ability to return calculated eloss
 * values.
 */
//================================================================================================//

class CPEloss {
public:
  CPEloss(CParticle target_in, CParticle projectile_in, CPModelType model_type_in,
          CPModelAngleCutoff model_angle_cutoff_in)
      : projectile(projectile_in), target(target_in), model_type(model_type_in),
        model_angle_cutoff(model_angle_cutoff_in) {}
  // ---------- //
  // Destructor //
  // ---------- //

  /*!
   * \brief Default CPEloss() destructor.
   *
   * This is required to correctly release memory when any object derived from CPEloss is destroyed.
   */
  virtual ~CPEloss() = default;
  //! Disable copy construction
  CPEloss(CPEloss const &rhs) = delete;
  //! Disable move construction
  CPEloss(CPEloss &&rhs) noexcept = delete;
  //! Disable assignment
  CPEloss &operator=(CPEloss const &rhs) = delete;
  //! Disable move-assignment
  CPEloss &operator=(CPEloss &&rhs) noexcept = delete;

  // --------- //
  // Accessors //
  // --------- //

  /*!
   * \brief Eloss accessor that returns a single eloss that corresponds to the provided temperature,
   *     density, and incident particle speed.
   *
   * \param targetTemperature The temperature value for which an eloss value is being requested
   *     (keV).
   *
   * \param targetDensity The density value for which an eloss value is being requested (g/cm^3).
   *
   * \param partSpeed The incident speed of the particle (cm/shk).
   *
   * \return A single Eloss rate coefficient (shk^-1).
   */
  virtual double getEloss(const double targetTemperature, const double targetDensity,
                          const double partSpeed) const = 0;
  /*!
   * \brief Query whether the data is in tables or functional form.
   */
  //virtual static constexpr bool is_data_in_tabular_form() = 0;

  /*!
   * \brief Returns a vector of temperatures that define the cached Eloss data table (keV).
   */
  virtual std::vector<double> getTemperatureGrid() const = 0;

  /*!
   * \brief Returns a vector of densities that define the cached Eloss data table (g/cm^3).
   */
  virtual std::vector<double> getDensityGrid() const = 0;

  /*!
   * \brief Returns a vector of energies that define the cached Eloss data table (keV).
   */
  virtual std::vector<double> getEnergyGrid() const = 0;

  /*!
   * \brief Returns the size of the temperature grid.
   */
  virtual size_t getNumTemperatures() const = 0;

  /*!
   * \brief Returns the size of the density grid.
   */
  virtual size_t getNumDensities() const = 0;

  /*!
   * \brief Returns the size of the energy grid.
   */
  virtual size_t getNumEnergies() const = 0;

  /*!
   * \brief Query to determine the transporting particle type.
   */
  CParticle getProjectile() const { return projectile; }

  /*!
   * \brief Query to determine the target species type.
   */
  CParticle getTarget() const { return target; }

  /*!
   * \brief Query to determine the CP Model type.
   */
  CPModelType getModelType() const { return model_type; }

  /*!
   * \brief Query to determine the angle cutoff used by the CP model.
   */
  CPModelAngleCutoff getModelAngleCutoff() const { return model_angle_cutoff; }

protected:
  // Particle being transported i.e. subject to energy loss
  const CParticle projectile;

  // Target particle
  const CParticle target;

  // CP Model type
  const CPModelType model_type;

  // CP Model angle cutoff
  const CPModelAngleCutoff model_angle_cutoff;
};

} // namespace rtt_cdi

#endif

//------------------------------------------------------------------------------------------------//
// End cdi/CPEloss.hh
//------------------------------------------------------------------------------------------------//
