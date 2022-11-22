//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   cdi_CPEloss/Analytic_Eloss_Model.hh
 * \author Kendra P. Long, Ben R. Ryan
 * \date   Feb 21 2020
 * \brief  Analytic_Eloss_Model class definition.
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_cdi_CPEloss_Analytic_Eloss_Model_hh
#define rtt_cdi_CPEloss_Analytic_Eloss_Model_hh

#include "cdi/CPCommon.hh"
#include "units/PhysicalConstexprs.hh"
#include <cmath>
#include <vector>

namespace rtt_cdi_cpeloss {

//================================================================================================//
/*!
 * \class Analytic_Eloss_Model
 * \brief Analytic_Eloss_Model base class.
 *
 * This is a base class that defines the interface given to Analytic_Eloss_Model constructors.  The
 * user can define any derived model class that will work with these analytic opacity generation
 * classes as long as it implements the functions required, namely:
 *
 * \arg double calculate_eloss(double T, double rho)
 *
 * This class is a pure virtual base class.
 *
 * The returned eloss coefficient is a power with units of keV shk^-1.
 */
//================================================================================================//

class Analytic_Eloss_Model {
public:
  //! Virtual destructor for proper inheritance destruction.
  virtual ~Analytic_Eloss_Model() = default;

  // Disable copy/move ctors and assignment operators.
  Analytic_Eloss_Model(Analytic_Eloss_Model const &rhs) = delete;
  Analytic_Eloss_Model(Analytic_Eloss_Model &&rhs) noexcept = delete;
  Analytic_Eloss_Model &operator=(Analytic_Eloss_Model const &rhs) = delete;
  Analytic_Eloss_Model &operator=(Analytic_Eloss_Model &&rhs) noexcept = delete;

  //! Interface for derived analytic eloss models.
  virtual double calculate_eloss(const double T, const double rho, const double v0) const = 0;

protected:
  //! Constructor initializes generically useful physical quantities.
  Analytic_Eloss_Model(const rtt_cdi::CParticle &target, const rtt_cdi::CParticle &projectile)
      : zaidt(target.get_zaid()), mt(target.get_mass()), qtabs(std::abs(target.get_z()) * pc.e()),
        zaidp(target.get_zaid()), mp(projectile.get_mass()),
        qpabs(std::abs(projectile.get_z()) * pc.e()) {}

protected:
  //! Unit system (use cgs internally)
  rtt_units::PhysicalConstexprs<rtt_units::CGS> pc;

  const int zaidt;    //!< Target ZAID
  const double mt;    //!< Target mass
  const double qtabs; //!< Absolute target electric charge
  const int zaidp;    //!< Projectile ZAID
  const double mp;    //!< Projectile mass
  const double qpabs; //!< Absolute projectile electric charge
};

} // namespace rtt_cdi_cpeloss

#endif // rtt_cdi_CPEloss_Analytic_Eloss_Model_hh

//------------------------------------------------------------------------------------------------//
// End cdi_CPEloss/Analytic_Eloss_Model.hh
//------------------------------------------------------------------------------------------------//
