//--------------------------------------------*-C++-*---------------------------------------------//
/*! \file   UnitSystem.hh
 *  \author Kelly Thompson
 *  \brief  Provide a definition of a unit system (7 dimensions: length,  mass, time, temperature,
 *          current, angle, quantity).
 *  \date   Fri Oct 24 15:07:43 2003
 *  \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_units_UnitSystem_hh
#define rtt_units_UnitSystem_hh

#include "UnitSystemType.hh"
#include "ds++/Assert.hh"
#include <cmath>

namespace rtt_units {

//================================================================================================//
/*! \class UnitSystem
 *
 * \brief Provide a units system object for Draco.
 *
 * This unit system consists of 7 unit types (Length, Mass, time, Temperature, Current, Angle and
 * Quantity).  The class then provides simple accessors that allow conversion between unit systems
 * and unit label names.  The general format for an accessor is u.L() or u.Lname().  L can be
 * replaced with M, t, T, I, A or Q to access the confersion factor for different fundamental units.
 * The accessor Xname() returns the label name for a unit.
 *
 * \verbatim
 * SI units are:
 *    length  (L) - meters
 *    mass    (M) - kilograms
 *    time    (t) - seconds
 *    temp    (T) - Kelvin
 *    current (I) - Amp
 *    angle (rad) - Radian
 *    quantity    - mole
 *
 * Derived units are:
 *    area             :         L^2
 *    volume           :         L^3
 *    velocity (v)     :         L   / t
 *    acceleration (a) :         L   / t^2
 *    density (rho)    :     M / L^3
 *    angular velocity : rad     / t
 *    force (f)        :     M * L   / t^2
 *    power (P)        :     M * L^2 / t^3
 *    torque           :     M * L^2 / t^2
 *    pressure (p)     :     M / L   / t^2
 *    energy           :     M * L^2 / t^2
 *
 *    charge                (q) : I * t     (Coulombs = Amp * sec)
 *    electrical potential  (V) : P / I     (Volt = Watt/Amp)
 *    electrical resistance (R) : V / I     (Ohm  = Volt/Amp)
 *    electrical capacitance (C): q / V     (Farads = Coulomb / Volt)
 *    electrical inductance (H) : V * t / I (Henry  = Volt * sec / Amp)
 *
 *    P = V * I
 *    V = I * R
 *    H = I / V / t
 *    C = I * t / V
 *
 *    luminous flux - lumen (lm) illuminance - lux (lx) \endverbatim
 *
 * \sa FundUnit
 * \sa PhysicalConstants
 *
 * \example test/tstUnitSystemType.cc
 * \example test/tstUnitSystem.cc
 *
 * Different ways to construct a UnitSystem
 *
 * \verbatim
 * using rtt_units::UnitSystemType;
 * using rtt_units::UnitSystem;
 * typedef rtt_units::UnitSystemType::
 *
 * UnitSystem uSI( UnitSystemType().SI() );
 * UnitSystem uX4( UnitSystemType().X4() );
 * UnitSystem uCGS( UnitSystemType().L( rtt_units::L_cm  )
 *                                  .M( rtt_units::M_g   )
 *                                  .t( rtt_units::t_s   )
 *                                  .T( rtt_units::T_keV )
 *                                  .I( rtt_units::I_amp )
 *                                  .A( rtt_units::A_rad )
 *                                  .Q( rtt_units::Q_mol ) );
 * UnitSystem uLM_only( UnitSystemType().L( rtt_units::L_cm )
 *                                      .M( rtt_units::M_g  );
 * \endverbatim
 */
//================================================================================================//

class UnitSystem {
public:
  // FRIENDS

  //! Define the equality operator for Units==Units.
  friend bool operator==(UnitSystem const &op1, UnitSystem const &op2);

  //! Define the equality operator for Units==Units.
  friend bool operator!=(UnitSystem const &op1, UnitSystem const &op2);

  // CREATORS

  //! Prefered constructor
  UnitSystem(UnitSystemType ust) // NOLINT [hicpp-explicit-conversions]
      : d_ust(std::move(ust)) {
    Require(validUnits());
  }

  //! Default constructor provides SI Units
  UnitSystem() : d_ust(UnitSystemType().SI()) { Require(validUnits()); }

  // ACCESSORS

  //! Return the conversion factor for length that when multiplied against values in user-units
  //! yields the value in SI-units.
  double L() const { return d_ust.L().cf(); }
  //! Return the label for this unit (e.g. cm, m, etc.).
  std::string Lname() const { return d_ust.L().label(); }

  //! Return the conversion factor for mass that when multiplied against values in user-units yields
  //! the value in SI-units.
  double M() const { return d_ust.M().cf(); }
  //! Return the label for this unit (e.g. kg, g, etc.).
  std::string Mname() const { return d_ust.M().label(); }

  //! Return the conversion factor for time that when multiplied against values in user-units yields
  //! the value in SI-units.
  double t() const { return d_ust.t().cf(); }
  //! Return the label for this unit (e.g. s, us, etc.).
  std::string tname() const { return d_ust.t().label(); }

  //! Return the conversion factor for temperature that when multiplied against values in user-units
  //! yields the value in SI-units.
  double T() const { return d_ust.T().cf(); }
  //! Return the label for this unit (e.g. K, keV, etc.).
  std::string Tname() const { return d_ust.T().label(); }

  //! Return the conversion factor for electric current that when multiplied against values in
  //! user-units yields the value in SI-units.
  double I() const { return d_ust.I().cf(); }
  //! Return the label for this unit (e.g. amp, etc.).
  std::string Iname() const { return d_ust.I().label(); }

  //! Return the conversion factor for angle that when multiplied against values in user-units
  //! yields the value in SI-units.
  double A() const { return d_ust.A().cf(); }
  //! Return the label for this unit (e.g. rad, deg, etc.).
  std::string Aname() const { return d_ust.A().label(); }

  //! Return the conversion factor for quantity that when multiplied against values in user-units
  //! yields the value in SI-units.
  double Q() const { return d_ust.Q().cf(); }
  //! Return the label for this unit (e.g. mol, etc.).
  std::string Qname() const { return d_ust.Q().label(); }

  //! Return the conversion factor for velocity
  double v() const { return this->L() / this->t(); }
  //! Return the conversion factor for acceleration
  double a() const { return this->L() / std::pow(this->t(), 2); }
  //! Return the conversion factor for force
  double f() const { return this->M() * this->L() / std::pow(this->t(), 2); }
  //! Return the conversion factor for energy
  double e() const { return this->M() * std::pow(this->L(), 2) / std::pow(this->t(), 2); }
  //! Return the conversion factor for power
  double p() const { return this->M() * std::pow(this->L(), 2) / std::pow(this->t(), 3); }

  // CLASS IMPLEMENTATION

  //! Check whether conversion units are in an acceptable range.
  bool validUnits() const;

private:
  // DATA

  UnitSystemType d_ust;

}; // end of class UnitSystem

} // end namespace rtt_units

#endif // rtt__units_UnitSystem_hh

//------------------------------------------------------------------------------------------------//
// end of UnitSystem.hh
//------------------------------------------------------------------------------------------------//
