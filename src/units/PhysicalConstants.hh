//--------------------------------------------*-C++-*---------------------------------------------//
/*! \file   PhysicalConstants.hh
 *  \author Kelly Thompson
 *  \brief  Provide a single place where physical constants (pi, speed of light, etc) are defined.
 *  \date   Fri Nov 07 10:04:52 2003
 *  \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_units_PhysicalConstants_hh
#define rtt_units_PhysicalConstants_hh

#include "MathConstants.hh"
#include "PhysicalConstantsSI.hh"
#include "UnitSystem.hh"

//! \namespace rtt_units Namespace for units and physical constants
namespace rtt_units {
//==============================================================================
/*!
 * \class PhysicalConstants
 *
 * \brief A class to define and encapsulate physical and mathematical constants in the current
 *        UnitSystem.
 *
 * \sa rtt_units::UnitSystem
 * \sa rtt_units::UnitSystemType
 * \sa UnitSystemEnums.hh
 * \sa rtt_units::FundUnit
 * \sa Nuclide Chart
 *
 * Code Sample:
 *
 * \verbatim
 * UnitSystem u( UnitSystem().getAstroUnits() );
 * PhyscialConstants pc( u );
 * double mass   = 1.0; // grams (AstroPhysics mass units).
 * double energy = mass * pc.c() * pc.c();  // E = mc^2 where energy is in
 *                                          // AstroPhysics units.
 * \endverbatim
 *
 * \example test/tstPhysicalConstants.cc
 * This is the unit regression test for the PhysicalConstants class.  It demonstrates typical usage.
 */
//================================================================================================//

class PhysicalConstants {
public:
  // Constructors.
  //----------------------------------------------------------------------------------------------//
  /*!
   * \brief Default constructor provides physical constants with SI units (kg, m, seconds, degree K,
   *        amp, radian, mole).
   */
  constexpr PhysicalConstants()
      : d_avogadro(AVOGADRO), d_planck(planckSI), d_gasConstant(gasConstantSI),
        d_boltzmann(boltzmannSI), d_electronCharge(electronChargeSI), d_cLight(cLightSI),
        d_stefanBoltzmann(stefanBoltzmannSI), d_gravitationalConstant(gravitationalConstantSI),
        d_accelerationFromGravity(accelerationFromGravitySI), d_faradayConstant(faradayConstantSI),
        d_permeabilityOfVacuum(permeabilityOfVacuumSI),
        d_permittivityOfFreeSpace(permittivityOfFreeSpaceSI),
        d_classicalElectronRadius(classicalElectronRadiusSI), d_electronMass(electronMassSI),
        d_protonMass(protonMassSI) {
    // empty
  }
  explicit PhysicalConstants(UnitSystem const &u);

  //! \todo Make electronCharge and Avaragodo adjustable based on units.

  //! accesses Avogadro's number (1/mole)
  constexpr double avogadro() { return d_avogadro; }

  //! see avogadro()
  constexpr double Na() { return avogadro(); }

  //! access the Planck constant (units of energy-time)
  constexpr double planck() const { return d_planck; }
  //! see planck()
  constexpr double h() const { return planck(); }

  //! access the Gas constant (units of energy/mol/temp)
  constexpr double gasConstant() const { return d_gasConstant; }
  //! see gasConstant()
  constexpr double R() const { return gasConstant(); }

  //! accesses the Boltzmann constant (Energy/Temp)
  constexpr double boltzmann() const { return d_boltzmann; }
  //! see boltzmann()
  constexpr double k() const { return boltzmann(); }

  //! accesses the electron charge (Charge)
  constexpr double electronCharge() const { return d_electronCharge; }
  //! see electronCharge()
  constexpr double e() const { return electronCharge(); }

  //! accesses the speed of light (units of velocity)
  constexpr double speedOfLight() const { return d_cLight; }
  //! see speedOfLight()
  constexpr double c() const { return speedOfLight(); }

  //! accesses the StefanBoltzmann constant (Work/Area/Temp^4 )
  constexpr double stefanBoltzmann() const { return d_stefanBoltzmann; }
  //! see stefanBoltzmann()
  constexpr double sigma() const { return stefanBoltzmann(); }

  //! accesses the gravitational constant
  constexpr double gravitationalConstant() const { return d_gravitationalConstant; }
  //! see gravitationalConstant()
  constexpr double G() const { return gravitationalConstant(); }

  //! access the acceleration due to gravity (standard).
  constexpr double accelerationFromGravity() const { return d_accelerationFromGravity; }
  //! see accelerationFromGravity()
  constexpr double g() const { return accelerationFromGravity(); }

  //! access the Faraday constant
  constexpr double faradayConstant() const { return d_faradayConstant; }
  //! see faradayConstant()
  constexpr double F() const { return faradayConstant(); }

  //! access the Permeability of vacuum (free space)
  constexpr double permeabilityOfVacuum() const { return d_permeabilityOfVacuum; }
  //! see permeabilityOfVacuum()
  constexpr double mu0() const { return permeabilityOfVacuum(); }

  //! accesses the permittivity of free space (units of force/length)
  constexpr double permittivityOfFreeSpace() const { return d_permittivityOfFreeSpace; }
  //! see permittivityOfFreeSpace()
  constexpr double epsi0() const { return permittivityOfFreeSpace(); }

  //! accesses the classical electron radius (units of length)
  constexpr double classicalElectronRadius() const { return d_classicalElectronRadius; }
  //! see classicalElectronRadius()
  constexpr double re() const { return classicalElectronRadius(); }

  //! accesses the electron mass (units of mass)
  constexpr double electronMass() const { return d_electronMass; }
  //! see electronMass()
  constexpr double me() const { return electronMass(); }

  //! accesses the proton mass (units of mass)
  constexpr double protonMass() const { return d_protonMass; }
  //! see protonMass()
  constexpr double mp() const { return protonMass(); }

  double unit_convert_boltzmann(UnitSystem const &u) const;

private:
  // Base physical constants in SI units:

  //    m - meters, kg - kilograms, s - seconds, K - kelvin
  //    W - watts, J - joules, C - coulombs, F - farads
  //    mol - mole

  //! [mol] Avogadro's number
  double const d_avogadro;

  //! [h] Planck constant ( Energy seconds )
  double const d_planck;

  //! [R] Gas constant
  double const d_gasConstant;

  //! [k] BOLTZMANN'S CONSTANT == R/Na (JOULES/K)
  double const d_boltzmann;

  //! [e] ELECTRON CHARGE (COULOMBS)
  double const d_electronCharge;

  //! [c] SPEED OF LIGHT (M/S)
  double const d_cLight;

  //! [sigma] STEFAN-BOLTZMANN CONSTANT (WATTS/(M**2-K**4)
  double const d_stefanBoltzmann;

  //! [G] Gravitational Constant
  double const d_gravitationalConstant; // SI units are N / m^2 / kg^2

  //! [g] Acceleration due to gravity (standard)
  double const d_accelerationFromGravity;

  //! [F] Faraday constant == Na * e
  double const d_faradayConstant; // SI units are coloumbs/mol

  //! [mu0] Permeability of vacuum (free space)
  double const d_permeabilityOfVacuum;

  //! [epsi0] PERMITTIVITY OF FREE SPACE (F/M)
  double const d_permittivityOfFreeSpace;

  //! [re] Classical electron radius (M)
  double const d_classicalElectronRadius;

  //! [me] ELECTRON REST MASS (KG)
  double const d_electronMass;

  //! [mp] PROTON REST MASS (KG)
  double const d_protonMass;

}; // end class PhysicalConstants

} // end namespace rtt_units

#endif // rtt_units_PhysicalConstants_hh

//------------------------------------------------------------------------------------------------//
// end of units/PhysicalConstants.hh
//------------------------------------------------------------------------------------------------//
