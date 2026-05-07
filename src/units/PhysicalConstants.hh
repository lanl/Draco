//--------------------------------------------*-C++-*---------------------------------------------//
/*! \file   PhysicalConstants.hh
 *  \author Kelly Thompson
 *  \brief  Provide a single place where physical constants (pi, speed of light, etc) are defined.
 *  \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_units_PhysicalConstants_hh
#define rtt_units_PhysicalConstants_hh

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
 * PhysicalConstants pc( u );
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
  [[nodiscard]] constexpr double planck() const { return d_planck; }
  //! see planck()
  [[nodiscard]] constexpr double h() const { return planck(); }

  //! access the Gas constant (units of energy/mol/temp)
  [[nodiscard]] constexpr double gasConstant() const { return d_gasConstant; }
  //! see gasConstant()
  [[nodiscard]] constexpr double R() const { return gasConstant(); }

  //! accesses the Boltzmann constant (Energy/Temp)
  [[nodiscard]] constexpr double boltzmann() const { return d_boltzmann; }
  //! see boltzmann()
  [[nodiscard]] constexpr double k() const { return boltzmann(); }

  //! accesses the electron charge (Charge)
  [[nodiscard]] constexpr double electronCharge() const { return d_electronCharge; }
  //! see electronCharge()
  [[nodiscard]] constexpr double e() const { return electronCharge(); }

  //! accesses the speed of light (units of velocity)
  [[nodiscard]] constexpr double speedOfLight() const { return d_cLight; }
  //! see speedOfLight()
  [[nodiscard]] constexpr double c() const { return speedOfLight(); }

  //! accesses the StefanBoltzmann constant (Work/Area/Temp^4 )
  [[nodiscard]] constexpr double stefanBoltzmann() const { return d_stefanBoltzmann; }
  //! see stefanBoltzmann()
  [[nodiscard]] constexpr double sigma() const { return stefanBoltzmann(); }

  //! accesses the gravitational constant
  [[nodiscard]] constexpr double gravitationalConstant() const { return d_gravitationalConstant; }
  //! see gravitationalConstant()
  [[nodiscard]] constexpr double G() const { return gravitationalConstant(); }

  //! access the acceleration due to gravity (standard).
  [[nodiscard]] constexpr double accelerationFromGravity() const {
    return d_accelerationFromGravity;
  }
  //! see accelerationFromGravity()
  [[nodiscard]] constexpr double g() const { return accelerationFromGravity(); }

  //! access the Faraday constant
  [[nodiscard]] constexpr double faradayConstant() const { return d_faradayConstant; }
  //! see faradayConstant()
  [[nodiscard]] constexpr double F() const { return faradayConstant(); }

  //! access the Permeability of vacuum (free space)
  [[nodiscard]] constexpr double permeabilityOfVacuum() const { return d_permeabilityOfVacuum; }
  //! see permeabilityOfVacuum()
  [[nodiscard]] constexpr double mu0() const { return permeabilityOfVacuum(); }

  //! accesses the permittivity of free space (units of force/length)
  [[nodiscard]] constexpr double permittivityOfFreeSpace() const {
    return d_permittivityOfFreeSpace;
  }
  //! see permittivityOfFreeSpace()
  [[nodiscard]] constexpr double epsi0() const { return permittivityOfFreeSpace(); }

  //! accesses the classical electron radius (units of length)
  [[nodiscard]] constexpr double classicalElectronRadius() const {
    return d_classicalElectronRadius;
  }
  //! see classicalElectronRadius()
  [[nodiscard]] constexpr double re() const { return classicalElectronRadius(); }

  //! accesses the electron mass (units of mass)
  [[nodiscard]] constexpr double electronMass() const { return d_electronMass; }
  //! see electronMass()
  [[nodiscard]] constexpr double me() const { return electronMass(); }

  //! accesses the proton mass (units of mass)
  [[nodiscard]] constexpr double protonMass() const { return d_protonMass; }
  //! see protonMass()
  [[nodiscard]] constexpr double mp() const { return protonMass(); }

  [[nodiscard]] double unit_convert_boltzmann(UnitSystem const &u) const;

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
