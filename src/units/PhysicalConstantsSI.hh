//--------------------------------------------*-C++-*---------------------------------------------//
/*! \file   PhysicalConstantsSI.hh
 *  \author Kelly Thompson, Kent G. Budge
 *  \brief  Provide a single place where physical constants (pi, speed of light, etc) are defined in
 *          SI units.
 *  \date   Fri Nov 07 10:04:52 2003
 *  \note   Copyright (C) 2014-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_units_PhysicalConstantsSI_hh
#define rtt_units_PhysicalConstantsSI_hh

#include "MathConstants.hh"
#include "UnitSystem.hh"

//! \namespace rtt_units Namespace for units and physical constants
namespace rtt_units {

// Base physical constants in SI units:

//    m - meters, kg - kilograms, s - seconds, K - kelvin
//    W - watts, J - joules, C - coulombs, F - farads
//    mol - mole

//------------------------------------------------------------------------------------------------//
/*!
 * FUNDAMENTAL CONSTANTS
 *
 * NIST Reference on Constants, Units and Uncertainty CODATA internationally recommended values of
 * the Fundamental Physical Constants, http://physics.nist.gov/cuu/Constants/
 *
 * The units of these fundamental constants should be factors of 10X different from the official
 * NIST 2010 CODATA report data to allow for easy comparison between these values and the NIST data.
 *
 * - Fundamental constants are listed first.
 * - Derived constants are listed second.
 * - Actual data is placed in a user-defined type for C-interoperability.
 */
//------------------------------------------------------------------------------------------------//

//! [c] SPEED OF LIGHT (M/S)
// exact value by NIST definition
static double constexpr cLightSI = 2.99792458e8; // m s^-1

//! [Na] AVOGADRO'S NUMBER ("entities"/MOLE)
// Wikipedia (2013-12-3) == NIST Codata 2010
static double constexpr AVOGADRO = 6.02214129e23; // entities/mol

//! [h] Planck constant ( Energy seconds )
// Wikipedia (2013-12-3) == NIST Codata 2010 (eps = 4.4e-8)
static double constexpr planckSI = 6.62606957e-34; // J s

//! [R] Molar Gas constant
// Wikipedia (2013-12-3) == NIST Codata 2010 (eps 9.1e-7)
static double constexpr gasConstantSI = 8.3144621; // J/mol/K

/*!
 * \brief [k] BOLTZMANN'S CONSTANT == R/Na (JOULES/K)
 *
 * \note If this changes you msut update the Enumerated Temperature Type in UnitSystemEnusm.hh!
 */
static double constexpr boltzmannSI = 1.380648800E-23; // J K^-1

/*!
 * \brief [e] ELECTRON CHARGE (COULOMBS)
 *
 * Wikipedia (2013-12-3) == NIST Codata 2010 (eps = 2.2e-8)
 * \note If this changes you must update the Enumerated Temperature Type in UnitSystemEnusm.hh!
 */
static double constexpr electronChargeSI = 1.602176565e-19; // Amp / sec

//! [me] ELECTRON REST MASS (KG)	 s
// Wikipedia (2013-12-3) == NIST Codata 2010 (eps = 4e-8)
static double constexpr electronMassSI = 9.10938291e-31; // kg

//! [G] Gravitational Constant
// Wikipedia (2013-12-3) == NIST Codata 2010 (eps = 1.2e-4)
static double constexpr gravitationalConstantSI = 6.67384e-11; // m^3/kg/s^2

//! [g] Acceleration due to gravity (standard)
// Wikipedia (2013-12-3)
static double constexpr accelerationFromGravitySI = 9.80665; // m/s^2

//! [mp] PROTON REST MASS (KG)
// Wikipedia (2013-12-3) == NIST Codata 2010 (eps = 4.4e-8)
static double constexpr protonMassSI = 1.672621777e-27; // kg

//------------------------------------------------------------------------------------------------//
// DERIVED CONSTANTS
//  - constants derived from fundamental constants listed above
//  - constants with alternate units
//  - conversion factors
//------------------------------------------------------------------------------------------------//

/*! \brief EV2K CONVERSION FACTOR FROM ELECTRON-VOLTS TO KELVIN (K/eV)
 *
 * (TEMPERATURE IN eV) * EV2K = (TEMPERATURE IN KELVIN)
 *
 * \note If this number is changed, you must also update the conversion factor found in
 *       UniSystemUnums.hh.
 */
static double constexpr EV2K = electronChargeSI / boltzmannSI;

/*!
 * \brief [sigma] STEFAN-BOLTZMANN CONSTANT (WATTS/(M**2-K**4)
 *
 * \f[
 * \sigma_{SB} = \frac{2 \pi^5 k^4} {15 h^3 c^2}
 *             = 5.670373e-8
 * \f]
 */
static double constexpr stefanBoltzmannSI =
    2.0 * PI * PI * PI * PI * PI * boltzmannSI * boltzmannSI * boltzmannSI * boltzmannSI /
    (15.0 * planckSI * planckSI * planckSI * cLightSI * cLightSI);

//! [F] Faraday constant == Na * e
static double constexpr faradayConstantSI = AVOGADRO * electronChargeSI;

//! [mu0] Permeability of vacuum (free space)
static double constexpr permeabilityOfVacuumSI = 4.0 * PI * 1.0e-7; // Henry/m

//! [epsi0] PERMITTIVITY OF FREE SPACE (F/M)
static double constexpr permittivityOfFreeSpaceSI =
    1.0 / permeabilityOfVacuumSI / cLightSI / cLightSI; // Coulomb^2/J/m

//! [re] Classical electron radius (M)
static double constexpr classicalElectronRadiusSI =
    electronChargeSI * electronChargeSI /
    (4 * PI * permittivityOfFreeSpaceSI * electronMassSI * cLightSI * cLightSI);

} // end namespace rtt_units

#endif // rtt_units_PhysicalConstantsSI_hh

//------------------------------------------------------------------------------------------------//
// end of units/PhysicalConstantsSI.hh
//------------------------------------------------------------------------------------------------//
