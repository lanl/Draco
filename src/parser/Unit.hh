//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/Unit.hh
 * \author Kent Budge
 * \brief  Definition the Unit struct
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_parser_Unit_hh
#define rtt_parser_Unit_hh

#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstantsSI.hh"

namespace rtt_parser {

//================================================================================================//
/*!
 * \class Unit
 * \brief Define units and conversions to SI
 *
 * Unit is a POD struct describing a physical unit.  It gives the dimensions of the unit in terms of
 * the nine fundamental SI units as well as the conversion factor to SI. The dimensions are
 * specified as real numbers because there are some physical quantities whose unit description
 * requires noninteger powers of the basic units.
 *
 * Several examples follow the struct definition.
 */
//================================================================================================//

struct Unit {
  // ----------------------------------------
  // DATA
  // ----------------------------------------

  double m;   //!< Length dimension
  double kg;  //!< Mass dimension
  double s;   //!< Time dimension
  double A;   //!< Current dimension
  double K;   //!< Temperature dimension
  double mol; //!< Quantity dimension
  double cd;  //!< Luminous intensity dimension
  double rad; //!< Plane angle dimension
  double sr;  //!< Solid angle dimension

  double conv; //!< Conversion factor
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute product of two units
 *
 * \param a First factor
 * \param b Second factor
 * \return Product of the two factors
 */

inline Unit operator*(Unit const &a, Unit const &b) {
  Unit Result{};

  Result.m = a.m + b.m;
  Result.kg = a.kg + b.kg;
  Result.s = a.s + b.s;
  Result.A = a.A + b.A;
  Result.K = a.K + b.K;
  Result.mol = a.mol + b.mol;
  Result.cd = a.cd + b.cd;
  Result.rad = a.rad + b.rad;
  Result.sr = a.sr + b.sr;

  Result.conv = a.conv * b.conv;

  return Result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute product of unit and scalar
 *
 * \param a Scalar factor
 * \param b Unit factor
 * \return Product of the two factors
 */
inline Unit operator*(double const a, Unit const &b) {
  Unit Result = b;

  Result.conv = a * b.conv;

  return Result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute product of unit and scalar
 *
 * \param b Unit factor
 * \param a Scalar factor
 * \return Product of the two factors
 */
inline Unit operator*(Unit const &b, double const a) { return operator*(a, b); }

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute ratio of two units
 *
 * \param a Numerator
 * \param b Denominator
 * \return Ratio of the two factors
 */
inline Unit operator/(Unit const &a, Unit const &b) {
  Unit Result{};

  Result.m = a.m - b.m;
  Result.kg = a.kg - b.kg;
  Result.s = a.s - b.s;
  Result.A = a.A - b.A;
  Result.K = a.K - b.K;
  Result.mol = a.mol - b.mol;
  Result.cd = a.cd - b.cd;
  Result.rad = a.rad - b.rad;
  Result.sr = a.sr - b.sr;

  Result.conv = a.conv / b.conv;

  return Result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute reciprocal of a unit (times a scalar)
 *
 * \param a Scalar numerator
 * \param b Unit denominator
 * \return Reciprocal of the unit times the scalar
 */
inline Unit operator/(double const a, Unit const &b) {
  Unit Result{};

  Result.m = -b.m;
  Result.kg = -b.kg;
  Result.s = -b.s;
  Result.A = -b.A;
  Result.K = -b.K;
  Result.mol = -b.mol;
  Result.cd = -b.cd;
  Result.rad = -b.rad;
  Result.sr = -b.sr;

  Result.conv = a / b.conv;

  return Result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute unit divided by a scalar.
 *
 * \param b Unit numerator
 * \param a Scalar denominator
 * \return unit divided by the scalar
 */
inline Unit operator/(Unit const &b, double const a) { return b * (1 / a); }

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute unit raised to the power of a scalar.
 *
 * \param b Unit base
 * \param a Scalar power
 * \return unit raised to power of the scalar
 */
inline Unit pow(Unit const &b, double const a) {
  Unit result = b;
  result.m *= a;
  result.kg *= a;
  result.s *= a;
  result.A *= a;
  result.K *= a;
  result.mol *= a;
  result.cd *= a;
  result.rad *= a;
  result.sr *= a;
  result.conv = std::pow(b.conv, a);

  return result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Test two units for equality
 *
 * Compare two Units for strict equality, in both dimensions and conversion factor.  To test only
 * the dimensions, ignoring the conversion factor, use the \c is_compatible function.
 *
 * \param a First factor
 * \param b Second factor
 * \return \c true if the units are identical; \c false otherwise
 */
inline bool operator==(Unit const &a, Unit const &b) {
  using rtt_dsxx::soft_equiv;
  return soft_equiv(a.m, b.m) && soft_equiv(a.kg, b.kg) && soft_equiv(a.s, b.s) &&
         soft_equiv(a.A, b.A) && soft_equiv(a.K, b.K) && soft_equiv(a.mol, b.mol) &&
         soft_equiv(a.cd, b.cd) && soft_equiv(a.rad, b.rad) && soft_equiv(a.sr, b.sr) &&
         soft_equiv(a.conv, b.conv);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Test two units for inequality
 *
 * Compare two Units for inequality.
 *
 * \param a First factor
 * \param b Second factor
 * \return \c !(a==b)
 */
inline bool operator!=(Unit const &a, Unit const &b) { return !(a == b); }

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Test two units for compatibility
 *
 * Compare two Units for equality in dimensions only.  Ignore the conversion factor.  For exact test
 * of equality, use the \c operator== function.
 *
 * \param a First factor
 * \param b Second factor
 * \return \c true if the units are compatible; \c false otherwise
 */
inline bool is_compatible(Unit const &a, Unit const &b) {
  return rtt_dsxx::soft_equiv(a.m, b.m) && rtt_dsxx::soft_equiv(a.kg, b.kg) &&
         rtt_dsxx::soft_equiv(a.s, b.s) && rtt_dsxx::soft_equiv(a.A, b.A) &&
         rtt_dsxx::soft_equiv(a.K, b.K) && rtt_dsxx::soft_equiv(a.mol, b.mol) &&
         rtt_dsxx::soft_equiv(a.cd, b.cd) && rtt_dsxx::soft_equiv(a.rad, b.rad) &&
         rtt_dsxx::soft_equiv(a.sr, b.sr);
}

//------------------------------------------------------------------------------------------------//
//! Write out the unit in text form.

std::ostream &operator<<(std::ostream &str, const Unit &u);

// Some useful examples

// No units

//! dimensionless quantity (pure number)
Unit const dimensionless = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

// Fundamental SI units

Unit const m = {1, 0, 0, 0, 0, 0, 0, 0, 0, 1};   //!< meter
Unit const kg = {0, 1, 0, 0, 0, 0, 0, 0, 0, 1};  //!< kilogram
Unit const s = {0, 0, 1, 0, 0, 0, 0, 0, 0, 1};   //!< second
Unit const A = {0, 0, 0, 1, 0, 0, 0, 0, 0, 1};   //!< ampere
Unit const K = {0, 0, 0, 0, 1, 0, 0, 0, 0, 1};   //!< Kelvin
Unit const mol = {0, 0, 0, 0, 0, 1, 0, 0, 0, 1}; //!< mole
Unit const cd = {0, 0, 0, 0, 0, 0, 1, 0, 0, 1};  //!< candela
Unit const rad = {0, 0, 0, 0, 0, 0, 0, 1, 0, 1}; //!< radian
Unit const sr = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1};  //!< steradian

// Derived SI units

Unit const C = {0, 0, 1, 1, 0, 0, 0, 0, 0, 1};     //!< coulomb
Unit const Hz = {0, 0, -1, 0, 0, 0, 0, 0, 0, 1};   //!< hertz
Unit const N = {1, 1, -2, 0, 0, 0, 0, 0, 0, 1};    //!< newton
Unit const J = {2, 1, -2, 0, 0, 0, 0, 0, 0, 1};    //!< joule
Unit const Pa = {-1, 1, -2, 0, 0, 0, 0, 0, 0, 1};  //!< pascal
Unit const W = {2, 1, -3, 0, 0, 0, 0, 0, 0, 1};    //!< watt
Unit const V = {2, 1, -3, -1, 0, 0, 0, 0, 0, 1};   //!< volt
Unit const F = {-2, -1, 4, 2, 0, 0, 0, 0, 0, 1};   //!< farad
Unit const ohm = {2, 1, -3, -2, 0, 0, 0, 0, 0, 1}; //!< ohm
Unit const S = {-2, -1, 3, 2, 0, 0, 0, 0, 0, 1};   //!< siemens
Unit const Wb = {2, 1, -2, -1, 0, 0, 0, 0, 0, 1};  //!< weber
Unit const T = {0, 1, -2, -1, 0, 0, 0, 0, 0, 1};   //!< tesla
Unit const H = {2, 1, -2, -2, 0, 0, 0, 0, 0, 1};   //!< henry
Unit const lm = {0, 0, 0, 0, 0, 0, 1, 0, 1, 1};    //!< lumen
Unit const lx = {-2, 0, 0, 0, 0, 0, 1, 0, 1, 1};   //!< lux

// CGS units

Unit const cm = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0.01}; //!< centimeter
Unit const g = {0, 1, 0, 0, 0, 0, 0, 0, 0, 1e-3};  //!< gram

Unit const dyne = {1, 1, -2, 0, 0, 0, 0, 0, 0, 1e-5}; //!< dyne
Unit const erg = {2, 1, -2, 0, 0, 0, 0, 0, 0, 1e-7};  //!< erg

// English units

Unit const inch = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0.0254};        //!< inch
Unit const foot = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0.3048};        //!< foot
Unit const lbm = {0, 1, 0, 0, 0, 0, 0, 0, 0, 0.45359237};     //!< pound mass
Unit const pound = {1, 1, -2, 0, 0, 0, 0, 0, 0, 4.448221615}; //!< pound force

// Miscellaneous units

//! Electron volts
Unit const eV = {2, 1, -2, 0, 0, 0, 0, 0, 0, rtt_units::electronChargeSI};
//! Thousands of electron volts
Unit const keV = {2, 1, -2, 0, 0, 0, 0, 0, 0, 1e3 * rtt_units::electronChargeSI};

Unit const sh = {0, 0, 1, 0, 0, 0, 0, 0, 0, 1e-8};   //!< shake
Unit const jerk = {2, 1, -2, 0, 0, 0, 0, 0, 0, 1e9}; //!< jerk

// Numbers for which no conversion is requested
Unit const constant = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0}; //!< used for numbers with no units
Unit const raw = {0, 0, 0, 0, 0,
                  0, 0, 0, 0, 1.0}; //!< another name for numbers with no units, i.e.,

//------------------------------------------------------------------------------------------------//
/*! Systems of units
 *
 * Here we let the "dimension" be the unit conversion factor.
 *
 * double m;            //!< Length dimension
 * double kg;           //!< Mass dimension
 * double s;            //!< Time dimension
 * double A;            //!< Current dimension
 * double K;            //!< Temperature dimension
 * double mol;          //!< Quantity dimension
 * double cd;           //!< Luminous intensity dimension
 * double rad;          //!< Plane angle dimension
 * double sr;           //!< Solid angle dimension
 */

Unit const MKS = {1., 1., 1., 1., 1., 1., 1., 1., 1., 0.};
Unit const CGS = {0.01, 0.001, 1., 1., 1., 1., 1., 1., 1., 0.};
Unit const CGMU = {0.01, 0.001, 1e-6, 1., 1., 1., 1., 1., 1., 0.};
Unit const CGSH = {0.01, 0.001, 1e-8, 1., 1e3 * rtt_units::EV2K, 1., 1., 1., 1., 0.};

//------------------------------------------------------------------------------------------------//
/*! Calculate conversion factor to a system of units. Assumes the units are
* initially MKS.
*/

double conversion_factor(Unit const &units, Unit const &unit_system);

//------------------------------------------------------------------------------------------------//
/*! Calculate conversion factor to a system of units. Assumes the units are
 * initially MKS.
 */

double conversion_factor(Unit const &units, rtt_units::UnitSystem const &unit_system);

} // end namespace rtt_parser

#endif // rtt_parser_Unit_hh

//------------------------------------------------------------------------------------------------//
// end of parser/Unit.hh
//------------------------------------------------------------------------------------------------//
