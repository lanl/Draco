//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   parser/Unit.hh
 * \author Kent Budge
 * \brief  Definition the Unit struct
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_parser_Unit_hh
#define rtt_parser_Unit_hh

#include "dsxx/Soft_Equivalence.hh"
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
 * the dimensions, ignoring the conversion factor, use the is_compatible function.
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

// Some useful examples --

// 1. No units

//! dimensionless quantity (pure number)
Unit constexpr dimensionless = {
    .m = 0, .kg = 0, .s = 0, .A = 0, .K = 0, .mol = 0, .cd = 0, .rad = 0, .sr = 0, .conv = 1};

// 2. Fundamental SI units

Unit constexpr m = {.m = 1,
                    .kg = 0,
                    .s = 0,
                    .A = 0,
                    .K = 0,
                    .mol = 0,
                    .cd = 0,
                    .rad = 0,
                    .sr = 0,
                    .conv = 1}; //!< meter
Unit constexpr kg = {.m = 0,
                     .kg = 1,
                     .s = 0,
                     .A = 0,
                     .K = 0,
                     .mol = 0,
                     .cd = 0,
                     .rad = 0,
                     .sr = 0,
                     .conv = 1}; //!< kilogram
Unit constexpr s = {.m = 0,
                    .kg = 0,
                    .s = 1,
                    .A = 0,
                    .K = 0,
                    .mol = 0,
                    .cd = 0,
                    .rad = 0,
                    .sr = 0,
                    .conv = 1}; //!< second
Unit constexpr A = {.m = 0,
                    .kg = 0,
                    .s = 0,
                    .A = 1,
                    .K = 0,
                    .mol = 0,
                    .cd = 0,
                    .rad = 0,
                    .sr = 0,
                    .conv = 1}; //!< ampere
Unit constexpr K = {.m = 0,
                    .kg = 0,
                    .s = 0,
                    .A = 0,
                    .K = 1,
                    .mol = 0,
                    .cd = 0,
                    .rad = 0,
                    .sr = 0,
                    .conv = 1}; //!< Kelvin
Unit constexpr mol = {.m = 0,
                      .kg = 0,
                      .s = 0,
                      .A = 0,
                      .K = 0,
                      .mol = 1,
                      .cd = 0,
                      .rad = 0,
                      .sr = 0,
                      .conv = 1}; //!< mole
Unit constexpr cd = {.m = 0,
                     .kg = 0,
                     .s = 0,
                     .A = 0,
                     .K = 0,
                     .mol = 0,
                     .cd = 1,
                     .rad = 0,
                     .sr = 0,
                     .conv = 1}; //!< candela
Unit constexpr rad = {.m = 0,
                      .kg = 0,
                      .s = 0,
                      .A = 0,
                      .K = 0,
                      .mol = 0,
                      .cd = 0,
                      .rad = 1,
                      .sr = 0,
                      .conv = 1}; //!< radian
Unit constexpr sr = {.m = 0,
                     .kg = 0,
                     .s = 0,
                     .A = 0,
                     .K = 0,
                     .mol = 0,
                     .cd = 0,
                     .rad = 0,
                     .sr = 1,
                     .conv = 1}; //!< steradian

// 3 - Derived SI units

Unit constexpr C = {.m = 0,
                    .kg = 0,
                    .s = 1,
                    .A = 1,
                    .K = 0,
                    .mol = 0,
                    .cd = 0,
                    .rad = 0,
                    .sr = 0,
                    .conv = 1}; //!< coulomb
Unit constexpr Hz = {.m = 0,
                     .kg = 0,
                     .s = -1,
                     .A = 0,
                     .K = 0,
                     .mol = 0,
                     .cd = 0,
                     .rad = 0,
                     .sr = 0,
                     .conv = 1}; //!< hertz
Unit constexpr N = {.m = 1,
                    .kg = 1,
                    .s = -2,
                    .A = 0,
                    .K = 0,
                    .mol = 0,
                    .cd = 0,
                    .rad = 0,
                    .sr = 0,
                    .conv = 1}; //!< newton
Unit constexpr J = {.m = 2,
                    .kg = 1,
                    .s = -2,
                    .A = 0,
                    .K = 0,
                    .mol = 0,
                    .cd = 0,
                    .rad = 0,
                    .sr = 0,
                    .conv = 1}; //!< joule
Unit constexpr Pa = {.m = -1,
                     .kg = 1,
                     .s = -2,
                     .A = 0,
                     .K = 0,
                     .mol = 0,
                     .cd = 0,
                     .rad = 0,
                     .sr = 0,
                     .conv = 1}; //!< pascal
Unit constexpr W = {.m = 2,
                    .kg = 1,
                    .s = -3,
                    .A = 0,
                    .K = 0,
                    .mol = 0,
                    .cd = 0,
                    .rad = 0,
                    .sr = 0,
                    .conv = 1}; //!< watt
Unit constexpr V = {.m = 2,
                    .kg = 1,
                    .s = -3,
                    .A = -1,
                    .K = 0,
                    .mol = 0,
                    .cd = 0,
                    .rad = 0,
                    .sr = 0,
                    .conv = 1}; //!< volt
Unit constexpr F = {.m = -2,
                    .kg = -1,
                    .s = 4,
                    .A = 2,
                    .K = 0,
                    .mol = 0,
                    .cd = 0,
                    .rad = 0,
                    .sr = 0,
                    .conv = 1}; //!< farad
Unit constexpr ohm = {.m = 2,
                      .kg = 1,
                      .s = -3,
                      .A = -2,
                      .K = 0,
                      .mol = 0,
                      .cd = 0,
                      .rad = 0,
                      .sr = 0,
                      .conv = 1}; //!< ohm
Unit constexpr S = {.m = -2,
                    .kg = -1,
                    .s = 3,
                    .A = 2,
                    .K = 0,
                    .mol = 0,
                    .cd = 0,
                    .rad = 0,
                    .sr = 0,
                    .conv = 1}; //!< siemens
Unit constexpr Wb = {.m = 2,
                     .kg = 1,
                     .s = -2,
                     .A = -1,
                     .K = 0,
                     .mol = 0,
                     .cd = 0,
                     .rad = 0,
                     .sr = 0,
                     .conv = 1}; //!< weber
Unit constexpr T = {.m = 0,
                    .kg = 1,
                    .s = -2,
                    .A = -1,
                    .K = 0,
                    .mol = 0,
                    .cd = 0,
                    .rad = 0,
                    .sr = 0,
                    .conv = 1}; //!< tesla
Unit constexpr H = {.m = 2,
                    .kg = 1,
                    .s = -2,
                    .A = -2,
                    .K = 0,
                    .mol = 0,
                    .cd = 0,
                    .rad = 0,
                    .sr = 0,
                    .conv = 1}; //!< henry
Unit constexpr lm = {.m = 0,
                     .kg = 0,
                     .s = 0,
                     .A = 0,
                     .K = 0,
                     .mol = 0,
                     .cd = 1,
                     .rad = 0,
                     .sr = 1,
                     .conv = 1}; //!< lumen
Unit constexpr lx = {.m = -2,
                     .kg = 0,
                     .s = 0,
                     .A = 0,
                     .K = 0,
                     .mol = 0,
                     .cd = 1,
                     .rad = 0,
                     .sr = 1,
                     .conv = 1}; //!< lux

// 4 - CGS units

Unit constexpr cm = {.m = 1,
                     .kg = 0,
                     .s = 0,
                     .A = 0,
                     .K = 0,
                     .mol = 0,
                     .cd = 0,
                     .rad = 0,
                     .sr = 0,
                     .conv = 0.01}; //!< centimeter
Unit constexpr g = {.m = 0,
                    .kg = 1,
                    .s = 0,
                    .A = 0,
                    .K = 0,
                    .mol = 0,
                    .cd = 0,
                    .rad = 0,
                    .sr = 0,
                    .conv = 1e-3}; //!< gram

Unit constexpr dyne = {.m = 1,
                       .kg = 1,
                       .s = -2,
                       .A = 0,
                       .K = 0,
                       .mol = 0,
                       .cd = 0,
                       .rad = 0,
                       .sr = 0,
                       .conv = 1e-5}; //!< dyne
Unit constexpr erg = {.m = 2,
                      .kg = 1,
                      .s = -2,
                      .A = 0,
                      .K = 0,
                      .mol = 0,
                      .cd = 0,
                      .rad = 0,
                      .sr = 0,
                      .conv = 1e-7}; //!< erg

// 5 - English units

Unit constexpr inch = {.m = 1,
                       .kg = 0,
                       .s = 0,
                       .A = 0,
                       .K = 0,
                       .mol = 0,
                       .cd = 0,
                       .rad = 0,
                       .sr = 0,
                       .conv = 0.0254}; //!< inch
Unit constexpr foot = {.m = 1,
                       .kg = 0,
                       .s = 0,
                       .A = 0,
                       .K = 0,
                       .mol = 0,
                       .cd = 0,
                       .rad = 0,
                       .sr = 0,
                       .conv = 0.3048}; //!< foot
Unit constexpr lbm = {.m = 0,
                      .kg = 1,
                      .s = 0,
                      .A = 0,
                      .K = 0,
                      .mol = 0,
                      .cd = 0,
                      .rad = 0,
                      .sr = 0,
                      .conv = 0.45359237}; //!< pound mass
Unit constexpr pound = {.m = 1,
                        .kg = 1,
                        .s = -2,
                        .A = 0,
                        .K = 0,
                        .mol = 0,
                        .cd = 0,
                        .rad = 0,
                        .sr = 0,
                        .conv = 4.448221615}; //!< pound force

// 6 - Miscellaneous units

//! Electron volts
Unit constexpr eV = {.m = 2,
                     .kg = 1,
                     .s = -2,
                     .A = 0,
                     .K = 0,
                     .mol = 0,
                     .cd = 0,
                     .rad = 0,
                     .sr = 0,
                     .conv = rtt_units::electronChargeSI};
//! Thousands of electron volts
Unit constexpr keV = {.m = 2,
                      .kg = 1,
                      .s = -2,
                      .A = 0,
                      .K = 0,
                      .mol = 0,
                      .cd = 0,
                      .rad = 0,
                      .sr = 0,
                      .conv = 1e3 * rtt_units::electronChargeSI};

Unit constexpr sh = {.m = 0,
                     .kg = 0,
                     .s = 1,
                     .A = 0,
                     .K = 0,
                     .mol = 0,
                     .cd = 0,
                     .rad = 0,
                     .sr = 0,
                     .conv = 1e-8}; //!< shake
Unit constexpr jerk = {.m = 2,
                       .kg = 1,
                       .s = -2,
                       .A = 0,
                       .K = 0,
                       .mol = 0,
                       .cd = 0,
                       .rad = 0,
                       .sr = 0,
                       .conv = 1e9}; //!< jerk

// Numbers for which no conversion is requested
Unit constexpr constant = {.m = 0,
                           .kg = 0,
                           .s = 0,
                           .A = 0,
                           .K = 0,
                           .mol = 0,
                           .cd = 0,
                           .rad = 0,
                           .sr = 0,
                           .conv = 1.0}; //!< used for numbers with no units
Unit constexpr raw = {.m = 0,
                      .kg = 0,
                      .s = 0,
                      .A = 0,
                      .K = 0,
                      .mol = 0,
                      .cd = 0,
                      .rad = 0,
                      .sr = 0,
                      .conv = 1.0}; //!< another name for numbers with no units, i.e.,

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

Unit constexpr MKS = {.m = 1.,
                      .kg = 1.,
                      .s = 1.,
                      .A = 1.,
                      .K = 1.,
                      .mol = 1.,
                      .cd = 1.,
                      .rad = 1.,
                      .sr = 1.,
                      .conv = 0.};
Unit constexpr CGS = {.m = 0.01,
                      .kg = 0.001,
                      .s = 1.,
                      .A = 1.,
                      .K = 1.,
                      .mol = 1.,
                      .cd = 1.,
                      .rad = 1.,
                      .sr = 1.,
                      .conv = 0.};
Unit constexpr CGMU = {.m = 0.01,
                       .kg = 0.001,
                       .s = 1e-6,
                       .A = 1.,
                       .K = 1.,
                       .mol = 1.,
                       .cd = 1.,
                       .rad = 1.,
                       .sr = 1.,
                       .conv = 0.};
Unit constexpr CGSH = {.m = 0.01,
                       .kg = 0.001,
                       .s = 1e-8,
                       .A = 1.,
                       .K = 1e3 * rtt_units::EV2K,
                       .mol = 1.,
                       .cd = 1.,
                       .rad = 1.,
                       .sr = 1.,
                       .conv = 0.};

//------------------------------------------------------------------------------------------------//
//! Calculate conversion factor to a system of units. Assumes the units are initially MKS.
double conversion_factor(Unit const &units, Unit const &unit_system);

//------------------------------------------------------------------------------------------------//
//! Calculate conversion factor to a system of units. Assumes the units are initially MKS.
double conversion_factor(Unit const &units, rtt_units::UnitSystem const &unit_system);

} // end namespace rtt_parser

#endif // rtt_parser_Unit_hh

//------------------------------------------------------------------------------------------------//
// end of parser/Unit.hh
//------------------------------------------------------------------------------------------------//
