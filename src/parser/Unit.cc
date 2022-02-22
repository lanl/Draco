//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   Unit.cc
 * \author Kent G. Budge
 * \date   Wed Jan 22 15:18:23 MST 2003
 * \brief  Definitions of Unit methods.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Unit.hh"
#include <iostream>
#include <limits>
#include <sstream>

namespace {
using namespace std;

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Write text for a component of a unit.
 *
 * This is a helper function for operator<<(ostream, const Unit &).  It has
 * the logic for ensuring that dashes are inserted where necessary, e.g., the
 * dash in kg-m^2/s^2.
 *
 * \param str  Stream to which to write the unit component.
 * \param dash A dash is required after the last component written to the
 *             stream.  May be modified on return.
 * \param value Unit dimension value
 * \param name Unit name
 */
void dash_insert(ostream &str, bool &dash, double const value, char const *const name) {
  double const eps = std::numeric_limits<double>::epsilon();
  double const mrv = std::numeric_limits<double>::min();
  if (std::abs(value) > mrv) {
    if (dash) {
      str << '-';
    }
    if (!rtt_dsxx::soft_equiv(value, 1.0, eps)) {
      str << name << '^' << value;
    } else {
      str << name;
    }
    dash = true;
  }
}

} // end anonymous namespace

namespace rtt_parser {
using namespace std;

//------------------------------------------------------------------------------------------------//
/*!
 * \param str Stream to which to write the text description.
 * \param u Unit to write the text description for.
 * \return A reference to s.
 */
std::ostream &operator<<(std::ostream &str, const Unit &u) {
  str << u.conv << ' ';
  bool dash = false;

  dash_insert(str, dash, u.m, "m");
  dash_insert(str, dash, u.kg, "kg");
  dash_insert(str, dash, u.s, "s");
  dash_insert(str, dash, u.A, "A");
  dash_insert(str, dash, u.K, "K");
  dash_insert(str, dash, u.mol, "mol");
  dash_insert(str, dash, u.cd, "cd");
  dash_insert(str, dash, u.rad, "rad");
  dash_insert(str, dash, u.sr, "sr");
  return str;
}

//------------------------------------------------------------------------------------------------//
double conversion_factor(Unit const &units, Unit const &unit_system) {
  using std::pow;

  double const conv = pow(unit_system.m, units.m) * pow(unit_system.kg, units.kg) *
                      pow(unit_system.s, units.s) * pow(unit_system.A, units.A) *
                      pow(unit_system.K, units.K) * pow(unit_system.mol, units.mol) *
                      pow(unit_system.cd, units.cd) * pow(unit_system.rad, units.rad) *
                      pow(unit_system.sr, units.sr);

  return conv;
}

//------------------------------------------------------------------------------------------------//
double conversion_factor(Unit const &units, rtt_units::UnitSystem const &unit_system) {
  using std::pow;

  double const conv = pow(unit_system.L(), units.m) * pow(unit_system.M(), units.kg) *
                      pow(unit_system.t(), units.s) * pow(unit_system.T(), units.K) *
                      pow(unit_system.I(), units.A) * pow(unit_system.A(), units.rad) *
                      pow(unit_system.Q(), units.mol);

  return conv;
}

} // namespace rtt_parser

//------------------------------------------------------------------------------------------------//
// end of Unit.cc
//------------------------------------------------------------------------------------------------//
