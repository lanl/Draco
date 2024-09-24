//--------------------------------------------*-C++-*---------------------------------------------//
/*! \file   UnitSystem.cc
 *  \author Kelly Thompson
 *  \date   Thu Oct 24 15:10:32 2003
 *  \note   Copyright (C) 2003-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "UnitSystem.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "units/UnitSystemEnums.hh"
#include <limits>

namespace rtt_units {

//! Ensure that all unit conversions are valid (larger than some minimum value.)
bool UnitSystem::validUnits() const {
  double const minConversion(std::numeric_limits<double>::min());

  // Only perform check if the fundamental unit type has been defined.
  if (d_ust.L().enumVal() != L_null && d_ust.L().cf() < minConversion)
    return false;
  if (d_ust.M().enumVal() != M_null && d_ust.M().cf() < minConversion)
    return false;
  if (d_ust.t().enumVal() != t_null && d_ust.t().cf() < minConversion)
    return false;
  if (d_ust.T().enumVal() != T_null && d_ust.T().cf() < minConversion)
    return false;
  if (d_ust.I().enumVal() != I_null && d_ust.I().cf() < minConversion)
    return false;
  if (d_ust.A().enumVal() != A_null && d_ust.A().cf() < minConversion)
    return false;
  if (d_ust.Q().enumVal() != Q_null && d_ust.Q().cf() < minConversion)
    return false;

  // If we get here then the units are ok.
  return true;
} // end validUnits()

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return true if op1 and op2 are identical.
 *
 * \return true if conversion data members are the same between op1 and op2. Otherwise return false.
 *
 * For example:
 * \verbatim
 * UnitSystem UserUnits(34.0, 60.0, 0.0003, 99);
 * UnitSystem SIUnits(1.0, 1.0, 1.0 1.0 );
 *
 * Units NewUnits = UserUnits/SIUnits
 * Ensure( NewUnits == UserUnits );
 * \endverbatim
 */
bool operator==(UnitSystem const &op1, UnitSystem const &op2) {
  return rtt_dsxx::soft_equiv(op1.L(), op2.L()) && rtt_dsxx::soft_equiv(op1.M(), op2.M()) &&
         rtt_dsxx::soft_equiv(op1.t(), op2.t()) && rtt_dsxx::soft_equiv(op1.T(), op2.T()) &&
         rtt_dsxx::soft_equiv(op1.I(), op2.I()) && rtt_dsxx::soft_equiv(op1.A(), op2.A()) &&
         rtt_dsxx::soft_equiv(op1.Q(), op2.Q());
}

bool operator!=(UnitSystem const &op1, UnitSystem const &op2) { return !(op1 == op2); }

} // end namespace rtt_units

//------------------------------------------------------------------------------------------------//
// end of UnitSystem.cc
//------------------------------------------------------------------------------------------------//
