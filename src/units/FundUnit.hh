//--------------------------------------------*-C++-*---------------------------------------------//
/*! \file   FundUnit.hh
 *  \author Kelly Thompson
 *  \brief  This file defines a fundamental unit type.
 *  \date   Mon Oct 27 16:24:31 2003
 *  \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_units_FundUnit_hh
#define rtt_units_FundUnit_hh

#include "UnitSystemEnums.hh"

namespace rtt_units {

//================================================================================================//
/*!
 * \class FundUnit
 * \brief Fundamental Unit Type
 *
 * \sa UnitSystem
 * \sa UnitSystemType
 */
//================================================================================================//

template <typename F> // T is one of { Ltype, Mtype, etc. }
class FundUnit        // Length, Mass, time, etc...
{
public:
  //! default constructor
  FundUnit(F const enumVal, double const *cf, std::string const &labels)
      : d_definingEnum(enumVal), d_cf(cf[d_definingEnum]),
        d_label(setUnitLabel(d_definingEnum, labels)) { /* empty */
  }

  // ACCESSORS

  //! return defining enumeration as specified in units/UnitSystemEnums.hh
  F enumVal() const { return d_definingEnum; }
  //! return conversion factor.
  //! Multiply by this number to get a value in SI units.
  double cf() const { return d_cf; }
  //! return unit label (e.g.: cm or keV)
  std::string label() const { return d_label; }

private:
  // DATA

  F d_definingEnum;
  double d_cf;
  std::string d_label;
};

} // end namespace rtt_units

#endif // rtt_units_FundUnit_hh

//------------------------------------------------------------------------------------------------//
// end of FundUnit.hh
//------------------------------------------------------------------------------------------------//
