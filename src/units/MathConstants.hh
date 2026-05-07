//--------------------------------------------*-C++-*---------------------------------------------//
/*! \file   units/MathConstants.hh
 *  \author Kelly Thompson
 *  \brief  Provide a single place where physical constants (pi, speed of light, etc) are defined.
 *  \note   Copyright (C) 2014-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_units_MathConstants_hh
#define rtt_units_MathConstants_hh

#include <numbers>

//! \namespace rtt_units Namespace for units and physical constants
namespace rtt_units {

//! pi the ratio of a circle's circumference to its diameter (dimensionless)
double constexpr PI = std::numbers::pi;

// Euler's number (dimensionless)
double constexpr N_EULER = std::numbers::e;

} // end namespace rtt_units

#endif // rtt_units_MathConstants_hh

//------------------------------------------------------------------------------------------------//
// end of units/MathConstants.hh
//------------------------------------------------------------------------------------------------//
