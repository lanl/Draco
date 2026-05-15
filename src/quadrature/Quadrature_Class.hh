//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Quadrature_Class.hh
 * \author Kent G. Budge
 * \brief  Define Quadrature_Class enumeration
 * \note   Copyright (C) 2012-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#pragma once
namespace rtt_quadrature {

enum Quadrature_Class : uint8_t {
  INTERVAL_QUADRATURE, //!< 1-D quadratures

  TRIANGLE_QUADRATURE, //!< 3-D triangular quadrature
  SQUARE_QUADRATURE,   //!< 3-D square quadrature
  OCTANT_QUADRATURE,   //!< 3-D octant quadrature, not triangular nor square

  END_QUADRATURE
};

} // end namespace rtt_quadrature

//------------------------------------------------------------------------------------------------//
// end of quadrature/Quadrature_Class.hh
//------------------------------------------------------------------------------------------------//
