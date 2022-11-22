//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Quadrature_Interface.hh
 * \author Jae Chang
 * \date   Tue Jan 27 08:51:19 2004
 * \brief  Quadrature interface definitions
 * \note   Copyright (C) 2015-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//
#ifndef rtt_quadrature_Quadrature_Interface_hh
#define rtt_quadrature_Quadrature_Interface_hh

extern "C" {

//================================================================================================//
/*!
 * \class quadrature_data
 * \brief Flattened quadrature data used for communicating with Fortran routines.
 */
//================================================================================================//

struct quadrature_data {
  int dimension{0};
  int type{0};
  int order{0};
  int azimuthal_order{0};
  int geometry{0};
  double *mu{nullptr};
  double *eta{nullptr};
  double *xi{nullptr};
  double *weights{nullptr};
};

//! An extern "C" interface to default constructor
void init_quadrature(quadrature_data &quad);

//! Get quadrature data (e.g.: weights and cosines)
void get_quadrature(quadrature_data &quad);

//! Ensure quadrature data is meaningful
void check_quadrature_validity(const quadrature_data &quad);

} // end extern "C" block

#endif // quadrature_Quadrature_Interface_hh

//------------------------------------------------------------------------------------------------//
// end of quadrature/Quadrature_Interface.hh
//------------------------------------------------------------------------------------------------//
