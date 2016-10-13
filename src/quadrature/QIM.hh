//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/QIM.hh
 * \author Kent Budge
 * \date   Mon Mar 26 16:11:19 2007
 * \brief  Definition of QIM enumeration
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------------------//
// $Id: QIM.hh 6718 2012-08-30 20:03:01Z warsa $
//---------------------------------------------------------------------------------------//

#ifndef quadrature_QIM_hh
#define quadrature_QIM_hh

#include "parser/Token_Stream.hh"

namespace rtt_quadrature {
using rtt_parser::Token_Stream;

//! Quadrature Interpolation Model: specifies how to compute the
//! Discrete-to-Moment operator.

enum QIM {
  SN,  /*!< Use the standard SN method. */
  GQ1, /*!< Use Morel's Galerkin Quadrature method. */
  GQ2, /*!< Use Warsa/Prinja Galerkin Quadrature method. */
  GQF, /*!< Use Morel's Galerkin Quadrature method and retain all moments. */
  SVD, /*!< Let M be an approximate inverse of D. */

  END_QIM //!< Sentinel value
};

DLL_PUBLIC_quadrature void parse_quadrature_interpolation_model(Token_Stream &,
                                                                QIM &);

std::string quadrature_interpolation_model_as_text(QIM);

} // end namespace rtt_quadrature

#endif // quadrature_QIM_hh

//---------------------------------------------------------------------------------------//
// end of quadrature/QIM.hh
//---------------------------------------------------------------------------------------//
