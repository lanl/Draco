//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Ordinate_Set_Factory.hh
 * \author Allan Wollaber
 * \date   Mon Mar  7 10:42:56 EST 2016
 * \brief  Builds an Ordinate_Set using a Quadrature_Interface struct
 * \note   Copyright (C) 2016-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef quadrature_Ordinate_Set_Factory_hh
#define quadrature_Ordinate_Set_Factory_hh

#include "Ordinate_Set.hh"
#include "Quadrature_Interface.hh"
#include <memory>

namespace rtt_quadrature {

//================================================================================================//
/*!
 * \class Ordinate_Set_Factory
 *
 * \brief Factory class to create an Ordinate_Set from quadrature_data
 */
//================================================================================================//

class Ordinate_Set_Factory {
public:
  // CREATORS

  explicit Ordinate_Set_Factory(const quadrature_data &quad_in) : quad_(quad_in) {
    check_quadrature_validity(quad_in);
  }

  // SERVICES

  //! Returns a smart pointer to the Quadrature object
  std::shared_ptr<Ordinate_Set> get_Ordinate_Set() const;

private:
  // DATA

  // Ordinate set data
  const quadrature_data quad_;
};

} // end namespace rtt_quadrature

#endif // quadrature_Ordinate_Set_Factory

//------------------------------------------------------------------------------------------------//
// end of quadrature/Ordinate_Set_Factory.hh
//------------------------------------------------------------------------------------------------//
