//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/General_Octant_Quadrature.hh
 * \author Kelly Thompson
 * \date   Wed Sep  1 10:19:52 2004
 * \brief  A class to encapsulate a 3D Level Symmetric quadrature set.
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef quadrature_General_Octant_Quadrature_hh
#define quadrature_General_Octant_Quadrature_hh

#include "Octant_Quadrature.hh"

namespace rtt_quadrature {

//================================================================================================//
/*!
 * \class General_Octant_Quadrature
 * \brief A class to encapsulate a client-defined ordinate set.
 */
//================================================================================================//

class General_Octant_Quadrature : public Octant_Quadrature {
public:
  // CREATORS
  General_Octant_Quadrature(unsigned const sn_order, std::vector<double> const &mu,
                            std::vector<double> const &eta, std::vector<double> const &xi,
                            std::vector<double> const &wt, unsigned number_of_levels,
                            Quadrature_Class quadradure_class);

  General_Octant_Quadrature() = delete;

  // ACCESSORS

  std::vector<double> const &mu() const { return mu_; }
  std::vector<double> const &eta() const { return eta_; }
  std::vector<double> const &xi() const { return xi_; }
  std::vector<double> const &wt() const { return wt_; }

  // SERVICES

  // These functions override the virtual member functions specifed in the parent class Quadrature.

  std::string name() const override;
  std::string parse_name() const override;
  Quadrature_Class quadrature_class() const override;
  unsigned number_of_levels() const override;
  std::string as_text(std::string const &indent) const override;
  bool is_open_interval() const override;
  bool check_class_invariants() const;

  // STATICS

  static std::shared_ptr<Quadrature> parse(Token_Stream &tokens);

private:
  // IMPLEMENTATION

  //! Virtual hook for create_ordinate_set
  void create_octant_ordinates_(std::vector<double> &mu, std::vector<double> &eta,
                                std::vector<double> &wt) const override;

  // DATA
  std::vector<double> mu_, eta_, xi_, wt_;
  unsigned number_of_levels_;
  Quadrature_Class quadrature_class_;
  bool is_open_interval_;
};

} // end namespace rtt_quadrature

#endif // quadrature_General_Octant_Quadrature_hh

//------------------------------------------------------------------------------------------------//
// end of quadrature/General_Octant_Quadrature.hh
//------------------------------------------------------------------------------------------------//
