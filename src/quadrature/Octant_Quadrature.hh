//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Octant_Quadrature.hh
 * \author Kent Budge
 * \date   Friday, Nov 30, 2012, 08:28 am
 * \brief  A class to encapsulate a 3D Level Symmetric quadrature set.
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef quadrature_Octant_Quadrature_hh
#define quadrature_Octant_Quadrature_hh

#include "Quadrature.hh"

namespace rtt_quadrature {

//================================================================================================//
/*!
 * \class Octant_Quadrature
 * \brief This is an abstract class representing all quadratures over the unit sphere.
 *
 * At present, all our unit sphere quadratures are symmetric in octants, though we will likely relax
 * this restriction in the future.
 *
 * For level quadratures, the levels must be in the xi direction cosine. The user may override the
 * default axis assignments when he constructs an Ordinate_Set or an Ordinate_Space from the
 * Octant_Quadrature.
 */
//================================================================================================//

class Octant_Quadrature : public Quadrature {
public:
  // CREATORS

  explicit Octant_Quadrature(unsigned const sn_order)
      : Quadrature(sn_order), has_axis_assignments_(false), mu_axis_(), eta_axis_() {}

  Octant_Quadrature(unsigned const sn_order, unsigned const mu_axis, unsigned const eta_axis)
      : Quadrature(sn_order), has_axis_assignments_(true), mu_axis_(mu_axis), eta_axis_(eta_axis) {}

  // SERVICES
  bool has_axis_assignments() const override;

protected:
  std::string as_text(std::string const &indent) const override = 0;

  // IMPLEMENTATION

  //! Virtual hook for create_ordinate_set
  virtual void create_octant_ordinates_(std::vector<double> &mu, std::vector<double> &eta,
                                        std::vector<double> &wt) const = 0;

  // STATICS

  static void parse(Token_Stream &tokens, bool &has_axis_assignments, unsigned &mu_axis,
                    unsigned &eta_axis);

private:
  // IMPLEMENTATION

  using Quadrature::create_ordinates_;

  //! Virtual hook for create_ordinates
  std::vector<Ordinate> create_ordinates_(unsigned dimension, Geometry /*unused*/, double norm,
                                          bool include_starting_directions,
                                          bool include_extra_directions) const override;

  //! Virtual hook for create_ordinate_set
  std::vector<Ordinate> create_ordinates_(unsigned dimension, Geometry /*unused*/, double norm,
                                          unsigned mu_axis, unsigned eta_axis,
                                          bool include_starting_directions,
                                          bool include_extra_directions) const override;

  // DATA

  bool has_axis_assignments_;
  unsigned mu_axis_, eta_axis_;
};

} // end namespace rtt_quadrature

#endif // quadrature_Octant_Quadrature_hh

//------------------------------------------------------------------------------------------------//
// end of quadrature/Octant_Quadrature.hh
//------------------------------------------------------------------------------------------------//
