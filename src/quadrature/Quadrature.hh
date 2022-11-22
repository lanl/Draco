//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Quadrature.hh
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  Quadrature class header file.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_quadrature_Quadrature_hh
#define rtt_quadrature_Quadrature_hh

#include "Ordinate_Space.hh"

namespace rtt_quadrature {

//================================================================================================//
/*!
 * \class Quadrature
 *
 * \brief A class to encapsulate the angular discretization.
 *
 * The Quadrature class is an abstraction representation of an angular quadrature scheme. It can be
 * used to generate an Ordinate_Set containing the correct set of ordinate directions for a given
 * geometry. It can also be used to generate an Ordinate_Space describing both a discrete ordinate
 * and a truncated moment representation of ordinate space.
 *
 * All multi-dimensional quadratures (not interval quadratures) are expected to align any level
 * structure they may have with the xi direction cosine.
 *
 * When an Ordinate_Set is constructed from the Quadrature, the direction cosines in the Quadrature
 * must be mapped to coordinate axes in the problem geometry. By default, 1-D non-axisymmetric maps
 * xi to the coordinate axis; 1-D axisymmetric maps mu to the coordinate axis and xi to the
 * (non-represented) symmetry axis; 2-D maps mu to the first coordinate axis and xi to the second
 * coordinate axis; and 3-D maps mu to the first, eta to the second, and xi to the third coordinate
 * axes. This is to ensure that the levels are placed on the axis of symmetry in reduced geometries.
 *
 * The client may override these default assignments. However, if he assigns any direction cosine
 * other than xi to the axis of symmetry in axisymmetric geometry, Bad Things Will Happen with any
 * supported quadrature except Level_Symmetric (for which axis assignment is without effect anyway.)
 */
//================================================================================================//
class Quadrature {
public:
  // ENUMERATIONS AND TYPEDEFS

  using Geometry = rtt_mesh_element::Geometry;

  // CREATORS
  explicit Quadrature(unsigned const sn_order) : sn_order_(sn_order) { /* empty */
  }

  //! Virtual destructor.
  virtual ~Quadrature() = default;

  Quadrature(Quadrature const &rhs) = delete;
  Quadrature(Quadrature &&rhs) noexcept = delete;
  Quadrature &operator=(Quadrature const &rhs) = delete;
  Quadrature &operator=(Quadrature &&rhs) noexcept = delete;

  // ACCESSORS
  unsigned sn_order() const { return sn_order_; }

  // SERVICES

  //! Returns a string containing the name of the quadrature set.
  virtual std::string name() const = 0;

  //! Returns a string containing the input deck name of the set.
  virtual std::string parse_name() const = 0;

  //! Is this an interval or octant (1-D or multi-D) quadrature?
  virtual Quadrature_Class quadrature_class() const = 0;

  //! Number of level sets. A value of 0 indicates this is not a level set quadrature.
  virtual unsigned number_of_levels() const = 0;

  //! Produce a text representation of the object
  virtual std::string as_text(std::string const &indent) const = 0;

  //! Are the axes assigned?
  virtual bool has_axis_assignments() const = 0;

  //! Is the quadrature an open interval quadrature?
  virtual bool is_open_interval() const;

  std::vector<Ordinate> create_ordinates(unsigned dimension, Geometry geometry, double norm,
                                         unsigned mu_axis, unsigned eta_axis,
                                         bool include_starting_directions,
                                         bool include_extra_directions) const;

  std::vector<Ordinate> create_ordinates(unsigned dimension, Geometry geometry, double norm,
                                         bool include_starting_directions,
                                         bool include_extra_directions) const;

  std::shared_ptr<Ordinate_Set> create_ordinate_set(unsigned dimension, Geometry geometry,
                                                    double norm, bool include_starting_directions,
                                                    bool include_extra_directions,
                                                    Ordinate_Set::Ordering ordering) const;

  std::shared_ptr<Ordinate_Space> create_ordinate_space(unsigned dimension, Geometry geometry,
                                                        int moment_expansion_order,
                                                        bool include_extra_directions,
                                                        Ordinate_Set::Ordering ordering,
                                                        QIM qim) const;

  std::shared_ptr<Ordinate_Space>
  create_ordinate_space(unsigned dimension, Geometry geometry, unsigned moment_expansion_order,
                        unsigned mu_axis, unsigned eta_axis, bool include_extra_directions,
                        Ordinate_Set::Ordering ordering, QIM qim) const;

protected:
  // IMPLEMENTATION

  void add_1D_starting_directions_(Geometry geometry, bool add_starting_directions,
                                   bool add_extra_starting_directions,
                                   std::vector<Ordinate> &ordinates) const;

  void add_2D_starting_directions_(Geometry geometry, bool add_starting_directions,
                                   bool add_extra_starting_directions,
                                   std::vector<Ordinate> &ordinates) const;

  void map_axes_(unsigned mu_axis, unsigned eta_axis, std::vector<double> &mu,
                 std::vector<double> &eta, std::vector<double> &xi) const;

  //! Virtual hook for create_ordinates
  virtual std::vector<Ordinate> create_ordinates_(unsigned dimension, Geometry geometry,
                                                  double norm, unsigned mu_axis, unsigned eta_axis,
                                                  bool include_starting_directions,
                                                  bool include_extra_directions) const = 0;

  //! Virtual hook for create_ordinates
  virtual std::vector<Ordinate> create_ordinates_(unsigned dimension, Geometry geometry,
                                                  double norm, bool include_starting_directions,
                                                  bool include_extra_directions) const = 0;

  // data

  unsigned const sn_order_;
};

} // end namespace rtt_quadrature

#endif // rtt_quadrature_Quadrature_hh

//------------------------------------------------------------------------------------------------//
// end of quadrature/Quadrature.hh
//------------------------------------------------------------------------------------------------//
