//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Ordinate_Space.hh
 * \author Kent Budge
 * \date   Mon Mar 26 16:11:19 2007
 * \brief  Definition of class Ordinate_Space
 * \note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef quadrature_Ordinate_Space_hh
#define quadrature_Ordinate_Space_hh

#include "Moment.hh"
#include "Ordinate_Set.hh"
#include "QIM.hh"
#include "Quadrature_Class.hh"

namespace rtt_quadrature {
using std::ostream;

//================================================================================================//
/*!
 * \class Ordinate_Space
  * \brief Describes a choice of discrete ordinate and truncated moment representations of ordinate
 *        space.
 *
 * This class encapsulates descriptions of a discrete ordinate space and a truncated moment space
 * and provides representations of a set of operators operating on these spaces. The moment space is
 * diagonal in any scattering operator for an isotropic material and thus is preferred for coupling
 * to the physics.
 *
 * The discrete ordinate space is described by the methods and data inherited from Ordinate_Set. The
 * heart of this description is the vector<Ordinate> returned by the Ordinate_Set::ordinates method.
 *
 * The moment space is described by the moment to discrete transformation matrix returned by
 * Ordinate_Space::M() and the discrete to moment transformation matrix returned by
 * Ordinate_Space::D(). The moment rank of these matrices is returned by
 * Ordinate_Space::number_of_ordinates while the ordinate rank is given by
 * Ordinate_Set::ordinates().size(). The number of moments of each order in the moment space is
 * returned as a vector<unsigned> by Ordinate_Space::moments_per_order(), and so the size of
 * Ordinate_Space::moments_per_order() is one more than the moment expansion order. The actual L and
 * M of the spherical harmonic corresponding to each moment is returned by
 * Ordinate_Space::moments().
 *
 * There is a subtlety here: Ordinate_Space::moments.size() is greater than
 * Ordinate_Space::number_of_moments() because the former may contain additional moments, of higher
 * order than the specified expansion order, used to construct the M and D matrices but not included
 * in the actual scattering expansion. See Galerkin_Ordinate_Space for an example where this is done
 * with an explanation of why it is useful.
 *
 * The actual spherical harmonic basis used for the moment space is given by rtt_sf::Ylm, which is a
 * real representation in which m=0 is the axially symmetric moment, m<0 are the moments even in
 * phi, and m>0 are the moments odd in phi. For 1-D not axisymmetric, we do the obvious thing and
 * choose to align the polar axis with the coordinate axis. However, for all other geometries, we
 * assign mu = cos(theta)*cos(phi) to the first coordinate axis and xi=sin(theta) to the second
 * coordinate axis. This may seem a strange choice, but it simplifies the representation of
 * symmetries in reduced geometry, particularly for Galerkin_Ordinate_Space.
 *
 * The mu, eta, and xi reflection maps give, for each ordinate i, the index of the ordinate that is
 * the reflection of i in the specified coordinate plane. Thus, on a reflection plane reflecting the
 * first coordinate, the specific intensity of ordinate i is reflected into the specific intensity
 * of ordinate reflec_mu[i]. This greatly simplifies implementing reflection boundary conditions.
 *
 * In curvilinear geometry, the streaming operator includes a nontrivial angle derivative ithat
 * introduces dependencies between ordinates. We assume that an angle derivative can be cast in
 * block bidiagonal form, so that there is not more than one direct dependency per ordinate. The
 * Ordinate_Space may then order the ordinates by dependency, so the first ordinate can have no
 * dependencies, the second may be directly dependent only on the first, and so on. Thus a client
 * need only check whether an ordinate is dependent on the preceeding angle or not.
 *
 * Each of the blocks in the block bidiagonal form of the angle operator is referred to as a
 * "level." This is terminology held over from the particular case of 2-D axisymmetric geometry,
 * where the ordinate sets generally are organized on "levels" having the same z direction cosine
 * (xi) which are coupled by the omega derivative term. It is useful to know the number of such
 * blocks to optimize storage of intermediate results.
 *
 * We illustrate with an example. In axisymmetric geometry, Morel's discretization (J.E. Morel, "An
 * R-Z Geometry Triangular-Mesh Sn Spatial Differencing Scheme," Technical Memorandum
 * CCS-4:03-14(U), May 12, 2003.) of the angle derivative term in the streaming operator is
 *
 * \f$\frac{\partial (\eta \psi)}{\partial \omega}\approx
 * \frac{1}{w_m}(\alpha_{m+1/2}\psi_{m+1/2}-\alpha_{m-1/2}\psi_{m-1/2})\f$
 *
 * where
 *
 * \f$\psi_{m+1/2} = \frac{1}{\tau_m}(\psi_m-(1-\tau_m)\psi_{m-1/2})\f$
 *
 * Thus we define
 *
 * <code>Psi_Coefficient(m)</code> = \f$P_m =
 * \frac{\alpha_{m+1/2}}{w_m\tau_m}\f$
 *
 * <code>Source_Coefficient(m)</code> =
 \f$S_m = \frac{\alpha_{m+1/2}\frac{1-\tau_m}{\tau_m}+\alpha_{m-1/2}}{w_m}\f$
 *
 * <code>Bookkeeping_Coefficient(m)</code> = \f$B_m = \frac{1}{\tau_m}\f$
 *
 * The angle derivative can then be coded as
 *
 * \f$\frac{\partial (\eta \psi)}{\partial \omega} =
 * P_m\psi_m-S_m\psi_{m-1/2}\f$
 *
 * and the next midpoint intensity as
 *
 * \f$\psi_{m+1/2} = B_m\psi_m-(1-B_m)\psi_{m-1/2})\f$
 *
 * Similar expressions can be written for spherical geometry. The Ordinate_Space interface hides
 * these details, presenting only the \f$P_m\f$, \f$S_m\f$, and \f$B_m\f$ coefficients required for
 * actual computation.
 *
 * Note that this discretization of the angle derivative terms must still be substituted into the
 * transport equation, which is then further discretized in space. Thus the angle derivative term
 * will generally be multiplied by an additional factor arising from the spatial discretization.
 */
//================================================================================================//

class Ordinate_Space : public rtt_quadrature::Ordinate_Set {
public:
  // CREATORS

  //! Specify the ordinate quadrature with defaults.
  Ordinate_Space(unsigned dimension, rtt_mesh_element::Geometry geometry,
                 std::vector<Ordinate> const &ordinates, int expansion_order,
                 bool extra_starting_directions = false, Ordering ordering = LEVEL_ORDERED);

  // ACCESSORS

  int expansion_order() const { return expansion_order_; }

  bool has_extra_starting_directions() const { return has_extra_starting_directions_; }

  unsigned number_of_levels() const { return number_of_levels_; }

  std::vector<unsigned> const &levels() const { return levels_; }

  //! Return the angle index for the most positively outward-directed angle on every level.
  std::vector<unsigned> const &first_angles() const { return first_angles_; }

  //! Is an ordinate on the same level as the preceding ordinate?
  bool is_dependent(unsigned const ordinate) const {
    Require(ordinate < ordinates().size());

    return is_dependent_[ordinate];
  }

  //! Return \f$\alpha_{m+1/2}\f$ for ordinate \f$m\f$
  std::vector<double> const &alpha() const { return alpha_; }

  //! Return \f$P\tau_m\f$ for ordinate \f$m\f$
  std::vector<double> const &tau() const { return tau_; }

  //! Return \f$P_m\f$ for ordinate \f$m\f$
  double psi_coefficient(unsigned ordinate_index) const;

  //! Return \f$S_m\f$ for ordinate \f$m\f$
  double source_coefficient(unsigned ordinate_index) const;

  //! Return \f$B_m\f$ for ordinate \f$m\f$
  double bookkeeping_coefficient(unsigned ordinate_index) const;

  unsigned number_of_moments() const { return number_of_moments_; }

  //! Return the moment descriptions of the moment space.
  std::vector<Moment> const &moments() const { return moments_; }

  //! Return vector containing the number of moments for each L
  std::vector<unsigned> const &moments_per_order() const { return moments_per_order_; }

  //! Return mu reflection map
  std::vector<unsigned> const &reflect_mu() const { return reflect_mu_; }

  //! Return eta reflection map
  std::vector<unsigned> const &reflect_eta() const { return reflect_eta_; }

  //! Return xi reflection map
  std::vector<unsigned> const &reflect_xi() const { return reflect_xi_; }

  bool check_class_invariants() const;

  // SERVICES

  //! What was the quadrature interpolation model?
  virtual QIM quadrature_interpolation_model() const = 0;

  //! Return the discrete to moment transform matrix
  virtual std::vector<double> D() const = 0;

  //! Return the moment to discrete transform matrix
  virtual std::vector<double> M() const = 0;

  //! Should the moment space be pruned to the specified order?
  virtual bool prune() const {
    return true;
    // By default, prune any moments beyond the user-specified expansion order.  Overridden by
    // Galerkin_Ordinate_Space::prune().
  }

  //! Return the scattering moment to flux map.
  void moment_to_flux(std::array<unsigned, 3> &flux_map, std::array<double, 3> &flux_fact) const;

  //! Return the flux to scattering moment map.
  void flux_to_moment(std::array<unsigned, 3> &flux_map, std::array<double, 3> &flux_fact) const;

  // STATICS

  double compute_azimuthalAngle(double mu, double eta);

protected:
  // IMPLEMENTATION

  void compute_moments_(Quadrature_Class quadrature_class, int sn_order);

  std::vector<Moment> compute_n2lk_(Quadrature_Class quadrature_class, unsigned sn_order);

  virtual std::vector<Moment> compute_n2lk_1D_(Quadrature_Class, unsigned sn_order) = 0;
  virtual std::vector<Moment> compute_n2lk_1Da_(Quadrature_Class, unsigned sn_order) = 0;
  virtual std::vector<Moment> compute_n2lk_2D_(Quadrature_Class, unsigned sn_order) = 0;
  virtual std::vector<Moment> compute_n2lk_2Da_(Quadrature_Class, unsigned sn_order) = 0;
  virtual std::vector<Moment> compute_n2lk_3D_(Quadrature_Class, unsigned sn_order) = 0;

private:
  // IMPLEMENTATION

  void compute_angle_operator_coefficients_();

  void compute_reflection_maps_();

  // DATA

  int expansion_order_;
  bool has_extra_starting_directions_;
  unsigned number_of_levels_;
  std::vector<unsigned> levels_;
  std::vector<unsigned> first_angles_;
  std::vector<bool> is_dependent_; //!< Is an ordinate dependent on the preceding ordinate?

  //! Reflection maps
  std::vector<unsigned> reflect_mu_, reflect_eta_, reflect_xi_;

  /*! Coefficients for angle derivative terms.  These are defined in Morel's research note of 12 May
   * 2003 for axisymmetric geometry. */
  std::vector<double> alpha_;
  std::vector<double> tau_;

  //! Number of moments up to the expansion order. Does not include Galerkin augments. This is the
  //! moment rank of the M and D matrices.
  unsigned number_of_moments_;
  //! Moments of moment expansion. Includes Galerkin augments, if any.
  std::vector<Moment> moments_;
  //! Moments per order. Does not include Galerkin augments.
  std::vector<unsigned> moments_per_order_;
};

} // end namespace rtt_quadrature

#endif // quadrature_Ordinate_Space_hh

//---------------------------------------------------------------------------------------//
// end of quadrature/Ordinate_Space.hh
//---------------------------------------------------------------------------------------//
