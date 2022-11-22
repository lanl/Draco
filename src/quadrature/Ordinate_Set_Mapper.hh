//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   quadrature/Ordinate_Set_Mapper.hh
 * \author Allan Wollaber
 * \date   Mon Mar  7 10:42:56 EST 2016
 * \brief  Declarations for the class rtt_quadrature::Ordinate_Set_Mapper.
 * \note   Copyright (C) 2016-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef quadrature_OrdinateSetMapper_hh
#define quadrature_OrdinateSetMapper_hh

#include "Ordinate_Set.hh"

namespace rtt_quadrature {

//================================================================================================//
/*!
 * \class Ordinate_Set_Mapper
 *
 * \brief Provides services to map an angle ordinate onto an ordinate set
 */
//================================================================================================//

class Ordinate_Set_Mapper {
public:
  // ENUMERATIONS

  //! Ordering of ordinates.
  enum Interpolation_Type {
    //! Associates an angle with its nearest ordinate
    NEAREST_NEIGHBOR,

    //! Reallocates weight to nearest three ordinates (at most, must be > 0) Note that this uses
    //! inverse weighting according to 1-dot product so that the closest ordinate is strongly
    //! preferred.
    NEAREST_THREE,

    //! Currently unimplemented, but would use an interpolating function based on a bandwidth
    //! depending on the dot product
    KERNEL_DENSITY_ESTIMATOR
  };

  // CREATORS

  explicit Ordinate_Set_Mapper(const Ordinate_Set &os_in) : os_(os_in) {
    Ensure(check_class_invariants());
  }

  // SERVICES

  //! Simple check of integrity of the private data
  bool check_class_invariants() const;

  //! Maps an ordinate and weight into the ordinate set
  void map_angle_into_ordinates(const Ordinate &ord_in, const Interpolation_Type &interp_in,
                                std::vector<double> &weights_in) const;

private:
  // DATA

  // Ordinate set data
  const Ordinate_Set os_;

  // SERVICE CLASSES
  // -------------------------------------------------------------------------
  // A simple functor to be used in computing a bunch of 3D dot products between a given ordinate
  // and all the ordinates in a container.
  // -------------------------------------------------------------------------
  struct dot_product_functor_3D {
    explicit dot_product_functor_3D(const Ordinate &o_in) : o1(o_in) {}

    // Returns the 3D dot product of the ordinate passed into the functor with the local ordinate
    double operator()(const Ordinate &o2) const {
      return o1.mu() * o2.mu() + o1.eta() * o2.eta() + o1.xi() * o2.xi();
    }
    const Ordinate o1;
  };

  // -------------------------------------------------------------------------
  // A simple functor to be used in computing a bunch of 1D dot products between a given ordinate
  // and all the 1D ordinates in a container.
  // -------------------------------------------------------------------------
  struct dot_product_functor_1D {
    explicit dot_product_functor_1D(const Ordinate &o_in) : o1(o_in) {}

    /*! Returns the dot product of the ordinate passed into the functor with the local ordinate. For
     *  1-D we only have the cosine of the polar axis for each ordinate, \f$ \theta_1 \f$ and \f$
     *  \theta_2 \f$.  To obtain the cosine of the angle between them, we need to calculate \f$
     *  \cos(\theta) = \cos(\theta_2 - \theta_1) \f$.  Instead of using a relatively expensive
     *  acos() function, we use the simple identity, \f$ \cos(\theta_2
     *  - \theta_1) = \cos(\theta_1)\cos(\theta_2) +
     *  \sin(\theta_1)\sin(\theta_2) \, ,\f$ where \f$ \sin(\theta_1) =
     *  \sqrt{1-\mu_1^2} . \f$
     */
    double operator()(const Ordinate &o2) const {
      if (rtt_dsxx::soft_equiv(o1.mu(), o2.mu()))
        return 1.0;

      const double &mu1(o1.mu());
      const double &mu2(o2.mu());
      const double eta1(sqrt(1.0 - mu1 * mu1));
      const double eta2(sqrt(1.0 - mu2 * mu2));
      double mu_btwn = mu1 * mu2 + eta1 * eta2;

      Ensure(-1.0 <= mu_btwn && mu_btwn <= 1.0);

      return mu_btwn;
    }
    const Ordinate o1;
  };

  //! Simple function to integrate the 0th moment
  double zeroth_moment(const std::vector<double> &weights) const;
};

} // end namespace rtt_quadrature

#endif // quadrature_OrdinateSetMapper_hh

//------------------------------------------------------------------------------------------------//
// end of quadrature/OrdinateSetMapper.hh
//------------------------------------------------------------------------------------------------//
