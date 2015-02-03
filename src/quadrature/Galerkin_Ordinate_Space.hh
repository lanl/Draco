//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/Galerkin_Ordinate_Space.hh
 * \author Kent Budge
 * \date   Mon Mar 26 16:11:19 2007
 * \brief  Definition of class Galerkin_Ordinate_Space
 * \note   Copyright (C) 2006-2015 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------------------//
// $Id: Galerkin_Ordinate_Space.hh 6718 2012-08-30 20:03:01Z warsa $
//---------------------------------------------------------------------------------------//

#ifndef quadrature_Galerkin_Ordinate_Space_hh
#define quadrature_Galerkin_Ordinate_Space_hh

#include "Ordinate_Space.hh"

namespace rtt_quadrature
{
using std::ostream;

//=======================================================================================//
/*!
 * \class Galerkin_Ordinate_Space
 * \brief Represents ordinate operators for a Galerkin moment space. 
 *
 * The moment space contains all moments (that are not identically zero due to
 * symmetry) up to the specified scattering order, but the moment to discrete
 * operator M and discrete to moment operator D are computed as if enough additional
 * higher moments are included in the moment space to make D and M square. The
 * higher moment terms are then discarded, but the non-square D and M retain
 * the property that DM is the identity. This stabilizes the moment to
 * discrete and discrete to moment operations at high scattering orders.
 *
 * When the additional moments are added, the SN quadrature order is provided,
 * and additional moments added based on the assumption that we are using
 * triangular quadrature sets. If an expansion order L < N-1 is
 * requested, both D and M wlll be appropriately truncated. If an expansion
 * order L > N-1 is requested, the computation cannot proceed and an exception
 * thrown.
 *
 * \todo For triangular sets, we assume L = N-1. The other two cases, L > N-1
 * and L < N-1 need to be handled.
 *
 * \todo The Galerkin quadrature is only currently implemented for triangular
 * quadratures, such that it will be necessary to extend the method to include
 * square (product) quadratures.
 *
 */
//=======================================================================================//

class Galerkin_Ordinate_Space : public Ordinate_Space
{
  public:

    // NESTED CLASSES AND TYPEDEFS

    // CREATORS

    //! Specify the ordinate quadrature with defaults.
    Galerkin_Ordinate_Space(unsigned dimension,
                            Geometry geometry,
                            vector<Ordinate> const &ordinates,
                            Quadrature_Class quadrature_class,
                            unsigned sn_order,
                            unsigned expansion_order,
                            QIM const method,
                            bool extra_starting_directions=false,
                            Ordering ordering=LEVEL_ORDERED);

    // MANIPULATORS

    // ACCESSORS

    bool check_class_invariants() const;

    // SERVICES
    
    virtual QIM quadrature_interpolation_model() const;

    //! Return the discrete to moment transform matrix
    virtual vector<double> D() const;

    //! Return the moment to discrete transform matrix
    virtual vector<double> M() const;

    bool prune() const { return method_ != GQF; }

    // STATICS


  private:

    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION

    virtual vector<Moment> compute_n2lk_1D_(Quadrature_Class,
                                            unsigned sn_order);
    
    virtual vector<Moment> compute_n2lk_1Da_(Quadrature_Class,
                                             unsigned sn_order);
    
    virtual vector<Moment> compute_n2lk_2D_(Quadrature_Class,
                                            unsigned sn_order);

    virtual vector<Moment> compute_n2lk_2Da_(Quadrature_Class,
                                             unsigned sn_order);
    
    virtual vector<Moment> compute_n2lk_3D_(Quadrature_Class,
                                            unsigned sn_order);

  private:

    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION
    
    void compute_operators();

    vector<double>  compute_M_SN(vector<Ordinate> const &ordinates);
    vector<double>  compute_D_SN(vector<Ordinate> const &ordinates,
                                 vector<double> const &Min);

    vector<double> compute_inverse(unsigned const m,
                                   unsigned const n,
                                   vector<double> const &Ain);

    vector<double> augment_D(vector<unsigned> const &indexes,
                             unsigned const numCartesianOrdinates,
                             vector<double> const &D);
    
    vector<double> augment_M(vector<unsigned> const &indexes,
                             vector<double> const &M);
    // DATA

    QIM const method_;
    
    //! Discrete to moment matrix
    vector<double> D_;
    //! Moment to discrete matrix
    vector<double> M_;
};

} // end namespace rtt_quadrature

#endif // quadrature_Galerkin_Ordinate_Space_hh

//---------------------------------------------------------------------------------------//
//              end of quadrature/Galerkin_Ordinate_Space.hh
//---------------------------------------------------------------------------------------//
