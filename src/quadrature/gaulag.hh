//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/gaulag.hh
 * \author Kent Budge
 * \date   Tue Sep 14 13:16:09 2004
 * \brief  Gauss-Laguerre quadrature
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_quadrature_gaulag_hh
#define rtt_quadrature_gaulag_hh

#include <limits>

#include <gsl/gsl_sf_gamma.h>

#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"
//#include "gammln.hh"

namespace rtt_utils
{
//---------------------------------------------------------------------------//
/*! 
 * \brief Gauss-Laguerre quadrature
 *
 * Calculate abscissae and weights for Gauss-Laguerre quadrature:
 * \f$\int_{0}^{\infty}{x^\alpha e^{-x} f(x) dx} = \sum_{j=0}^{N-1}w_j f(x_j)\f$
 * 
 * \arg \a FieldVector A random access container on a field type.
 *
 * \param x On return, contains abscissae \f$x_j\f$ for quadrature.
 * \param w On return, contains weights \f$w_j\f$ for quadrature.
 * \param alf Lagrange power parameter.
 * \param N Number of points in quadrature.
 * 
 */
template<class FieldVector>
void gaulag(FieldVector &x,
	    FieldVector &w,
	    double const alf,
	    unsigned const n)
{
    using namespace std;

    typedef typename FieldVector::value_type Field;

    const unsigned MAXITS=10;

    x.resize(n);
    w.resize(n);

    Field z;
    for (unsigned i=0; i<n; ++i)
    {
	if (i==0)
	{
	    z = (1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
	}
	else if (i==1)
	{
	    z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
	}
	else
	{
	    const Field ai = i-1;
	    z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
		  (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
	}
	unsigned its;
	Field pp, p2;
	for (its=0; its<MAXITS; ++its)
	{
	    Field p1 = 1.0;
	    p2 = 0.0;
	    for (unsigned j=0; j<n; ++j)
	    {
		const Field p3=p2;
		p2 = p1;
		p1 = ((2*j+1+alf-z)*p2 - (j+alf)*p3)/(j+1);
	    }
	    
	    // p1 is now the desired Legendre polynomial evaluated at z. We
	    // next compute pp, its derivative, by a standard relation
	    // involving also p2, the polynomial of one lower order.
	    pp = (n*p1-(n+alf)*p2)/z;

	    const Field z1 = z;

            // Newton's Method
            z = z1-p1/pp;
	    if (fabs(z-z1) < 100*numeric_limits<Field>::epsilon()) break;
	}
	if (its >= MAXITS) 
	{
	    throw std::range_error(std::string("too many iterations in gaulag"));
	}
	x[i] = z;
	w[i] = -exp( gsl_sf_lngamma(alf+n)
                     - gsl_sf_lngamma( static_cast<double>(n) ) )
               / (pp*n*p2);
    }
}

} // end namespace rtt_utils

#endif // utils_gaulag_hh

//---------------------------------------------------------------------------//
//              end of utils/gaulag.hh
//---------------------------------------------------------------------------//
