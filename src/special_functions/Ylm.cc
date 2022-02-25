//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   special_functions/Ylm.cc
 * \author Kent Budge
 * \date   Tue Sep 21 09:20:10 2004
 * \brief  Implementation of Ylm
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Ylm.hh"
#include "Factorial.hh"
#include "units/PhysicalConstants.hh"
#include <cstdlib>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <iostream>

namespace rtt_sf {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute the spherical harmonic coefficient multiplied by the Associated Legendre
 *        Polynomial
 *
 * \f[
 * {\sqrt {\frac{{2l + 1}}{{4\pi }}\frac{{(l - m)!}}{{(l + m)!}}}
 * P_l^m (\cos \theta )}.
 * \f]
 *
 * \param l Used to specify the degree \f$ \ell \f$ of \f$ P_{\ell,k}(\theta) \f$.
 * \param k Used to specify the order \f$ k \f$ of \f$ P_{\ell,k}(\theta) \f$.
 * \param mu The cosine of the azimuthal angle, \f$ \mu = \cos \theta \f$.
 * \return The spherical harmonic coefficient multiplied by the Associated Legendre Polynomial, \f$
 *         c_{l,k}P_{l,k}(\mu) \f$.
 *
 * For details on what is being computed by this routine see the Draco documentation on sf_sph.
 *
 * This function simply wraps the GSL function \c gsl_sf_legendre_sphPlm computes the value of \f$
 * {\sqrt {\frac{{2l + 1}}{{4\pi }}\frac{{(l - m)!}}{{(l + m)!}}} P_l^m ( \mu )} \f$.
 *
 * The Condon-Shortley Phase, \f$ (-1)^m \f$ is provided by this definition of \f$ P_{l,k} \f$.
 *
 * \sa <a href="http://mathworld.wolfram.com/Condon-ShortleyPhase.html"> Condon-Shortley Phase</a>
 *
 * \pre \f$ l \ge 0 \f$
 * \pre \f$ l \ge k \f$
 * \pre \f$ \mu \in [-1, 1] \f$
 */
double cPlk(unsigned const l, unsigned const k, double const mu) {
  Require(k <= l);
  Require(mu >= -1.0);
  Require(mu <= 1.0);

  // This routine computes the normalized associated legendre polynomial \f$ \sqrt{(2l+1)/(4\pi)}
  // \sqrt{(l-m)!/(l+m)!} P_l^m(x) \f$ suitable for use in spherical harmonics.

  return gsl_sf_legendre_sphPlm(l, k, mu);
}
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute the spherical harmonic coefficient multiplied by the Associated Legendre
 *        Polynomial as specified by Morel's Galerkin Quadrature paper.
 *
 * \param l Used to specify the degree \f$ \ell \f$ of \f$ P_{\ell,k}(\theta) \f$.
 * \param k Used to specify the order \f$ k \f$ of \f$ P_{\ell,k}(\theta) \f$.
 * \param mu The cosine of the azimuthal angle, \f$ \mu = \cos \theta \f$.
 * \param sumwt normalizing constant
 * \return Morel's spherical harmonic coefficient multiplied by the Associated Legendre Polynomial,
 * \f$ c_{l,k}P_{l,k}(\mu) \f$.
 *
 * Computes the coefficient \f$ \frac{{2l + 1}}{{\sum\limits_m {w_m } }}\sqrt {\frac{{(l - \left | k
 * \right |)!}}{{(l - \left | k \right |)!}}}P_{l,\left | k \right |} (\cos \theta ) \f$
 *
 * For details on what is being computed by this routine see the Draco documentation on sf_sph.
 *
 * This function uses the GSL function \c gsl_sf_legendre_Plm to compute the value of \f$ P_l^k (
 * \mu ) \f$.
 *
 * The Condon-Shortley Phase, \f$ (-1)^m \f$ is provided by this definition of \f$ P_{l,k} \f$.
 *
 * \sa <a href="http://mathworld.wolfram.com/Condon-ShortleyPhase.html"> Condon-Shortley Phase</a>
 *
 * \pre \f$ l \ge 0 \f$
 * \pre \f$ l \ge k \f$
 * \pre \f$ k \ge 0 \f$
 * \pre \f$ \mu \in [-1, 1] \f$
 */
double cPlkGalerkin(unsigned const l, unsigned const k, double const mu, double const sumwt) {
  using std::sqrt;

  Require(k <= l);
  Require(mu >= -1.0);
  Require(mu <= 1.0);

  double coeff((2 * l + 1) / sumwt);
  coeff *= sqrt(1.0 * factorial_fraction(l - k, l + k)); // ff = (l-k)! / (l+k)!
  coeff *= gsl_sf_legendre_Plm(l, k, mu);

  return coeff;
}
//------------------------------------------------------------------------------------------------//
/*!\brief Compute the normalized spherical harmonic \f$ y_{l,k}(\theta,\phi) \f$
 *
 * \param l Used to specify the degree \f$ \ell \f$ of \f$ Y_{\ell,k}(\theta,\phi) \f$.
 *
 * \param k Used to specify the order \f$ k \f$ of \f$ Y_{\ell,k}(\theta,\phi) \f$.
 *
 * \param theta The polar (colatitudinal) coordinate in (0,\f$\pi\f$).
 *
 * \param phi The azimuthal (longitudinal) coordinate in (0,\f$ 2\pi\f$).
 *
 * \return The normalized spherical harmonic value, \f$ y_l^m(\theta,\phi) \f$.
 *
 * For details on what is being computed by this routine see the Draco documentation on sf_sph.
 *
 * The GSL function \c gsl_sf_legendre_sphPlm computes the value of \f$ {\sqrt {\frac{{2l +
 * 1}}{{4\pi }}\frac{{(l - m)!}}{{(l + m)!}}} P_l^m (\cos \theta )} \f$.
 *
 * \pre \f$ l \ge 0 \f$
 * \pre \f$\left|m\right| \le l\f$
 * \pre \f$ \theta \in [0,\pi] \f$
 * \pre \f$ \phi \in [0, 2\pi] \f$
 *
 *  \warning We cannot correct for sum of quadrature weights not equal to 4 PI, \f$
 *  \sum\limits_m{w_m} \ne 4\pi \f$, because this package lives below quadrature.  This adjustment
 *  will need to be done if the quadrature packages uses Ylm.
 */
double normalizedYlk(unsigned const l, int const k, double const theta, double const phi) {
  int const absk(std::abs(k));
  double const mu(std::cos(theta));

  // The constant and the Associated Legendre Polynomial.
  double const cP(cPlk(l, absk, mu));

  // for k>0 and odd, the sign will be negative.
  double sign(std::pow(-1.0, absk));

  // As noted on the \ref special_functions_overview for this package, we are interested in the real
  // portion of the spherical harmonics function.

  double result(-9999.0);

  if (k > 0)
    result = cP * std::sqrt(2.0) * std::sin(k * phi);
  else if (k == 0)
    result = cP;
  else /* if k<0 */
    result = sign * cP * std::sqrt(2.0) * std::cos(absk * phi);
  return result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute the real portion of the spherical harmonic \f$ Y_{l,k}(\theta,\phi) \f$
 *
 * \param l Used to specify the degree \f$ \ell \f$ of \f$ Y_{\ell,k}(\theta,\phi) \f$.
 * \param k Used to specify the order \f$ k \f$ of \f$ Y_{\ell,k}(\theta,\phi) \f$.
 * \param theta The polar (colatitudinal) coordinate in (0,\f$\pi\f$).
 * \param phi The azimuthal (longitudinal) coordinate in (0,\f$ 2\pi\f$).
 * \return The real portion of the spherical harmonic value, \f$ Y_l^m(\theta,\phi) \f$.
 *
 * For details on what is being computed by this routine see the Draco documentation on sf_sph.
 *
 * The GSL function \c gsl_sf_legendre_sphPlm computes the value of \f$ {\sqrt {\frac{{2l +
 * 1}}{{4\pi }}\frac{{(l - m)!}}{{(l + m)!}}} P_l^m (\cos \theta )} \f$.
 *
 * Note that there is a possible sign change for k<0,
 * \f[
 *  Y_{l,-k}(\theta,\phi) = (-1)^m Y^*_{l,k}(\theta,phi)
 * \f]
 *
 * \pre \f$ l \ge 0 \f$
 * \pre \f$\left|m\right| \le l\f$
 * \pre \f$ \theta \in [0,\pi] \f$
 * \pre \f$ \phi \in [0, 2\pi] \f$
 *
 *  \warning We cannot correct for sum of quadrature weights not equal to 4 PI, \f$
 *     \sum\limits_m{w_m} \ne 4\pi \f$, because this package lives below quadrature.  This
 *     adjustment will need to be done if the quadrature packages uses Ylm.
 */
double realYlk(unsigned const l, int const k, double const theta, double const phi) {
  int const absk(std::abs(k));
  double const mu(std::cos(theta));
  double sign(1.0);

  // The constant and the Associated Legendre Polynomial.
  double const cP(cPlk(l, absk, mu));

  // Adjust the sign.
  if (k < 0)
    sign = std::pow(-1.0, absk);

  // As noted on the \ref special_functions_overview for this package, we are interested in the real
  // portion of the spherical harmonics function.

  return sign * cP * std::cos(absk * phi);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute the complex portion of the spherical harmonic \f$ Y_{l,k}(\theta,\phi) \f$
 *
 * \param l Used to specify the degree \f$ \ell \f$ of \f$ Y_{\ell,k}(\theta,\phi) \f$.
 * \param k Used to specify the order \f$ k \f$ of \f$ Y_{\ell,k}(\theta,\phi) \f$.
 * \param theta The polar (colatitudinal) coordinate in (0,\f$\pi\f$).
 * \param phi The azimuthal (longitudinal) coordinate in (0,\f$ 2\pi\f$).
 * \return The complex portion of the spherical harmonic value, \f$ Y_l^m(\theta,\phi) \f$.
 *
 * For details on what is being computed by this routine see the Draco documentation on sf_sph.
 *
 * The GSL function \c gsl_sf_legendre_sphPlm computes the value of \f$ {\sqrt
 * {\frac{{2l + 1}}{{4\pi }}\frac{{(l - m)!}}{{(l + m)!}}} P_l^m (\cos \theta
 * )} \f$.
 *
 * \pre \f$ l \ge 0 \f$
 * \pre \f$\left|m\right| \le l\f$
 * \pre \f$ \theta \in [0,\pi] \f$
 * \pre \f$ \phi \in [0, 2\pi] \f$
 *
 *  \warning We cannot correct for sum of quadrature weights not equal to 4 PI, \f$
 *     \sum\limits_m{w_m} \ne 4\pi \f$, because this package lives below quadrature.  This
 *     adjustment will need to be done if the quadrature packages uses Ylm.
 */
double complexYlk(unsigned const l, int const k, double const theta, double const phi) {
  int const absk(std::abs(k));
  double const mu(std::cos(theta));
  double sign(1.0);

  // The constant and the Associated Legendre Polynomial.
  double const cP(cPlk(l, absk, mu));

  // Adjust the sign.
  if (k < 0)
    sign = std::pow(-1.0, absk);

  return sign * cP * std::sin(absk * phi);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Compute the spherical harmonic as used by Morel's Galerkin Quadrature paper.
 *
 * \param l Used to specify the degree \f$ \ell \f$ of \f$ Y_{\ell,k}(\theta,\phi) \f$.
 * \param k Used to specify the order \f$ k \f$ of \f$ Y_{\ell,k}(\theta,\phi) \f$.
 * \param mu The cosine of the polar (colatitudinal) coordinate in (-1,1).
 * \param phi The azimuthal (longitudinal) coordinate in (0,\f$ 2\pi\f$).
 * \param sumwt normalizing coefficient
 * \return The spherical harmonic value as specified in the detailed documentation.
 *
 * We use a special form of spherical harmonics for creating the moment-to-discrete and
 * discrete-to-moment matrixes.  Morel describes the functional form in his paper \e A \e Hybrid \e
 * Collocation-Galerkin-Sn \e Method \e for \e Solving \e the \e Boltzmann \e Transport \e
 * Equation," Nuclear Science and Engineering, \b 101, 72-87, 1989.
 *
 * The principal difference is that the \f$ \frac{2\ell+1}{4\pi} \f$ coefficient does not appear
 * inside of a square root.  This change from the normal spherical harmonic coefficient is required
 * in order to provide the correct magnitude for the quadrature weights.  Morel's formulation is
 *
 * \f[
 * \begin{array}{l}
 *  k > 0:\;\;Y_{l,k} (\theta ,\phi ) = \frac{{2l + 1}}{{\sum\limits_m {w_m } }}\sqrt {2\frac{{(l - k)!}}{{(l - k)!}}} P_{l,k} (\cos \theta )\cos (\phi ) \ \
 *  k = 0:\;\;Y_{l,0} (\theta ,\phi ) = \frac{{2l + 1}}{{\sum\limits_m {w_m } }}P_l (\cos \theta ) \ \
 *  k > 0:\;\;Y_{l,k} (\theta ,\phi ) = \frac{{2l + 1}}{{\sum\limits_m {w_m } }}\sqrt {2\frac{{(l - \left| k \right|)!}}{{(l - \left| k \right|)!}}} P_{l,\left| k \right|} (\cos \theta )\sin (\phi ) \ \
 *  \end{array}
 * \f]
 *
 * Morel's scheme only includes specific order and moment combinations.  Please see the reference
 * for a full description.
 *
 * \pre \f$ l \ge 0 \f$
 * \pre \f$\left|m\right| \le l\f$
 * \pre \f$ \theta \in [0,\pi] \f$
 * \pre \f$ \phi \in [0, 2\pi] \f$
 */
double galerkinYlk(unsigned const l, int const k, double const mu, double const phi,
                   double const sumwt) {
  int const absk(std::abs(k));

  // The constant and the Associated Legendre Polynomial.
  double Ylk(cPlkGalerkin(l, absk, mu, sumwt));

  // Adjust the sign.
  if (k < 0) {
    //        Ylk*=std::pow( -1.0, absk );
    Ylk *= std::sqrt(2.0) * std::sin(absk * phi);
  } else if (k > 0)
    Ylk *= std::sqrt(2.0) * std::cos(absk * phi);

  return Ylk;
}

//------------------------------------------------------------------------------------------------//
double Ylm(unsigned const l, int const m, double const mu, double const phi, double const sumwt) {
  int const absm(std::abs(m));

  // The constant and the Associated Legendre Polynomial.

  double const alpha = (2.0 * l + 1.0) * (gsl_sf_fact(l - absm) / gsl_sf_fact(l + absm)) *
                       (m != 0 ? 2.0 : 1.0) / sumwt;

  double ylm = sqrt(alpha) * gsl_sf_legendre_Plm(l, absm, mu);

  // Adjust
  if (m < 0) {
    ylm *= std::sin(absm * phi);
  } else if (m > 0) {
    ylm *= std::cos(absm * phi);
  }

  return ylm;
}

} // end namespace rtt_sf

//------------------------------------------------------------------------------------------------//
// end of Ylm.cc
//------------------------------------------------------------------------------------------------//
