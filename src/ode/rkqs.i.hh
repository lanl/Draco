//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ode/rkqs.i.hh
 * \author Kent Budge
 * \date   Mon Sep 20 15:15:40 2004
 * \brief  Integrate an ordinary differential equation with local error
 *         control using fifth-order Cash-Karp Runge-Kutta steps.
 * \note   Copyright (C) 2004-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef ode_rkqs_i_hh
#define ode_rkqs_i_hh

#include "Function_Traits.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include <algorithm>
#include <cmath>

namespace rtt_ode {

//------------------------------------------------------------------------------------------------//
/*! Perform a single fifth-order Cash-Karp Runge-Kutta step.
 *
 * \arg \a Field A field type, such as \c double or \c std::complex<double>.

 * \arg \a Function A function type that maps \f$ N+1 \f$ variables to \f$ N \f$
 *                  variables, supporting <code>operator()(double, const
 *                  FieldVector &, FieldVector &)</code>.
 *
 * \param y Values of dependent variables at start of step.
 * \param dydx Values of derivatives of dependent variables at start of step.
 * \param x Dependent variable at start of step.
 * \param h Step size.
 * \param yout On return, contains values of dependent variables at end of step.
 * \param yerr On return, contains estimated truncation error.
 * \param derivs Function for computing derivatives of dependent variables.
 *
 * \pre \c dydx.size()==y.size()
 * \post \c yout.size()==y.size()
 * \post \c yerr.size()==y.size();
 */
template <typename Field, typename Function>
void rkck(std::vector<Field> const &y, std::vector<Field> const &dydx, double x, double h,
          std::vector<Field> &yout, std::vector<Field> &yerr, Function derivs) {
  Require(dydx.size() == y.size());

  using std::vector;

  static double const a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875, b21 = 0.2,
                      b31 = 3.0 / 40.0, b32 = 9.0 / 40.0, b41 = 0.3, b42 = -0.9, b43 = 1.2,
                      b51 = -11.0 / 54.0, b52 = 2.5, b53 = -70.0 / 27.0, b54 = 35.0 / 27.0,
                      b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0, b63 = 575.0 / 13824.0,
                      b64 = 44275.0 / 110592.0, b65 = 253.0 / 4096.0, c1 = 37.0 / 378.0,
                      c3 = 250.0 / 621.0, c4 = 125.0 / 594.0, c6 = 512.0 / 1771.0,
                      dc1 = c1 - 2825.0 / 27648.0, dc3 = c3 - 18575.0 / 48384.0,
                      dc4 = c4 - 13525.0 / 55296.0, dc5 = -277.0 / 14336.0, dc6 = c6 - 0.25;

  Check(y.size() < UINT_MAX);
  const auto n = static_cast<unsigned>(y.size());

  yout.resize(n);
  yerr.resize(n);

  vector<Field> ytemp(n);
  vector<Field> ak2(n);
  vector<Field> ak3(n);
  vector<Field> ak4(n);
  vector<Field> ak5(n);
  vector<Field> ak6(n);

  for (unsigned i = 0; i < n; i++) {
    ytemp[i] = y[i] + b21 * h * dydx[i];
  }
  derivs(x + a2 * h, ytemp, ak2);
  for (unsigned i = 0; i < n; i++) {
    ytemp[i] = y[i] + h * (b31 * dydx[i] + b32 * ak2[i]);
  }
  derivs(x + a3 * h, ytemp, ak3);
  for (unsigned i = 0; i < n; i++) {
    ytemp[i] = y[i] + h * (b41 * dydx[i] + b42 * ak2[i] + b43 * ak3[i]);
  }
  derivs(x + a4 * h, ytemp, ak4);
  for (unsigned i = 0; i < n; i++) {
    ytemp[i] = y[i] + h * (b51 * dydx[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i]);
  }
  derivs(x + a5 * h, ytemp, ak5);
  for (unsigned i = 0; i < n; i++) {
    ytemp[i] =
        y[i] + h * (b61 * dydx[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i] + b65 * ak5[i]);
  }
  derivs(x + a6 * h, ytemp, ak6);
  for (unsigned i = 0; i < n; i++) {
    yout[i] = y[i] + h * (c1 * dydx[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);
  }
  for (unsigned i = 0; i < n; i++) {
    yerr[i] = h * (dc1 * dydx[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i]);
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Integrate an ordinary differential equation with local error
 * control using fifth-order Cash-Karp Runge-Kutta steps.
 *
 * \arg \a Field A field type.
 * \arg \a Function A function type that maps \f$N+1\f$ variables to \f$N\f$
 * variables, supporting <code>operator()(double, const FieldVector &,
 * FieldVector &)</code>.
 *
 * \param y
 * Dependent variables at the start of the step.  On return, contains the
 * dependent variables at the end of the step.
 * \param dydx
 * Derivatives of the dependent variables at the start of the step.
 * \param x
 * The independent variable.
 * \param htry
 * The step size to attempt.
 * \param eps
 * The error tolerance for the integration.
 * \param yscal
 * Points to scaling factors for the error tolerance.
 * \param hdid
 * On return, contains the step actually taken.
 * \param hnext
 * On return, contains the recommended next step.
 * \param derivs
 * Function to calculate derivatives of the independent variables.
 *
 * \pre <code> yscal.size()==y.size() </code>
 * \pre <code> dydx.size()==y.size()</code>
 * \pre \c eps>=0
 *
 * \post \c dydx.size()==y.size()
 */

template <typename Field, typename Function>
void rkqs(std::vector<Field> &y, std::vector<Field> const &dydx, double &x, double htry, double eps,
          std::vector<Field> const &yscal, double &hdid, double &hnext, Function derivs) {
  Require(yscal.size() == y.size());
  Require(dydx.size() == y.size());
  Require(eps >= 0);

  using std::vector;

  Check(y.size() < UINT_MAX);
  const auto n = static_cast<unsigned>(y.size());

  double const SAFETY = 0.9;
  double const PSHRNK = -0.25;
  double const PGROW = -0.2;
  double const ERRCON = 1.89e-4;

  vector<double> yerr(n);
  vector<double> ytemp(n);

  double errmax;
  double h = htry;
  for (;;) {
    rkck(y, dydx, x, h, ytemp, yerr, derivs);
    errmax = 0.0;
    for (unsigned i = 0; i < n; i++) {
      errmax = std::max(errmax, std::abs(yerr[i] / yscal[i]));
    }
    errmax /= eps;
    if (errmax <= 1.0)
      break;
    double htemp = SAFETY * h * std::pow(errmax, PSHRNK);
    h = (h >= 0 ? std::max(htemp, 0.1 * h) : std::min(htemp, 0.1 * h));
    double xnew = x + h;
    if (rtt_dsxx::soft_equiv(xnew, x))
      throw std::range_error("stepsize underflow in rkqs");
  }
  if (errmax > ERRCON)
    hnext = SAFETY * h * std::pow(errmax, PGROW);
  else
    hnext = 5 * h;
  x += (hdid = h);
  for (unsigned i = 0; i < n; i++)
    y[i] = ytemp[i];
}

} // end namespace rtt_ode

#endif // ode_rkqs_i_hh

//------------------------------------------------------------------------------------------------//
// end of ode/rkqs.i.hh
//------------------------------------------------------------------------------------------------//
