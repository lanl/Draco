//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   min/mnbrak.hh
 * \author Kent Budge
 * \date   Tue Aug 17 15:30:23 2004
 * \brief  Bracket minimum of a function.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//

#ifndef min_mnbrak_hh
#define min_mnbrak_hh

#include "ds++/DracoMath.hh"
#include "ds++/dbc.hh"

namespace rtt_min {

//------------------------------------------------------------------------------------------------//
/*!
 * Bracket minimum of a function.
 *
 * \arg \a Function A function type supporting <code>double
 * operator()(double)</code>.
 *
 * \param[in,out] ax Lower bound of bracket interval.
 * \param[in,out] bx Upper bound of bracket interval.
 * \param[in,out] cx Estimate of minimum location within bracket.This will
 * \em not be a polished minimum; that is the job of rtt_min::brent.
 * \param[out] fa Function value at ax.
 * \param[out] fb Function value at bx.
 * \param[out] fc Function value at cx.
 * \param[in] func Function to be minimized.
 */

template <class Function>
void mnbrak(double &ax, double &bx, double &cx, double &fa, double &fb, double &fc, Function func) {
  using namespace std;
  using namespace rtt_dsxx;

  double const GOLD = 1.618034;
  double const GLIMIT = 100.0;
  double const TINY = 1.0e-20;

  fa = func(ax);
  fb = func(bx);
  if (fb > fa) {
    swap(ax, bx);
    swap(fb, fa);
  }
  cx = bx + GOLD * (bx - ax);
  fc = func(cx);
  while (fb > fc) {
    double const r = (bx - ax) * (fb - fc);
    double const q = (bx - cx) * (fb - fa);
    double u = bx - ((bx - cx) * q - (bx - ax) * r) / (2 * sign(max(fabs(q - r), TINY), q - r));
    double const ulim = bx + GLIMIT * (cx - bx);
    double fu;
    if ((bx - u) * (u - cx) > 0.0) {
      fu = func(u);
      if (fu < fc) {
        ax = bx;
        bx = u;
        fa = fb;
        fb = fu;
        return;
      } else if (fu > fb) {
        cx = u;
        fc = fu;
        return;
      }
      u = cx + GOLD * (cx - bx);
      fu = func(u);
    } else if ((cx - u) * (u - ulim) > 0) {
      u = ulim;
      fu = func(u);
    } else {
      u = cx + GOLD * (cx - bx);
      fu = func(u);
    }
    ax = bx;
    bx = cx;
    cx = u;
    fa = fb;
    fb = fc;
    fc = fu;
  }
}

} // end namespace rtt_min

#endif // min_mnbrak_hh

//------------------------------------------------------------------------------------------------//
// end of min/mnbrak.hh
//------------------------------------------------------------------------------------------------//
