//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   min/mrqmin.hh
 * \author Kent Budge
 * \date   Fri Aug 7 11:11:31 MDT 2009
 * \brief  Levenberg-Marquardt method for nonlinear data fitting
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef min_mrqmin_hh
#define min_mrqmin_hh

#include "ds++/config.h"

namespace rtt_min {
//! Perform a nonlinear least-squares fit using Levenberg-Marquardt method
template <class RandomContainer, class RandomBoolContainer, class ModelFunction>
void mrqmin(RandomContainer const &x, RandomContainer const &y, RandomContainer const &sig,
            unsigned n, unsigned m, RandomContainer &a, RandomBoolContainer &ia,
            RandomContainer &covar, RandomContainer &alpha, unsigned p, double &chisq,
            ModelFunction funcs, double &alamda);

} // end namespace rtt_min

#endif // min_mrqmin_hh

//------------------------------------------------------------------------------------------------//
// end of min/mrqmin.hh
//------------------------------------------------------------------------------------------------//
