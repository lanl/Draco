//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   min/mrqmin_pt.cc
 * \author Kent Budge
 * \date   Fri Aug 7 11:11:31 MDT 2009
 * \brief  Specializations of mrqmin
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "mrqmin.t.hh"
#include <vector>

namespace rtt_min {
using std::vector;

//------------------------------------------------------------------------------------------------//
// RandomContainer=vector<double>
//------------------------------------------------------------------------------------------------//

template void mrqmin(vector<double> const &x, vector<double> const &y, vector<double> const &sig,
                     unsigned n, unsigned m, vector<double> &a, vector<bool> &ia,
                     vector<double> &covar, vector<double> &alpha, unsigned p, double &chisq,
                     void funcs(vector<double> const &, vector<double> const &, double &,
                                vector<double> &),
                     double &alamda);

} // end namespace rtt_min

//------------------------------------------------------------------------------------------------//
// end of mrqmin_pt.cc
//------------------------------------------------------------------------------------------------//
