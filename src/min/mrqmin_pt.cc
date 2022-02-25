//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   min/mrqmin_pt.cc
 * \author Kent Budge
 * \date   Fri Aug 7 11:11:31 MDT 2009
 * \brief  Specializations of mrqmin
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "mrqmin.t.hh"
#include <vector>

namespace rtt_min {

//------------------------------------------------------------------------------------------------//
// RandomContainer=vector<double>
//------------------------------------------------------------------------------------------------//

template void mrqmin(std::vector<double> const &x, std::vector<double> const &y,
                     std::vector<double> const &sig, unsigned n, unsigned m, std::vector<double> &a,
                     std::vector<bool> &ia, std::vector<double> &covar, std::vector<double> &alpha,
                     unsigned p, double &chisq,
                     void funcs(std::vector<double> const &, std::vector<double> const &, double &,
                                std::vector<double> &),
                     double &alamda);

} // end namespace rtt_min

//------------------------------------------------------------------------------------------------//
// end of mrqmin_pt.cc
//------------------------------------------------------------------------------------------------//
