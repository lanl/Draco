//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/global_containers_pt.cc
 * \author Kent Budge
 * \date   Mon Mar 24 10:17:40 2008
 * \brief  Explicit template instatiations for class global_containers.
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/config.h"

#ifdef C4_MPI

#include "global_containers.t.hh"

namespace rtt_c4 {

template void global_merge(std::set<unsigned> &);
template void global_merge(std::map<unsigned, double> &);
template void global_merge(std::map<unsigned, int> &);
template void global_merge(std::map<unsigned, unsigned> &);
template void global_merge(std::map<unsigned, bool> &);

} // end namespace rtt_c4

#endif // C4_MPI

//------------------------------------------------------------------------------------------------//
// end of global_containers_pt.cc
//------------------------------------------------------------------------------------------------//
