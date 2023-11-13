//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/CAMR_Mesh_pt.cc
 * \author H. Park
 * \date   May 16 2023
 * \brief  Explicit template instantiations for class CAMR_Mesh.
 * \note   Copyright (C) 2023 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "CAMR_Mesh.t.hh"

namespace rtt_mesh {

template class CAMR_Mesh<2>;
template class CAMR_Mesh<3>;

} // end namespace rtt_mesh

//------------------------------------------------------------------------------------------------//
// end of CAMR_Mesh_pt.cc
//------------------------------------------------------------------------------------------------//
