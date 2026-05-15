//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/Draco_Mesh_Builder_pt.cc
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \brief  Draco_Mesh_Builder class header file.
 * \note   Copyright (C) 2018-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "mesh/Draco_Mesh_Builder.hh"
#include "mesh/Draco_Mesh_Builder.t.hh" // NOLINT(misc-include-cleaner)
#include "mesh/RTT_Draco_Mesh_Reader.hh"
#include "mesh/X3D_Draco_Mesh_Reader.hh"

namespace rtt_mesh {

template class Draco_Mesh_Builder<RTT_Draco_Mesh_Reader>;
template class Draco_Mesh_Builder<X3D_Draco_Mesh_Reader>;

} // end namespace rtt_mesh

//------------------------------------------------------------------------------------------------//
// end of mesh/Draco_Mesh_Builder_pt.cc
//------------------------------------------------------------------------------------------------//
