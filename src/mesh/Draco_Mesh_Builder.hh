//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/Draco_Mesh_Builder.hh
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \date   Friday, Jun 29, 2018, 09:55 am
 * \brief  Draco_Mesh_Builder class header file.
 * \note   Copyright (C) 2018-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_mesh_Draco_Mesh_Builder_hh
#define rtt_mesh_Draco_Mesh_Builder_hh

#include "mesh_element/Geometry.hh"
#include <memory>

namespace rtt_mesh {

// Forward declare mesh class
class Draco_Mesh;

//================================================================================================//
/*!
 * \class Draco_Mesh_Builder
 *
 * \brief Draco_Mesh unstructured mesh builder from reader.
 *
 * Reader is a template parameter ("FRT" = "Format Reader Type") to the builder class. Hence all
 * readers that instantiate the template must have a common set of accessors.  In particular the
 * reader must supply some functions in the RTT_Format_Reader class.
 */
//================================================================================================//

template <typename FRT> class Draco_Mesh_Builder {
private:
  // >>> DATA

  //! Pointer to reader
  std::shared_ptr<FRT> reader;

public:
  //! Constructor
  explicit Draco_Mesh_Builder(std::shared_ptr<FRT> reader_);

  // >>> SERVICES

  std::shared_ptr<Draco_Mesh> build_mesh(rtt_mesh_element::Geometry geometry);
};

} // end namespace rtt_mesh

#endif // rtt_mesh_Draco_Mesh_Builder_hh

//------------------------------------------------------------------------------------------------//
// end of mesh/Draco_Mesh_Builder.hh
//------------------------------------------------------------------------------------------------//
