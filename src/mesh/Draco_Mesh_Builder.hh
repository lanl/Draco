//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/Draco_Mesh_Builder.hh
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \brief  Draco_Mesh_Builder class header file.
 * \note   Copyright (C) 2018-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_mesh_Draco_Mesh_Builder_hh
#define rtt_mesh_Draco_Mesh_Builder_hh

#include "dsxx/Assert.hh"
#include "mesh_element/Geometry.hh"
#include <memory>
#include <vector>

namespace rtt_mesh {

// Forward declare mesh class
class Draco_Mesh;

//================================================================================================//
/*!
 * \class Draco_Mesh_Builder
 *
 * \brief Draco_Mesh unstructured mesh builder from reader.
 *
 * \tparam FRT Reader is a template parameter ("FRT" = "Format Reader Type") to the builder class.
 * Hence all readers that instantiate the template must have a common set of accessors.  In
 * particular the reader must supply some functions in the RTT_Format_Reader class. The mesh
 * constructor arguments are stored to permit rebuilding the mesh without re-reading it, if the user
 * wants to make a modification to the vertex coordinates (without changing topology) or geometry,
 * for instance. */
//================================================================================================//

template <typename FRT> class Draco_Mesh_Builder {
private:
  // >>> DATA

  //! Pointer to reader
  const std::shared_ptr<FRT> reader;

  // >>> MESH CONSTRUCTOR ARGUMENTS

  unsigned dimension{0};
  std::vector<unsigned> num_faces_per_cell;
  std::vector<unsigned> cell_to_node_linkage;
  std::vector<unsigned> side_set_flag;
  std::vector<unsigned> side_node_count;
  std::vector<unsigned> side_to_node_linkage;
  std::vector<double> coordinates;
  std::vector<unsigned> global_node_number;
  std::vector<unsigned> num_nodes_per_face_per_cell;

public:
  //! Constructor
  // NOLINTNEXTLINE(performance-unnecessary-value-param)
  explicit Draco_Mesh_Builder(std::shared_ptr<FRT> reader_);

  // >>> SERVICES

  std::shared_ptr<Draco_Mesh> build_mesh(rtt_mesh_element::Geometry geometry);

  /*!
   * \brief Rebuild mesh after resetting some of the mesh data (see reset functions below)
   * \param[in] geometry enumeration of mesh geometry
   * \return shared pointer to the Draco_Mesh object */
  std::shared_ptr<Draco_Mesh> rebuild_mesh(rtt_mesh_element::Geometry geometry);

  // >>> ACCESSORS

  [[nodiscard]] std::vector<double> get_coordinates() const { return coordinates; }
  [[nodiscard]] std::vector<unsigned> get_side_set_flag() const { return side_set_flag; }

  // >>> MESH DATA RESET FUNCTIONS

  //! Reset coordinates (preserving the topology of the mesh)
  void reset_coordinates(const std::vector<double> &new_coordinates) {
    Require(new_coordinates.size() == coordinates.size());
    coordinates = new_coordinates;
  }

  //! Reset boundary conditions on physical cell faces
  void reset_side_set_flag(const std::vector<unsigned> &new_side_set_flag) {
    Require(new_side_set_flag.size() == side_set_flag.size());
    side_set_flag = new_side_set_flag;
  }
};

} // end namespace rtt_mesh

#endif // rtt_mesh_Draco_Mesh_Builder_hh

//------------------------------------------------------------------------------------------------//
// end of mesh/Draco_Mesh_Builder.hh
//------------------------------------------------------------------------------------------------//
