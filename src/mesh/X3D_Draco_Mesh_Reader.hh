//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/X3D_Draco_Mesh_Reader.cc
 * \author Ryan Wollaeger <wollaeger@lanl.gov>, Kendra Long
 * \date   Wednesday, Jul 11, 2018, 14:24 pm
 * \brief  X3D_Draco_Mesh_Reader header file.
 * \note   Copyright (C) 2018-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_mesh_X3D_Draco_Mesh_Reader_hh
#define rtt_mesh_X3D_Draco_Mesh_Reader_hh

#include "Draco_Mesh_Reader.hh"
#include "ds++/Assert.hh"
#include "ds++/config.h"
#include <map>
#include <string>
#include <vector>

namespace rtt_mesh {

//================================================================================================//
/*!
 * \class X3D_Draco_Mesh_Reader
 *
 * \brief Read x3d file format and supply data to Draco_Mesh_Builder.
 *
 * This class parses mesh data from an x3d file format and manipulates it into a form that is
 * compatible to the Draco_Mesh_Builder class.  The parsing of the file stream to the map follows
 * work by Kendra Long on a parsing framework.
 *
 * For a description of the X3D mesh file layout, see:
 * https://xcp-confluence.lanl.gov/display/SIMC/Ingen->Flag+Data+Transfer
 *
 * \todo: Boundary condition data is evidently parsed separately in X3D, so meshes generated from
 *        this reader will not have side flag data (which ids boundary conditions.
 *
 * \todo: Consider using the Class_Parse_Table formalism developed by Kent Budge as an alternative.
 */
//================================================================================================//

class X3D_Draco_Mesh_Reader : public Draco_Mesh_Reader {
public:
  // >>> TYPEDEFS
  using Parsed_Element = std::pair<std::string, std::vector<std::string>>;
  using Parsed_Elements = std::vector<Parsed_Element>;

private:
  // >>> DATA

  //! File name
  const std::string filename;

  //! Boundary file names (optional data)
  const std::vector<std::string> bdy_filenames;

  //! Boundary conditions per boundary file (optional data)
  const std::vector<unsigned> bdy_flags;

  //! Vector of all parsed key-value data pairs (includes valueless delimiters)
  Parsed_Elements parsed_pairs;

  //! Header data map (header in x3d file)
  std::map<std::string, std::vector<size_t>> x3d_header_map;

  //! Node coordinates
  std::map<int, std::vector<double>> x3d_coord_map;

  //! Face-to-node map
  std::map<int, std::vector<int>> x3d_facenode_map;

  //! Cell-to-face map
  std::map<int, std::vector<int>> x3d_cellface_map;

  //! Matid map
  // As the x3d manual notes, matids are a little weird.
  std::vector<std::string> x3d_matids;

  //! Side-to-node map (0-based indices, unlike other maps)
  std::map<int, std::vector<unsigned>> x3d_sidenode_map;

  //! Side-to-flag map
  std::map<int, unsigned> x3d_sideflag_map;

  //! B.C. index-to-node map
  std::map<size_t, std::vector<unsigned>> bc_node_map;

public:
  //! Constructor
  explicit X3D_Draco_Mesh_Reader(const std::string &filename_,
                                 const std::vector<std::string> &bdy_filenames_ = {},
                                 const std::vector<unsigned> &bdy_flags_ = {});

  // >>> SERVICES

  void read_mesh() override;

  // >>> ACCESSORS

  // header data
  unsigned get_process() const {
    Check(x3d_header_map.at("process")[0] - 1 < UINT_MAX);
    return static_cast<unsigned>(x3d_header_map.at("process")[0] - 1);
  }
  unsigned get_numdim() const override {
    Check(x3d_header_map.at("numdim")[0] < UINT_MAX);
    return static_cast<unsigned>(x3d_header_map.at("numdim")[0]);
  }
  size_t get_numcells() const override { return x3d_header_map.at("elements")[0]; }
  size_t get_numnodes() const override { return x3d_header_map.at("nodes")[0]; }

  // coordinate data
  std::vector<double> get_nodecoord(size_t node) const override {
    Check(node + 1 < INT_MAX);
    return x3d_coord_map.at(static_cast<int>(node + 1));
  }

  // matid data
  std::string get_matid(size_t cell) const {
    Check(cell + 1 < INT_MAX);
    return x3d_matids.at(cell);
  }

  // accessors with deferred implementations
  unsigned get_celltype(size_t cell) const override;
  std::vector<unsigned> get_cellnodes(size_t cell) const override;
  std::vector<unsigned> get_cellfacenodes(size_t cell, size_t face) const;

  // data needed from x3d boundary file
  size_t get_numsides() const override { return x3d_sidenode_map.size(); }
  size_t get_sidetype(size_t side) const override {
    Check(side < INT_MAX);
    return x3d_sidenode_map.at(static_cast<int>(side)).size();
  }
  unsigned get_sideflag(size_t side) const override {
    Check(side < INT_MAX);
    return x3d_sideflag_map.at(static_cast<int>(side));
  }
  std::vector<unsigned> get_sidenodes(size_t side) const override {
    Check(side < INT_MAX);
    return x3d_sidenode_map.at(static_cast<int>(side));
  }
  const std::map<size_t, std::vector<unsigned>> &get_bc_node_map() const { return bc_node_map; }

private:
  // >>> SUPPORT FUNCTIONS

  Parsed_Elements::const_iterator find_iter_of_key(const Parsed_Elements &pairs, std::string key,
                                                   size_t start = 0);

  template <typename KT, typename VT>
  std::map<KT, std::vector<VT>> map_x3d_block(const std::string &block_name, size_t &dist);

  template <typename KT> KT convert_key(const std::string &skey);

  std::vector<unsigned> get_facenodes(size_t face) const;

  void read_bdy_files();
};

//------------------------------------------------------------------------------------------------//
// EXPLICIT SPECIALIZATIONS
//------------------------------------------------------------------------------------------------//
template <> std::string X3D_Draco_Mesh_Reader::convert_key<std::string>(const std::string &skey);

} // end namespace rtt_mesh

// implementation header file
#include "X3D_Draco_Mesh_Reader.i.hh"

#endif // rtt_mesh_X3D_Draco_Mesh_Reader_hh

//------------------------------------------------------------------------------------------------//
// end of mesh/X3D_Draco_Mesh_Reader.hh
//------------------------------------------------------------------------------------------------//
