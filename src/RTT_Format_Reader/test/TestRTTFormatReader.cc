//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/test/TestRTTFormatReader.cc
 * \author Thomas M. Evans
 * \date   Wed Mar 27 10:26:42 2002
 * \brief  RTT_Format_Reader test.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "TestRTTFormatReader.hh"
#include "RTT_Format_Reader/RTT_Mesh_Reader.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include <array>
#include <sstream>

using namespace std;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
map<Meshes, bool> Dims_validated;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void runTest(UnitTest &ut) {
  // Find the mesh file
  string const inpPath = ut.getTestSourcePath();

  // New meshes added to this test will have to be added to the enumeration Meshes in the header
  // file.
  int constexpr MAX_MESHES = 1;
  std::array<std::string, MAX_MESHES> filename = {inpPath + string("rttdef.mesh")};
  Meshes mesh_type = MESHES_LASTENTRY;

  for (int mesh_number = 0; mesh_number < MAX_MESHES; mesh_number++) {
    // Construct an RTT_Format_Reader class object from the data in the specified mesh file.
    RTT_Format_Reader mesh(filename[mesh_number]);
    {
      ostringstream m;
      m << "Read " << filename[mesh_number] << " without coreing in or firing an assertion."
        << std::endl;
      PASSMSG(m.str());
    }
    bool all_passed = true;
    // The following switch allows addition of other meshes for testing, with the "DEFINED" mesh
    // providing an example. Only the check_dims tests is required and it will be automatically
    // called by the other tests (with the exception of check header) if not invoked herein.  The
    // comparison data must also be provided for additional meshes within the switch structure
    // residing in the test functions.
    if (mesh_number == 0) {
      // Test all nested class accessor functions for a very simplistic mesh file (enum DEFINED).
      mesh_type = DEFINED;
      all_passed =
          check_header(mesh, mesh_type, ut) && check_dims(mesh, mesh_type, ut) &&
          check_node_flags(mesh, mesh_type, ut) && check_side_flags(mesh, mesh_type, ut) &&
          check_cell_flags(mesh, mesh_type, ut) && check_node_data_ids(mesh, mesh_type, ut) &&
          check_side_data_ids(mesh, mesh_type, ut) && check_cell_data_ids(mesh, mesh_type, ut) &&
          check_cell_defs(mesh, mesh_type, ut) && check_nodes(mesh, mesh_type, ut) &&
          check_sides(mesh, mesh_type, ut) && check_cells(mesh, mesh_type, ut) &&
          check_node_data(mesh, mesh_type, ut) && check_side_data(mesh, mesh_type, ut) &&
          check_cell_data(mesh, mesh_type, ut);
    } else {
      FAILMSG("Invalid mesh type encountered.");
      all_passed = false;
    }

    if (!all_passed) {
      ostringstream m;
      m << "Errors occured testing mesh number " << mesh_type << std::endl;
      FAILMSG(m.str());
    }
  }

  try {
    RTT_Format_Reader reader("no such");
    FAILMSG("did NOT detect nonexistent file");
  } catch (...) {
    PASSMSG("detected nonexistent file correctly");
  }

  // Report results of test.
  if (ut.numFails == 0 && ut.numPasses > 0)
    PASSMSG("All tests passed.");
  else
    FAILMSG("Some tests failed.");
  return;
}

//------------------------------------------------------------------------------------------------//
bool verify_Dims(RTT_Format_Reader const &mesh, Meshes const &meshtype, UnitTest &ut) {
  // Verify that the Dims data was previously validated.
  if (!Dims_validated.count(meshtype))
    check_dims(mesh, meshtype, ut);
  // Return the integrity state of the Dims data.
  return Dims_validated.find(meshtype)->second;
}

//------------------------------------------------------------------------------------------------//
bool check_header(RTT_Format_Reader const &mesh, Meshes const &meshtype, UnitTest &ut) {
  // Exercise the header accessor functions for this mesh.
  bool all_passed = true;
  std::string version, title, date;
  size_t cycle, ncomments;
  std::vector<std::string> comments;
  double time;

  switch (meshtype) {
  case DEFINED:
    version = "v1.0.0";
    title = "RTT_format mesh file definition, version 7.";
    date = "24 Jul 97";
    cycle = 1;
    time = 0.0;
    ncomments = 3;
    comments.emplace_back("One tet mesh in an RTT mesh file format.");
    comments.emplace_back("Date     : 24 Jul 97");
    comments.emplace_back("Author(s): H. Trease, J.McGhee");
    break;

  default:
    FAILMSG("Invalid mesh type encountered.");
    all_passed = false;
    return all_passed;
  }
  // Check version.
  if (mesh.get_header_version() != version) {
    FAILMSG("Header version not obtained.");
    all_passed = false;
  }
  // Check title.
  if (mesh.get_header_title() != title) {
    FAILMSG("Header title not obtained.");
    all_passed = false;
  }
  // Check date.
  if (mesh.get_header_date() != date) {
    FAILMSG("Header date not obtained.");
    all_passed = false;
  }
  // Check cycle.
  if (mesh.get_header_cycle() != cycle) {
    FAILMSG("Header cycle not obtained.");
    all_passed = false;
  }
  // Check time.
  if (!rtt_dsxx::soft_equiv(mesh.get_header_time(), time)) {
    FAILMSG("Header time not obtained.");
    all_passed = false;
  }
  // Check ncomments.
  if (mesh.get_header_ncomments() != ncomments) {
    FAILMSG("Header ncomments not obtained.");
    all_passed = false;
  }
  // Check comments.
  bool got_comments = true;
  for (size_t i = 0; i < ncomments; i++)
    if (comments[i] != mesh.get_header_comments(i))
      got_comments = false;
  if (!got_comments) {
    FAILMSG("Header comments not obtained.");
    all_passed = false;
  }
  // Check that all Header class accessors passed their tests.
  if (all_passed) {
    PASSMSG("Got all Header accessors.");
  } else {
    FAILMSG("Errors in some Header accessors.");
  }

  return all_passed;
}

//------------------------------------------------------------------------------------------------//
bool check_dims(RTT_Format_Reader const &mesh, Meshes const &meshtype, UnitTest &ut) {
  // Exercise the dims accessor functions for this mesh.
  bool all_passed = true;
  std::string coor_units, prob_time_units;
  size_t ncell_defs;
  size_t nnodes_max;
  size_t nsides_max;
  size_t nnodes_side_max;
  size_t ndim;
  size_t ndim_topo;
  size_t nnodes;
  size_t nnode_flag_types;
  std::vector<size_t> nnode_flags;
  size_t nnode_data;
  size_t nsides;
  size_t nside_types;
  std::vector<int> side_types;
  size_t nside_flag_types;
  std::vector<size_t> nside_flags;
  size_t nside_data;
  size_t ncells;
  size_t ncell_types;
  std::vector<int> cell_types;
  size_t ncell_flag_types;
  std::vector<size_t> ncell_flags;
  size_t ncell_data;

  switch (meshtype) {
  case DEFINED:
    coor_units = "cm";
    prob_time_units = "s";
    ncell_defs = 8;
    nnodes_max = 8;
    nsides_max = 6;
    nnodes_side_max = 4;
    ndim = 3;
    ndim_topo = 3;
    nnodes = 4;
    nnode_flag_types = 3;
    nnode_flags.push_back(3);
    nnode_flags.push_back(2);
    nnode_flags.push_back(2);
    nnode_data = 3;
    nsides = 4;
    nside_types = 1;
    // All side types are decremented relative to the value in the input file for zero indexing.
    side_types.push_back(2);
    nside_flag_types = 1;
    nside_flags.push_back(2);
    nside_data = 2;
    ncells = 1;
    ncell_types = 1;
    // All cell types are decremented relative to the value in the input file for zero indexing.
    cell_types.push_back(5);
    ncell_flag_types = 2;
    ncell_flags.push_back(2);
    ncell_flags.push_back(2);
    ncell_data = 1;
    break;

  default:
    FAILMSG("Invalid mesh type encountered.");
    all_passed = false;
    return all_passed;
  }
  // Check coordinate units.
  if (mesh.get_dims_coor_units() != coor_units) {
    FAILMSG("Dims coor_units not obtained.");
    all_passed = false;
  }
  // Check problem time units.
  if (mesh.get_dims_prob_time_units() != prob_time_units) {
    FAILMSG("Dims prob_time_units not obtained.");
    all_passed = false;
  }
  // Check number of cell definitions.
  if (mesh.get_dims_ncell_defs() != ncell_defs) {
    FAILMSG("Dims ncell_defs not obtained.");
    all_passed = false;
  }
  // Check maximum number of nodes for cells in the "cell_defs" block.
  if (mesh.get_dims_nnodes_max() != nnodes_max) {
    FAILMSG("Dims nnodes_max not obtained.");
    all_passed = false;
  }
  // Check maximum number of sides for cells in the "cell_defs" block.
  if (mesh.get_dims_nsides_max() != nsides_max) {
    FAILMSG("Dims nsides_max not obtained.");
    all_passed = false;
  }
  // Check maximum number of nodes/side for cells in the "cell_defs" block.
  if (mesh.get_dims_nnodes_side_max() != nnodes_side_max) {
    FAILMSG("Dims nnodes_side_max not obtained.");
    all_passed = false;
  }
  // Check number of spatial dimensions.
  if (mesh.get_dims_ndim() != ndim) {
    FAILMSG("Dims ndim not obtained.");
    all_passed = false;
  }
  // Check number of topological dimensions.
  if (mesh.get_dims_ndim_topo() != ndim_topo) {
    FAILMSG("Dims ndim_topo not obtained.");
    all_passed = false;
  }
  // Check total number of nodes in the mesh.
  if (mesh.get_dims_nnodes() != nnodes) {
    FAILMSG("Dims nnodes not obtained.");
    all_passed = false;
  }
  // Check number of node flag types.
  if (mesh.get_dims_nnode_flag_types() != nnode_flag_types) {
    FAILMSG("Dims nnode_flag_types not obtained.");
    all_passed = false;
  }
  // Check number of flags/node flag type.
  bool got_nnode_flags = true;
  for (size_t f = 0; f < nnode_flag_types; f++)
    if (mesh.get_dims_nnode_flags(f) != nnode_flags[f])
      got_nnode_flags = false;
  if (!got_nnode_flags) {
    FAILMSG("Dims nnode_flags not obtained.");
    all_passed = false;
  }
  // Check number of node data fields.
  if (mesh.get_dims_nnode_data() != nnode_data) {
    FAILMSG("Dims nnode_data not obtained.");
    all_passed = false;
  }
  // Check number of sides in the mesh.
  if (mesh.get_dims_nsides() != nsides) {
    FAILMSG("Dims nsides not obtained.");
    all_passed = false;
  }
  // Check number of side types actually present in "side" block.
  if (mesh.get_dims_nside_types() != nside_types) {
    FAILMSG("Dims nside_types not obtained.");
    all_passed = false;
  }
  // Check side type indexes used in "side" block.
  bool got_side_types = true;
  for (size_t s = 0; s < nside_types; s++)
    if (mesh.get_dims_side_types(s) != side_types[s])
      got_side_types = false;
  if (!got_side_types) {
    FAILMSG("Dims side_types not obtained.");
    all_passed = false;
  }
  // Check number of side flag types.
  if (mesh.get_dims_nside_flag_types() != nside_flag_types) {
    FAILMSG("Dims nside_flag_types not obtained.");
    all_passed = false;
  }
  // Check number of side flags/side flag type.
  bool got_nside_flags = true;
  for (size_t f = 0; f < nside_flag_types; f++)
    if (mesh.get_dims_nside_flags(f) != nside_flags[f])
      got_nside_flags = false;
  if (!got_nside_flags) {
    FAILMSG("Dims nside_flags not obtained.");
    all_passed = false;
  }
  // Check number of side data fields.
  if (mesh.get_dims_nside_data() != nside_data) {
    FAILMSG("Dims nside_data not obtained.");
    all_passed = false;
  }
  // Check total number of cells in the mesh.
  if (mesh.get_dims_ncells() != ncells) {
    FAILMSG("Dims ncells not obtained.");
    all_passed = false;
  }
  // Check number of cell types actually present in "cells" block.
  if (mesh.get_dims_ncell_types() != ncell_types) {
    FAILMSG("Dims ncell_types not obtained.");
    all_passed = false;
  }
  // Check cell type indexes used in "cells" block.
  bool got_ncell_types = true;
  for (size_t f = 0; f < ncell_types; f++)
    if (mesh.get_dims_cell_types(f) != cell_types[f])
      got_ncell_types = false;
  if (!got_ncell_types) {
    FAILMSG("Dims cell_types not obtained.");
    all_passed = false;
  }
  // Check number of cell flag types.
  if (mesh.get_dims_ncell_flag_types() != ncell_flag_types) {
    FAILMSG("Dims ncell_flag_types not obtained.");
    all_passed = false;
  }
  // Check number of flags/cell flag type.
  bool got_ncell_flags = true;
  for (size_t f = 0; f < ncell_flag_types; f++)
    if (mesh.get_dims_ncell_flags(f) != ncell_flags[f])
      got_ncell_flags = false;
  if (!got_ncell_flags) {
    FAILMSG("Dims ncell_flags not obtained.");
    all_passed = false;
  }
  // Check number of cell data fields.
  if (mesh.get_dims_ncell_data() != ncell_data) {
    FAILMSG("Dims ncell_data not obtained.");
    all_passed = false;
  }
  // Check that all Dims class accessors passed their tests.
  if (all_passed) {
    PASSMSG("Got all Dims accessors.");
  } else {
    FAILMSG("Errors in some Dims accessors.");
  }

  // Retain the result of testing the Dims integrity for this mesh type.
  Dims_validated.insert(std::make_pair(meshtype, all_passed));

  return all_passed;
}

//------------------------------------------------------------------------------------------------//
bool check_node_flags(RTT_Format_Reader const &mesh, Meshes const &meshtype, UnitTest &ut) {
  // Return if the Dims data is corrupt.
  if (!verify_Dims(mesh, meshtype, ut))
    return false;

  // Exercise the node_flags accessor functions for this mesh.
  bool all_passed = true;
  std::vector<std::string> flagTypes;
  std::vector<std::vector<std::pair<int, std::string>>> flag_num_name;
  std::vector<std::pair<int, std::string>> num_name;
  int ntype, bndry, src;

  switch (meshtype) {
  case DEFINED:
    flagTypes.emplace_back("node_type");
    num_name.emplace_back(std::make_pair(11, std::string("interior")));
    num_name.emplace_back(std::make_pair(21, std::string("dudded")));
    num_name.emplace_back(std::make_pair(6, std::string("parent")));
    flag_num_name.emplace_back(num_name);
    num_name.resize(0);
    ntype = 0;
    flagTypes.emplace_back("boundary");
    num_name.emplace_back(std::make_pair(1, std::string("reflective")));
    num_name.emplace_back(std::make_pair(4, std::string("vacuum")));
    flag_num_name.emplace_back(num_name);
    num_name.resize(0);
    bndry = 1;
    flagTypes.emplace_back("source");
    num_name.emplace_back(std::make_pair(101, std::string("no_source")));
    num_name.emplace_back(std::make_pair(22, std::string("rad_source")));
    flag_num_name.emplace_back(num_name);
    num_name.resize(0);
    src = 2;
    break;

  default:
    FAILMSG("Invalid mesh type encountered.");
    all_passed = false;
    return all_passed;
  }

  // Check node flag types.
  bool got_node_flag_types = true;
  for (size_t i = 0; i < mesh.get_dims_nnode_flag_types(); i++)
    if (flagTypes[i] != mesh.get_node_flags_flag_type(i))
      got_node_flag_types = false;
  if (!got_node_flag_types) {
    FAILMSG("Node Flags flag_types not obtained.");
    all_passed = false;
  }
  // Check node flag node_type flag number.
  if (ntype != mesh.get_node_flags_flag_type_index(flagTypes[ntype])) {
    FAILMSG("Node Flags node_type flag not obtained.");
    all_passed = false;
  }
  // Check node flag boundary flag number.
  if (bndry != mesh.get_node_flags_flag_type_index(flagTypes[bndry])) {
    FAILMSG("Node Flags boundary flag not obtained.");
    all_passed = false;
  }
  // Check node flag source flag number.
  if (src != mesh.get_node_flags_flag_type_index(flagTypes[src])) {
    FAILMSG("Node Flags source flag not obtained.");
    all_passed = false;
  }
  // Check node flag numbers for each of the flag types.
  bool got_node_flag_numbers = true;
  for (size_t i = 0; i < mesh.get_dims_nnode_flag_types(); i++) {
    num_name = flag_num_name[i];
    for (size_t j = 0; j < mesh.get_dims_nnode_flags(i); j++)
      if (num_name[j].first != mesh.get_node_flags_flag_number(i, j))
        got_node_flag_numbers = false;
  }
  if (!got_node_flag_numbers) {
    FAILMSG("Node Flags flag_numbers not obtained.");
    all_passed = false;
  }
  // Check number of flags for each node flag type.
  bool got_node_flag_size = true;
  for (size_t i = 0; i < mesh.get_dims_nnode_flag_types(); i++) {
    if (flag_num_name[i].size() != mesh.get_node_flags_flag_size(i))
      got_node_flag_size = false;
  }
  if (!got_node_flag_size) {
    FAILMSG("Node Flags flag_size not obtained.");
    all_passed = false;
  }
  // Check node flag names for each of the flag types.
  bool got_node_flag_name = true;
  for (size_t i = 0; i < mesh.get_dims_nnode_flag_types(); i++) {
    num_name = flag_num_name[i];
    for (size_t j = 0; j < mesh.get_dims_nnode_flags(i); j++)
      if (num_name[j].second != mesh.get_node_flags_flag_name(i, j))
        got_node_flag_name = false;
  }
  if (!got_node_flag_name) {
    FAILMSG("Node Flags flag_name not obtained.");
    all_passed = false;
  }

  if (all_passed) {
    PASSMSG("Got all NodeFlags accessors.");
  } else {
    FAILMSG("Errors in some NodeFlags accessors.");
  }

  return all_passed;
}

//------------------------------------------------------------------------------------------------//
bool check_side_flags(RTT_Format_Reader const &mesh, Meshes const &meshtype, UnitTest &ut) {
  // Return if the Dims data is corrupt.
  if (!verify_Dims(mesh, meshtype, ut))
    return false;

  // Exercise the side_flags accessor functions for this mesh.
  bool all_passed = true;
  std::vector<std::string> flagTypes;
  std::vector<std::vector<std::pair<int, std::string>>> flag_num_name;
  std::vector<std::pair<int, std::string>> num_name;
  int bndry;

  switch (meshtype) {
  case DEFINED:
    flagTypes.emplace_back("boundary");
    num_name.emplace_back(std::make_pair(1, std::string("reflective")));
    num_name.emplace_back(std::make_pair(2, std::string("vacuum")));
    flag_num_name.emplace_back(num_name);
    num_name.resize(0);
    bndry = 0;
    break;

  default:
    FAILMSG("Invalid mesh type encountered.");
    all_passed = false;
    return all_passed;
  }

  // Check side flag types.
  bool got_side_flag_types = true;
  for (size_t i = 0; i < mesh.get_dims_nside_flag_types(); i++)
    if (flagTypes[i] != mesh.get_side_flags_flag_type(i))
      got_side_flag_types = false;
  if (!got_side_flag_types) {
    FAILMSG("Side Flags flag_types not obtained.");
    all_passed = false;
  }
  // Check side flag boundary flag number.
  if (bndry != mesh.get_side_flags_flag_type_index(flagTypes[bndry])) {
    FAILMSG("Side Flags boundary flag not obtained.");
    all_passed = false;
  }
  // Check side flag numbers for each of the flag types.
  bool got_side_flag_numbers = true;
  for (size_t i = 0; i < mesh.get_dims_nside_flag_types(); i++) {
    num_name = flag_num_name[i];
    for (size_t j = 0; j < mesh.get_dims_nside_flags(i); j++)
      if (num_name[j].first != mesh.get_side_flags_flag_number(i, j))
        got_side_flag_numbers = false;
  }
  if (!got_side_flag_numbers) {
    FAILMSG("Side Flags flag_numbers not obtained.");
    all_passed = false;
  }
  // Check number of flags for each side flag type.
  bool got_side_flag_size = true;
  for (size_t i = 0; i < mesh.get_dims_nside_flag_types(); i++) {
    if (flag_num_name[i].size() != mesh.get_side_flags_flag_size(i))
      got_side_flag_size = false;
  }
  if (!got_side_flag_size) {
    FAILMSG("Side Flags flag_size not obtained.");
    all_passed = false;
  }
  // Check side flag names for each of the flag types.
  bool got_side_flag_name = true;
  for (size_t i = 0; i < mesh.get_dims_nside_flag_types(); i++) {
    num_name = flag_num_name[i];
    for (size_t j = 0; j < mesh.get_dims_nside_flags(i); j++)
      if (num_name[j].second != mesh.get_side_flags_flag_name(i, j))
        got_side_flag_name = false;
  }
  if (!got_side_flag_name) {
    FAILMSG("Side Flags flag_name not obtained.");
    all_passed = false;
  }

  if (all_passed) {
    PASSMSG("Got all SideFlags accessors.");
  } else {
    FAILMSG("Errors in some SideFlags accessors.");
  }

  return all_passed;
}

//------------------------------------------------------------------------------------------------//
bool check_cell_flags(RTT_Format_Reader const &mesh, Meshes const &meshtype, UnitTest &ut) {
  // Return if the Dims data is corrupt.
  if (!verify_Dims(mesh, meshtype, ut))
    return false;

  // Exercise the cell_flags accessor functions for this mesh.
  bool all_passed = true;
  std::vector<std::string> flagTypes;
  std::vector<std::vector<std::pair<int, std::string>>> flag_num_name;
  std::vector<std::pair<int, std::string>> num_name;
  int matl, rsrc;

  switch (meshtype) {
  case DEFINED:
    flagTypes.emplace_back("material");
    num_name.emplace_back(std::make_pair(1, std::string("control_rod")));
    num_name.emplace_back(std::make_pair(2, std::string("shield")));
    flag_num_name.emplace_back(num_name);
    num_name.resize(0);
    matl = 0;
    flagTypes.emplace_back("rad_source");
    num_name.emplace_back(std::make_pair(1, std::string("src_name1")));
    num_name.emplace_back(std::make_pair(2, std::string("src_name2")));
    flag_num_name.emplace_back(num_name);
    num_name.resize(0);
    rsrc = 1;
    break;

  default:
    FAILMSG("Invalid mesh type encountered.");
    all_passed = false;
    return all_passed;
  }

  // Check cell flag types.
  bool got_cell_flag_types = true;
  for (size_t i = 0; i < mesh.get_dims_ncell_flag_types(); i++)
    if (flagTypes[i] != mesh.get_cell_flags_flag_type(i))
      got_cell_flag_types = false;
  if (!got_cell_flag_types) {
    FAILMSG("Cell Flags flag_types not obtained.");
    all_passed = false;
  }
  // Check cell flag material flag number.
  if (matl != mesh.get_cell_flags_flag_type_index(flagTypes[matl])) {
    FAILMSG("Cell Flags material flag not obtained.");
    all_passed = false;
  }
  // Check cell flag radiation source flag number.
  if (rsrc != mesh.get_cell_flags_flag_type_index(flagTypes[rsrc])) {
    FAILMSG("Cell Flags volume source flag not obtained.");
    all_passed = false;
  }
  // Check cell flag numbers for each of the flag types.
  bool got_cell_flag_numbers = true;
  for (size_t i = 0; i < mesh.get_dims_ncell_flag_types(); i++) {
    num_name = flag_num_name[i];
    for (size_t j = 0; j < mesh.get_dims_ncell_flags(i); j++)
      if (num_name[j].first != mesh.get_cell_flags_flag_number(i, j))
        got_cell_flag_numbers = false;
  }
  if (!got_cell_flag_numbers) {
    FAILMSG("Cell Flags flag_numbers not obtained.");
    all_passed = false;
  }
  // Check number of flags for each cell flag type.
  bool got_cell_flag_size = true;
  for (size_t i = 0; i < mesh.get_dims_ncell_flag_types(); i++) {
    if (flag_num_name[i].size() != mesh.get_cell_flags_flag_size(i))
      got_cell_flag_size = false;
  }
  if (!got_cell_flag_size) {
    FAILMSG("Cell Flags flag_size not obtained.");
    all_passed = false;
  }
  // Check cell flag names for each of the flag types.
  bool got_cell_flag_name = true;
  for (size_t i = 0; i < mesh.get_dims_ncell_flag_types(); i++) {
    num_name = flag_num_name[i];
    for (size_t j = 0; j < mesh.get_dims_ncell_flags(i); j++)
      if (num_name[j].second != mesh.get_cell_flags_flag_name(i, j))
        got_cell_flag_name = false;
  }
  if (!got_cell_flag_name) {
    FAILMSG("Cell Flags flag_name not obtained.");
    all_passed = false;
  }

  if (all_passed) {
    PASSMSG("Got all CellFlags accessors.");
  } else {
    FAILMSG("Errors in some CellFlags accessors.");
  }

  return all_passed;
}

//------------------------------------------------------------------------------------------------//
bool check_node_data_ids(RTT_Format_Reader const &mesh, Meshes const &meshtype, UnitTest &ut) {
  // Return if the Dims data is corrupt.
  if (!verify_Dims(mesh, meshtype, ut))
    return false;

  // Exercise the node_data_ids accessor functions for this mesh.
  bool all_passed = true;
  std::vector<std::string> names;
  std::vector<std::string> units;

  switch (meshtype) {
  case DEFINED:
    names.emplace_back("density");
    units.emplace_back("gm/cm**3");
    names.emplace_back("ion_temp");
    units.emplace_back("keV");
    names.emplace_back("x_vel");
    units.emplace_back("cm/sec");
    break;

  default:
    FAILMSG("Invalid mesh type encountered.");
    all_passed = false;
    return all_passed;
  }

  // Check node data id names.
  bool got_node_data_id_names = true;
  for (size_t i = 0; i < mesh.get_dims_nnode_data(); i++) {
    if (names[i] != mesh.get_node_data_id_name(i))
      got_node_data_id_names = false;
  }
  if (!got_node_data_id_names) {
    FAILMSG("NodeDataIDs names not obtained.");
    all_passed = false;
  }

  // Check node data id units.
  bool got_node_data_id_units = true;
  for (size_t i = 0; i < mesh.get_dims_nnode_data(); i++) {
    if (units[i] != mesh.get_node_data_id_units(i))
      got_node_data_id_units = false;
  }
  if (!got_node_data_id_units) {
    FAILMSG("NodeDataIDs units not obtained.");
    all_passed = false;
  }

  if (all_passed) {
    PASSMSG("Got all NodeDataIDs accessors.");
  } else {
    FAILMSG("Errors in some NodeDataIDs accessors.");
  }

  return all_passed;
}

//------------------------------------------------------------------------------------------------//
bool check_side_data_ids(RTT_Format_Reader const &mesh, Meshes const &meshtype, UnitTest &ut) {
  // Return if the Dims data is corrupt.
  if (!verify_Dims(mesh, meshtype, ut))
    return false;

  // Exercise the side_data_ids accessor functions for this mesh.
  bool all_passed = true;
  std::vector<std::string> names;
  std::vector<std::string> units;

  switch (meshtype) {
  case DEFINED:
    names.emplace_back("density");
    units.emplace_back("gm/cm**3");
    names.emplace_back("ion_temp");
    units.emplace_back("keV");
    break;

  default:
    FAILMSG("Invalid mesh type encountered.");
    all_passed = false;
    return all_passed;
  }

  // Check side data id names.
  bool got_side_data_id_names = true;
  for (size_t i = 0; i < mesh.get_dims_nside_data(); i++) {
    if (names[i] != mesh.get_side_data_id_name(i))
      got_side_data_id_names = false;
  }
  if (!got_side_data_id_names) {
    FAILMSG("SideDataIDs names not obtained.");
    all_passed = false;
  }

  // Check side data id units.
  bool got_side_data_id_units = true;
  for (size_t i = 0; i < mesh.get_dims_nside_data(); i++) {
    if (units[i] != mesh.get_side_data_id_units(i))
      got_side_data_id_units = false;
  }
  if (!got_side_data_id_units) {
    FAILMSG("SideDataIDs units not obtained.");
    all_passed = false;
  }

  if (all_passed) {
    PASSMSG("Got all SideDataIDs accessors.");
  } else {
    FAILMSG("Errors in some SideDataIDs accessors.");
  }

  return all_passed;
}

//------------------------------------------------------------------------------------------------//
bool check_cell_data_ids(RTT_Format_Reader const &mesh, Meshes const &meshtype, UnitTest &ut) {
  // Return if the Dims data is corrupt.
  if (!verify_Dims(mesh, meshtype, ut))
    return false;

  // Exercise the cell_data_ids functions for this mesh.
  bool all_passed = true;
  std::vector<std::string> names;
  std::vector<std::string> units;

  switch (meshtype) {
  case DEFINED:
    names.emplace_back("density");
    units.emplace_back("gm/cm**3");
    break;

  default:
    FAILMSG("Invalid mesh type encountered.");
    all_passed = false;
    return all_passed;
  }

  // Check cell data id names.
  bool got_cell_data_id_names = true;
  for (size_t i = 0; i < mesh.get_dims_ncell_data(); i++) {
    if (names[i] != mesh.get_cell_data_id_name(i))
      got_cell_data_id_names = false;
  }
  if (!got_cell_data_id_names) {
    FAILMSG("CellDataIDs names not obtained.");
    all_passed = false;
  }

  // Check cell data id units.
  bool got_cell_data_id_units = true;
  for (size_t i = 0; i < mesh.get_dims_ncell_data(); i++) {
    if (units[i] != mesh.get_cell_data_id_units(i))
      got_cell_data_id_units = false;
  }
  if (!got_cell_data_id_units) {
    FAILMSG("CellDataIDs units not obtained.");
    all_passed = false;
  }

  if (all_passed) {
    PASSMSG("Got all CellDataIDs accessors.");
  } else {
    FAILMSG("Errors in some CellDataIDs accessors.");
  }

  return all_passed;
}

//------------------------------------------------------------------------------------------------//
//!\brief Test reading of the cell defs block from an RTT file
//------------------------------------------------------------------------------------------------//

bool check_cell_defs(RTT_Format_Reader const &mesh, Meshes const &meshtype, UnitTest &ut) {
  // Return if the Dims data is corrupt.
  if (!verify_Dims(mesh, meshtype, ut))
    return false;

  // Exercise the cell_defs accessor functions for this mesh.
  bool all_passed = true;
  std::vector<std::string> names;
  std::vector<size_t> nnodes;
  std::vector<size_t> nsides;
  std::vector<std::vector<int>> side_types;
  std::vector<std::vector<std::vector<unsigned>>> sides;
  std::vector<std::vector<std::vector<unsigned>>> ordered_sides;

  switch (meshtype) {
  case DEFINED:
    // Size std::vectors and load data consistent between both the standard and sorted cell
    //  definitions (names, nnodes, nsides).
    side_types.resize(8);
    sides.resize(8);
    ordered_sides.resize(8);
    // Size and load point data
    names.emplace_back("point");
    nnodes.emplace_back(1);
    nsides.emplace_back(0);
    side_types[0].resize(0);
    sides[0].resize(0);
    ordered_sides[0].resize(0);
    // Size and load line data.
    names.emplace_back("line");
    nnodes.emplace_back(2);
    nsides.emplace_back(2);
    side_types[1].resize(2, 0);
    sides[1].resize(2);
    ordered_sides[1].resize(2);
    ordered_sides[1][0].emplace_back(0);
    ordered_sides[1][1].emplace_back(1);
    // Size triangle data.
    names.emplace_back("triangle");
    nnodes.emplace_back(3);
    nsides.emplace_back(3);
    side_types[2].resize(3, 1);
    sides[2].resize(3);
    ordered_sides[2].resize(3);
    // Size quad data.
    names.emplace_back("quad");
    nnodes.emplace_back(4);
    nsides.emplace_back(4);
    side_types[3].resize(4, 1);
    sides[3].resize(4);
    ordered_sides[3].resize(4);
    // Size quad pyramid data.
    names.emplace_back("quad_pyr");
    nnodes.emplace_back(5);
    nsides.emplace_back(5);
    side_types[4].resize(5, 2);
    side_types[4][0] = 3;
    sides[4].resize(5);
    ordered_sides[4].resize(5);
    // Size tetrahedron data.
    names.emplace_back("tetrahedron");
    nnodes.emplace_back(4);
    nsides.emplace_back(4);
    side_types[5].resize(4, 2);
    sides[5].resize(4);
    ordered_sides[5].resize(4);
    // Size tri_prism data.
    names.emplace_back("tri_prism");
    nnodes.emplace_back(6);
    nsides.emplace_back(5);
    side_types[6].resize(5, 3);
    sides[6].resize(5);
    ordered_sides[6].resize(5);
    // Size hexahedron data.
    names.emplace_back("hexahedron");
    nnodes.emplace_back(8);
    nsides.emplace_back(6);
    side_types[7].resize(6, 3);
    sides[7].resize(6);
    ordered_sides[7].resize(6);
    // triangle
    ordered_sides[2][0].emplace_back(1);
    ordered_sides[2][0].emplace_back(2);
    ordered_sides[2][1].emplace_back(2);
    ordered_sides[2][1].emplace_back(0);
    ordered_sides[2][2].emplace_back(0);
    ordered_sides[2][2].emplace_back(1);
    // quad
    ordered_sides[3][0].emplace_back(0);
    ordered_sides[3][0].emplace_back(1);
    ordered_sides[3][1].emplace_back(1);
    ordered_sides[3][1].emplace_back(2);
    ordered_sides[3][2].emplace_back(2);
    ordered_sides[3][2].emplace_back(3);
    ordered_sides[3][3].emplace_back(3);
    ordered_sides[3][3].emplace_back(0);
    //quad_pyr
    ordered_sides[4][0].emplace_back(0);
    ordered_sides[4][0].emplace_back(3);
    ordered_sides[4][0].emplace_back(2);
    ordered_sides[4][0].emplace_back(1);
    ordered_sides[4][1].emplace_back(0);
    ordered_sides[4][1].emplace_back(1);
    ordered_sides[4][1].emplace_back(4);
    ordered_sides[4][2].emplace_back(1);
    ordered_sides[4][2].emplace_back(2);
    ordered_sides[4][2].emplace_back(4);
    ordered_sides[4][3].emplace_back(2);
    ordered_sides[4][3].emplace_back(3);
    ordered_sides[4][3].emplace_back(4);
    ordered_sides[4][4].emplace_back(3);
    ordered_sides[4][4].emplace_back(0);
    ordered_sides[4][4].emplace_back(4);
    //tetrahedron
    ordered_sides[5][0].emplace_back(1);
    ordered_sides[5][0].emplace_back(2);
    ordered_sides[5][0].emplace_back(3);
    ordered_sides[5][1].emplace_back(0);
    ordered_sides[5][1].emplace_back(3);
    ordered_sides[5][1].emplace_back(2);
    ordered_sides[5][2].emplace_back(0);
    ordered_sides[5][2].emplace_back(1);
    ordered_sides[5][2].emplace_back(3);
    ordered_sides[5][3].emplace_back(0);
    ordered_sides[5][3].emplace_back(2);
    ordered_sides[5][3].emplace_back(1);
    // tri_prism
    side_types[6][0] = side_types[6][1] = 2;
    ordered_sides[6][0].emplace_back(0);
    ordered_sides[6][0].emplace_back(2);
    ordered_sides[6][0].emplace_back(1);
    ordered_sides[6][1].emplace_back(3);
    ordered_sides[6][1].emplace_back(4);
    ordered_sides[6][1].emplace_back(5);
    ordered_sides[6][2].emplace_back(0);
    ordered_sides[6][2].emplace_back(1);
    ordered_sides[6][2].emplace_back(4);
    ordered_sides[6][2].emplace_back(3);
    ordered_sides[6][3].emplace_back(0);
    ordered_sides[6][3].emplace_back(3);
    ordered_sides[6][3].emplace_back(5);
    ordered_sides[6][3].emplace_back(2);
    ordered_sides[6][4].emplace_back(1);
    ordered_sides[6][4].emplace_back(2);
    ordered_sides[6][4].emplace_back(5);
    ordered_sides[6][4].emplace_back(4);
    // hexahedron
    ordered_sides[7][0].emplace_back(0);
    ordered_sides[7][0].emplace_back(3);
    ordered_sides[7][0].emplace_back(2);
    ordered_sides[7][0].emplace_back(1);
    ordered_sides[7][1].emplace_back(4);
    ordered_sides[7][1].emplace_back(5);
    ordered_sides[7][1].emplace_back(6);
    ordered_sides[7][1].emplace_back(7);
    ordered_sides[7][2].emplace_back(0);
    ordered_sides[7][2].emplace_back(1);
    ordered_sides[7][2].emplace_back(5);
    ordered_sides[7][2].emplace_back(4);
    ordered_sides[7][3].emplace_back(1);
    ordered_sides[7][3].emplace_back(2);
    ordered_sides[7][3].emplace_back(6);
    ordered_sides[7][3].emplace_back(5);
    ordered_sides[7][4].emplace_back(2);
    ordered_sides[7][4].emplace_back(3);
    ordered_sides[7][4].emplace_back(7);
    ordered_sides[7][4].emplace_back(6);
    ordered_sides[7][5].emplace_back(0);
    ordered_sides[7][5].emplace_back(4);
    ordered_sides[7][5].emplace_back(7);
    ordered_sides[7][5].emplace_back(3);
    // Load the (sorted) sides data from the ordered_sides data specified above. The first element
    // of ordered_sides and sides std::vectors corresponds to a point (no sides, and thus the second
    // "array" size is zero) so start at cell definition one.
    for (size_t c = 1; c < ordered_sides.size(); c++)
      for (size_t s = 0; s < ordered_sides[c].size(); s++)
        for (size_t n = 0; n < ordered_sides[c][s].size(); n++)
          sides[c][s].emplace_back(ordered_sides[c][s][n]);
    break;

  default:
    FAILMSG("Invalid mesh type encountered.");
    all_passed = false;
    return all_passed;
  }

  // Check cell definition names.
  bool got_cell_defs_names = true;
  for (size_t i = 0; i < mesh.get_dims_ncell_defs(); i++) {
    if (names[i] != mesh.get_cell_defs_name(i))
      got_cell_defs_names = false;
  }
  if (!got_cell_defs_names) {
    FAILMSG("CellDefs names not obtained.");
    all_passed = false;
  }

  // Check cell definition access
  {
    rtt_RTT_Format_Reader::CellDef const myCellDef(mesh.get_cell_defs_cell_def(0)); // NOLINT
    if (myCellDef.get_name() == std::string("point") && myCellDef.get_nnodes() == 1 &&
        myCellDef.get_nsides() == 0 && myCellDef.get_all_side_types().size() == 0 &&
        myCellDef.get_all_sides().size() == 0 && myCellDef.get_all_ordered_sides().size() == 0) {
      PASSMSG("mesh.get_cell_defs_cell_def() works.");
    } else {
      FAILMSG("mesh.get_cell_defs_cell_def() failed.");
    }
  }

  // Check get_cell_defs_node_map(int)
  {
    std::vector<unsigned> const myNodes = mesh.get_cell_defs_node_map(0); // NOLINT
    size_t mySize = myNodes.size();
    // std::cout << "mySize = " << mySize << std::endl;
    if (mySize == 0) {
      PASSMSG("get_cell_defs_node_map(int) returned an empty vector.");
    } else {
      PASSMSG("get_cell_defs_node_map(int) did not return an empty vector.");
      std::cout << "myNodes = { ";
      if (mySize > 0) {
        for (size_t i = 0; i < myNodes.size() - 1; ++i)
          std::cout << myNodes[i] << ", ";
        std::cout << myNodes[myNodes.size() - 1];
      }
      std::cout << " }." << std::endl;
    }
  }

  // Check get_cell_defs_redefined()
  {
    bool myBool = mesh.get_cell_defs_redefined();
    if (myBool)
      FAILMSG("Unexpected value for get_cell_defs_redefined(): Cells are "
              "redefined.");
    else
      PASSMSG("Expected value for get_cell_defs_redefined(): Cells are not "
              "redefined.");
  }

  // Check cell definition number of nodes.
  bool got_cell_defs_nnodes = true;
  for (size_t i = 0; i < mesh.get_dims_ncell_defs(); i++) {
    if (nnodes[i] != mesh.get_cell_defs_nnodes(i))
      got_cell_defs_nnodes = false;
  }
  if (!got_cell_defs_nnodes) {
    FAILMSG("CellDefs nnodes not obtained.");
    all_passed = false;
  }
  // Check cell definition number of sides.
  bool got_cell_defs_nsides = true;
  for (size_t i = 0; i < mesh.get_dims_ncell_defs(); i++) {
    if (nsides[i] != mesh.get_cell_defs_nsides(i))
      got_cell_defs_nsides = false;
  }
  if (!got_cell_defs_nsides) {
    FAILMSG("CellDefs nsides not obtained.");
    all_passed = false;
  }
  // Check cell definition side types.
  bool got_cell_defs_side_types = true;
  for (size_t i = 0; i < mesh.get_dims_ncell_defs(); i++) {
    for (size_t s = 0; s < mesh.get_cell_defs_nsides(i); s++)
      if (side_types[i][s] != mesh.get_cell_defs_side_types(i, s))
        got_cell_defs_side_types = false;
  }
  if (!got_cell_defs_side_types) {
    FAILMSG("CellDefs side_types not obtained.");
    all_passed = false;
  }

  // Check cell definition side sets.
  bool got_cell_defs_side = true;
  for (size_t i = 0; i < mesh.get_dims_ncell_defs(); i++) {
    for (size_t s = 0; s < mesh.get_cell_defs_nsides(i); s++)
      if (sides[i][s] != mesh.get_cell_defs_side(i, s))
        got_cell_defs_side = false;
  }
  if (!got_cell_defs_side) {
    FAILMSG("CellDefs side not obtained.");
    all_passed = false;
  }
  // Check cell definition ordered_side sets.
  bool got_cell_defs_ordered_side = true;
  for (size_t i = 0; i < mesh.get_dims_ncell_defs(); i++) {
    for (size_t s = 0; s < mesh.get_cell_defs_nsides(i); s++)
      if (ordered_sides[i][s] != mesh.get_cell_defs_ordered_side(i, s))
        got_cell_defs_ordered_side = false;
  }
  if (!got_cell_defs_ordered_side) {
    FAILMSG("CellDefs ordered_side not obtained.");
    all_passed = false;
  }

  if (all_passed) {
    PASSMSG("Got all CellDefs accessors.");
  } else {
    FAILMSG("Errors in some CellDefs accessors.");
  }

  return all_passed;
}

//------------------------------------------------------------------------------------------------//
// check_nodes
//------------------------------------------------------------------------------------------------//

bool check_nodes(RTT_Format_Reader const &mesh, Meshes const &meshtype, UnitTest &ut) {
  // Return if the Dims data is corrupt.
  if (!verify_Dims(mesh, meshtype, ut))
    return false;

  // Exercise the nodes accessor functions for this mesh.
  bool all_passed = true;
  std::vector<std::vector<double>> coords(mesh.get_dims_nnodes(),
                                          std::vector<double>(mesh.get_dims_ndim(), 0.0));
  std::vector<size_t> parents(mesh.get_dims_nnodes());
  std::vector<std::vector<int>> flags(mesh.get_dims_nnodes());

  switch (meshtype) {
  case DEFINED:
    // set node coords per the input deck.
    coords[2][1] = 2.0;
    // set node parents per the input deck.
    for (size_t i = 0; i < mesh.get_dims_nnodes(); i++)
      parents[i] = i;
    // load the node flags.
    flags[0].emplace_back(11);
    flags[0].emplace_back(1);
    flags[0].emplace_back(101);
    flags[2].emplace_back(21);
    flags[2].emplace_back(4);
    flags[2].emplace_back(101);
    coords[1][2] = 3.0;
    coords[3][0] = 1.0;
    flags[1].emplace_back(21);
    flags[1].emplace_back(1);
    flags[1].emplace_back(101);
    flags[3].emplace_back(6);
    flags[3].emplace_back(4);
    flags[3].emplace_back(22);
    break;

  default:
    FAILMSG("Invalid mesh type encountered.");
    all_passed = false;
    return all_passed;
  }

  // Check all of the node coords.
  if (coords != mesh.get_nodes_coords()) {
    FAILMSG("Nodes coordinates not obtained.");
    all_passed = false;
  }
  // Check all of the coordinate directions for a single node.
  bool got_node_coords = true;
  for (size_t i = 0; i < mesh.get_dims_nnodes(); i++)
    if (coords[i] != mesh.get_nodes_coords(i))
      got_node_coords = false;
  if (!got_node_coords) {
    FAILMSG("Node coordinates not obtained.");
    all_passed = false;
  }
  // Check a single coordinate direction for a single node.
  bool got_node_coord = true;
  for (size_t i = 0; i < mesh.get_dims_nnodes(); i++) {
    for (size_t d = 0; d < mesh.get_dims_ndim(); d++)
      if (!rtt_dsxx::soft_equiv(coords[i][d], mesh.get_nodes_coords(i, d)))
        got_node_coord = false;
  }
  if (!got_node_coord) {
    FAILMSG("Node coordinate not obtained.");
    all_passed = false;
  }
  // Check the node parents.
  bool got_nodes_parents = true;
  for (size_t i = 0; i < mesh.get_dims_nnodes(); i++) {
    Check(mesh.get_nodes_parents(i) >= 0);
    Check(parents[i] < INT_MAX);
    if (static_cast<int>(parents[i]) != mesh.get_nodes_parents(i))
      got_nodes_parents = false;
  }
  if (!got_nodes_parents) {
    FAILMSG("Nodes parents not obtained.");
    all_passed = false;
  }
  // Check the node flags.
  bool got_nodes_flags = true;
  for (size_t i = 0; i < mesh.get_dims_nnodes(); i++) {
    for (size_t f = 0; f < mesh.get_dims_nnode_flag_types(); f++)
      if (flags[i][f] != mesh.get_nodes_flags(i, f))
        got_nodes_flags = false;
  }
  if (!got_nodes_flags) {
    FAILMSG("Nodes flags not obtained.");
    all_passed = false;
  }

  if (all_passed) {
    PASSMSG("Got all Nodes accessors.");
  } else {
    FAILMSG("Errors in some Nodes accessors.");
  }

  return all_passed;
}

//------------------------------------------------------------------------------------------------//
bool check_sides(RTT_Format_Reader const &mesh, Meshes const &meshtype, UnitTest &ut) {
  // Return if the Dims data is corrupt.
  if (!verify_Dims(mesh, meshtype, ut))
    return false;

  // Exercise the sides accessor functions for this mesh.
  bool all_passed = true;
  std::vector<int> sideType;
  std::vector<std::vector<unsigned>> nodes;
  std::vector<std::vector<unsigned>> flags;

  switch (meshtype) {
  case DEFINED:
    // set side sideType per the input deck.
    sideType.resize(mesh.get_dims_nsides(), 2);
    // load the side nodes.
    nodes.resize(mesh.get_dims_nsides());
    nodes[0].emplace_back(1);
    nodes[0].emplace_back(2);
    nodes[0].emplace_back(3);
    nodes[1].emplace_back(0);
    nodes[1].emplace_back(3);
    nodes[1].emplace_back(2);
    nodes[2].emplace_back(0);
    nodes[2].emplace_back(1);
    nodes[2].emplace_back(3);
    nodes[3].emplace_back(0);
    nodes[3].emplace_back(2);
    nodes[3].emplace_back(1);
    // load the side flags.
    flags.resize(mesh.get_dims_nsides());
    flags[0].emplace_back(2);
    flags[1].emplace_back(1);
    flags[2].emplace_back(1);
    flags[3].emplace_back(1);

    break;

  default:
    FAILMSG("Invalid mesh type encountered.");
    all_passed = false;
    return all_passed;
  }

  // Check the side sideType.
  bool got_sides_type = true;
  for (size_t i = 0; i < mesh.get_dims_nsides(); i++)
    if (sideType[i] != mesh.get_sides_type(i))
      got_sides_type = false;
  if (!got_sides_type) {
    FAILMSG("Side type not obtained.");
    all_passed = false;
  }
  // Check all of the side nodes.
  if (nodes != mesh.get_sides_nodes()) {
    FAILMSG("Sides nodes not obtained.");
    all_passed = false;
  }
  // Check all of the nodes for a single side.
  bool got_side_nodes = true;
  for (size_t i = 0; i < mesh.get_dims_nsides(); i++)
    if (nodes[i] != mesh.get_sides_nodes(i))
      got_side_nodes = false;
  if (!got_side_nodes) {
    FAILMSG("Side nodes not obtained.");
    all_passed = false;
  }
  // Check a single node for a single side.
  bool got_side_node = true;
  for (size_t i = 0; i < mesh.get_dims_nsides(); i++) {
    for (size_t n = 0; n < mesh.get_cell_defs_nnodes(mesh.get_sides_type(i)); n++)
      if (nodes[i][n] != mesh.get_sides_nodes(i, n))
        got_side_node = false;
  }
  if (!got_side_node) {
    FAILMSG("Side node not obtained.");
    all_passed = false;
  }
  // Check the side flags.
  bool got_sides_flags = true;
  for (size_t i = 0; i < mesh.get_dims_nsides(); i++) {
    for (size_t f = 0; f < mesh.get_dims_nside_flag_types(); f++)
      if (flags[i][f] != mesh.get_sides_flags(i, f))
        got_sides_flags = false;
  }
  if (!got_sides_flags) {
    FAILMSG("Side flags not obtained.");
    all_passed = false;
  }

  if (all_passed) {
    PASSMSG("Got all Sides accessors.");
  } else {
    FAILMSG("Errors in some Sides accessors.");
  }

  return all_passed;
}

//------------------------------------------------------------------------------------------------//
bool check_cells(RTT_Format_Reader const &mesh, Meshes const &meshtype, UnitTest &ut) {
  // Return if the Dims data is corrupt.
  if (!verify_Dims(mesh, meshtype, ut))
    return false;

  // Exercise the cells functions for this mesh.
  bool all_passed = true;
  std::vector<int> cellType;
  std::vector<std::vector<unsigned>> nodes;
  std::vector<std::vector<unsigned>> flags;

  switch (meshtype) {
  case DEFINED:
    // set cell cellType per the input deck.
    cellType.resize(mesh.get_dims_ncells(), 5);
    // load the cell flags.
    flags.resize(mesh.get_dims_ncells());
    flags[0].emplace_back(1);
    flags[0].emplace_back(2);
    // Load the cell nodes.
    nodes.resize(mesh.get_dims_ncells());
    nodes[0].emplace_back(0);
    nodes[0].emplace_back(1);
    nodes[0].emplace_back(2);
    nodes[0].emplace_back(3);
    break;

  default:
    FAILMSG("Invalid mesh type encountered.");
    all_passed = false;
    return all_passed;
  }

  // Check the cell cellType.
  bool got_cells_type = true;
  for (size_t i = 0; i < mesh.get_dims_ncells(); i++)
    if (cellType[i] != mesh.get_cells_type(i))
      got_cells_type = false;
  if (!got_cells_type) {
    FAILMSG("Cell type not obtained.");
    all_passed = false;
  }
  // Check all of the cell nodes.
  if (nodes != mesh.get_cells_nodes()) {
    FAILMSG("Cells nodes not obtained.");
    all_passed = false;
  }
  // Check all of the nodes for a single cell.
  bool got_cell_nodes = true;
  for (size_t i = 0; i < mesh.get_dims_ncells(); i++)
    if (nodes[i] != mesh.get_cells_nodes(i))
      got_cell_nodes = false;
  if (!got_cell_nodes) {
    FAILMSG("Cell nodes not obtained.");
    all_passed = false;
  }
  // Check a single node for a single cell.
  bool got_cell_node = true;
  for (size_t i = 0; i < mesh.get_dims_ncells(); i++) {
    for (size_t n = 0; n < mesh.get_cell_defs_nnodes(mesh.get_cells_type(i)); n++)
      if (nodes[i][n] != mesh.get_cells_nodes(i, n))
        got_cell_node = false;
  }
  if (!got_cell_node) {
    FAILMSG("Cell node not obtained.");
    all_passed = false;
  }
  // Check the cell flags.
  bool got_cells_flags = true;
  for (size_t i = 0; i < mesh.get_dims_ncells(); i++) {
    for (size_t f = 0; f < mesh.get_dims_ncell_flag_types(); f++) {
      Check(mesh.get_cells_flags(i, f) >= 0);
      Check(flags[i][f] < INT_MAX);
      if (static_cast<int>(flags[i][f]) != mesh.get_cells_flags(i, f))
        got_cells_flags = false;
    }
  }
  if (!got_cells_flags) {
    FAILMSG("Cell flags not obtained.");
    all_passed = false;
  }

  if (all_passed) {
    PASSMSG("Got all Cells accessors.");
  } else {
    FAILMSG("Errors in some Cells accessors.");
  }

  return all_passed;
}

//------------------------------------------------------------------------------------------------//
bool check_node_data(RTT_Format_Reader const &mesh, Meshes const &meshtype, UnitTest &ut) {
  // Return if the Dims data is corrupt.
  if (!verify_Dims(mesh, meshtype, ut))
    return false;

  // Exercise the node_data functions for this mesh.
  bool all_passed = true;
  std::vector<std::vector<double>> data(mesh.get_dims_nnodes(),
                                        std::vector<double>(mesh.get_dims_nnode_data(), 0.0));

  switch (meshtype) {
  case DEFINED:
    // set node data per the input deck.
    data[1][2] = 3.0;
    data[2][1] = 2.0;
    data[3][0] = 1.0;
    break;

  default:
    FAILMSG("Invalid mesh type encountered.");
    all_passed = false;
    return all_passed;
  }

  // Check all of the node data.
  if (data != mesh.get_node_data()) {
    FAILMSG("NodeData not obtained for all nodes/fields.");
    all_passed = false;
  }
  // Check all of the data fields for a single node.
  bool got_node_data_fields = true;
  for (size_t i = 0; i < mesh.get_dims_nnodes(); i++)
    if (data[i] != mesh.get_node_data(i))
      got_node_data_fields = false;
  if (!got_node_data_fields) {
    FAILMSG("NodeData fields not obtained for a node.");
    all_passed = false;
  }
  // Check a single data field for a single node.
  bool got_node_data = true;
  for (size_t i = 0; i < mesh.get_dims_nnodes(); i++) {
    for (size_t d = 0; d < mesh.get_dims_nnode_data(); d++)
      if (!rtt_dsxx::soft_equiv(data[i][d], mesh.get_node_data(i, d)))
        got_node_data = false;
  }
  if (!got_node_data) {
    FAILMSG("NodeData value not obtained.");
    all_passed = false;
  }

  if (all_passed) {
    PASSMSG("Got all NodeData accessors.");
  } else {
    FAILMSG("Errors in some NodeData accessors.");
  }

  return all_passed;
}

//------------------------------------------------------------------------------------------------//
bool check_side_data(RTT_Format_Reader const &mesh, Meshes const &meshtype, UnitTest &ut) {
  // Return if the Dims data is corrupt.
  if (!verify_Dims(mesh, meshtype, ut))
    return false;

  // Exercise the side_data functions for this mesh.
  bool all_passed = true;
  std::vector<std::vector<double>> data(mesh.get_dims_nsides(),
                                        std::vector<double>(mesh.get_dims_nside_data(), 0.0));

  switch (meshtype) {
  case DEFINED:
    // set side data per the input deck.
    data[2][1] = 2.0;
    data[3][0] = 1.0;
    break;

  default:
    FAILMSG("Invalid mesh type encountered.");
    all_passed = false;
    return all_passed;
  }

  // Check all of the side data.
  if (data != mesh.get_side_data()) {
    FAILMSG("SideData not obtained for all sides/fields.");
    all_passed = false;
  }
  // Check all of the data fields for a single side.
  bool got_side_data_fields = true;
  for (size_t i = 0; i < mesh.get_dims_nsides(); i++)
    if (data[i] != mesh.get_side_data(i))
      got_side_data_fields = false;
  if (!got_side_data_fields) {
    FAILMSG("SideData fields not obtained for a side.");
    all_passed = false;
  }
  // Check a single data field for a single side.
  bool got_side_data = true;
  for (size_t i = 0; i < mesh.get_dims_nsides(); i++) {
    for (size_t d = 0; d < mesh.get_dims_nside_data(); d++)
      if (!rtt_dsxx::soft_equiv(data[i][d], mesh.get_side_data(i, d)))
        got_side_data = false;
  }
  if (!got_side_data) {
    FAILMSG("SideData value not obtained.");
    all_passed = false;
  }

  if (all_passed) {
    PASSMSG("Got all SideData accessors.");
  } else {
    FAILMSG("Errors in some SideData accessors.");
  }

  return all_passed;
}

//------------------------------------------------------------------------------------------------//
bool check_cell_data(RTT_Format_Reader const &mesh, Meshes const &meshtype, UnitTest &ut) {
  // Return if the Dims data is corrupt.
  if (!verify_Dims(mesh, meshtype, ut))
    return false;

  // Exercise the cell_data functions for this mesh.
  bool all_passed = true;
  std::vector<std::vector<double>> data;
  std::vector<double> cell_data(mesh.get_dims_ncell_data(), 0.0);

  switch (meshtype) {
  case DEFINED:
    // set cell data per the input deck.
    data.emplace_back(cell_data);
    break;

  default:
    FAILMSG("Invalid mesh type encountered.");
    all_passed = false;
    return all_passed;
  }

  // Check all of the cell data.
  if (data != mesh.get_cell_data()) {
    FAILMSG("CellData not obtained for all cells/fields.");
    all_passed = false;
  }
  // Check all of the data fields for a single cell.
  bool got_cell_data_fields = true;
  for (size_t i = 0; i < mesh.get_dims_ncells(); i++)
    if (data[i] != mesh.get_cell_data(i))
      got_cell_data_fields = false;
  if (!got_cell_data_fields) {
    FAILMSG("CellData fields not obtained for a cell.");
    all_passed = false;
  }
  // Check a single data field for a single cell.
  bool got_cell_data = true;
  for (size_t i = 0; i < mesh.get_dims_ncells(); i++) {
    for (size_t d = 0; d < mesh.get_dims_ncell_data(); d++)
      if (!rtt_dsxx::soft_equiv(data[i][d], mesh.get_cell_data(i, d)))
        got_cell_data = false;
  }
  if (!got_cell_data) {
    FAILMSG("CellData value not obtained.");
    all_passed = false;
  }

  if (all_passed) {
    PASSMSG(std::string("Got all CellData accessors."));
  } else {
    FAILMSG("Errors in some CellData accessors.");
  }

  return all_passed;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    runTest(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of TestRTTFormatReader.cc
//------------------------------------------------------------------------------------------------//
