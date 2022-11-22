//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/test/TestRTTMeshReader.cc
 * \author Thomas M. Evans
 * \date   Wed Mar 27 10:41:12 2002
 * \brief  RTT_Mesh_Reader test.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "RTT_Format_Reader/RTT_Mesh_Reader.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/path.hh"
#include <array>
#include <sstream>

//------------------------------------------------------------------------------------------------//
// Enum definitions and forward declarations
//------------------------------------------------------------------------------------------------//

using namespace std;
using namespace rtt_dsxx;
using rtt_mesh_element::Element_Definition;
using rtt_RTT_Format_Reader::RTT_Mesh_Reader;

enum Meshes { DEFINED, MESHES_LASTENTRY };

bool check_virtual(rtt_dsxx::UnitTest &ut, RTT_Mesh_Reader const &mesh, Meshes const &meshtype);

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void runTest(rtt_dsxx::UnitTest &ut) {
  // Find the mesh file
  string const inpPath = ut.getTestSourcePath();

  // New meshes added to this test will have to be added to the enumeration Meshes in the header
  // file.
  constexpr int MAX_MESHES = 1;
  std::array<string, MAX_MESHES> filename = {inpPath + string("rttdef.mesh")};
  Meshes mesh_type = MESHES_LASTENTRY;

  for (int mesh_number = 0; mesh_number < MAX_MESHES; mesh_number++) {
    // Construct an RTT_Mesh_Reader class object from the data in the specified mesh file.
    RTT_Mesh_Reader mesh(filename[mesh_number]);
    {
      ostringstream m;
      m << "Read " << filename[mesh_number] << " without coreing in or firing an assertion."
        << endl;
      PASSMSG(m.str());
    }
    bool all_passed = true;
    // The following switch allows addition of other meshes for testing, with the "DEFINED" mesh
    // providing an example.
    if (mesh_number == 0) {
      // Test all nested class accessor functions for a very simplistic mesh file (enum DEFINED).
      mesh_type = DEFINED;
      all_passed = all_passed && check_virtual(ut, mesh, mesh_type);
    } else {
      FAILMSG("Invalid mesh type encountered.");
      all_passed = false;
    }
    if (!all_passed) {
      ostringstream m;
      m << "Errors occured testing mesh number " << mesh_type << endl;
      FAILMSG(m.str());
    }
  }

  // Report results of test.
  if (ut.numFails == 0 && ut.numPasses > 0) {
    PASSMSG("All tests passed.");
  } else {
    FAILMSG("Some tests failed.");
  }
  return;
}

//------------------------------------------------------------------------------------------------//
bool check_virtual(rtt_dsxx::UnitTest &ut, RTT_Mesh_Reader const &mesh, Meshes const &meshtype) {
  // Save and reset at end of function
  bool passed(true);

  // Exercise the virtual accessor functions for this mesh.
  vector<vector<double>> node_coords;
  string node_coord_units;
  vector<vector<unsigned>> element_nodes;
  vector<Element_Definition::Element_Type> element_types;
  vector<Element_Definition::Element_Type> unique_element_types;
  map<string, set<unsigned>> node_sets;
  map<string, set<unsigned>> element_sets;
  string title;
  vector<double> coords(3, 0.0);
  vector<unsigned> side_nodes;
  set<unsigned> flag_nodes;
  set<unsigned> flag_elements;

  switch (meshtype) {
  case DEFINED:
    // set node coords per the input deck.
    node_coords.push_back(coords);
    coords[2] = 3.0;
    node_coords.push_back(coords);
    coords[1] = 2.0;
    coords[2] = 0.0;
    node_coords.push_back(coords);
    coords[0] = 1.0;
    coords[1] = 0.0;
    node_coords.push_back(coords);
    // set the coordinate units used for the nodes.
    node_coord_units = "cm";
    // load the node numbers for the single tet cell defined in the input file (note that the node
    // numbers are zero indexed).
    side_nodes.push_back(1);
    side_nodes.push_back(2);
    side_nodes.push_back(3);
    element_nodes.push_back(side_nodes);
    side_nodes.resize(0);
    side_nodes.push_back(0);
    side_nodes.push_back(3);
    side_nodes.push_back(2);
    element_nodes.push_back(side_nodes);
    side_nodes.resize(0);
    side_nodes.push_back(0);
    side_nodes.push_back(1);
    side_nodes.push_back(3);
    element_nodes.push_back(side_nodes);
    side_nodes.resize(0);
    side_nodes.push_back(0);
    side_nodes.push_back(2);
    side_nodes.push_back(1);
    element_nodes.push_back(side_nodes);
    side_nodes.resize(0);
    side_nodes.push_back(0);
    side_nodes.push_back(1);
    side_nodes.push_back(2);
    side_nodes.push_back(3);
    element_nodes.push_back(side_nodes);
    side_nodes.resize(0);
    // load the element types defined for RTT_Format according to the corresponding
    // Element_Definition::Element_Type.
    element_types.push_back(Element_Definition::TRI_3);
    element_types.push_back(Element_Definition::TRI_3);
    element_types.push_back(Element_Definition::TRI_3);
    element_types.push_back(Element_Definition::TRI_3);
    element_types.push_back(Element_Definition::TETRA_4);
    // load the unique element types defined for RTT_Format according to the corresponding
    // Element_Definition::Element_Type.
    unique_element_types.push_back(Element_Definition::NODE);
    unique_element_types.push_back(Element_Definition::BAR_2);
    unique_element_types.push_back(Element_Definition::TRI_3);
    unique_element_types.push_back(Element_Definition::QUAD_4);
    unique_element_types.push_back(Element_Definition::PYRA_5);
    unique_element_types.push_back(Element_Definition::TETRA_4);
    unique_element_types.push_back(Element_Definition::PENTA_6);
    unique_element_types.push_back(Element_Definition::HEXA_8);
    // load the node sets
    flag_nodes.insert(0);
    node_sets.insert(make_pair(string("node_type/interior"), flag_nodes));
    flag_nodes.erase(flag_nodes.begin(), flag_nodes.end());
    flag_nodes.insert(1);
    flag_nodes.insert(2);
    node_sets.insert(make_pair(string("node_type/dudded"), flag_nodes));
    flag_nodes.erase(flag_nodes.begin(), flag_nodes.end());
    flag_nodes.insert(3);
    node_sets.insert(make_pair(string("node_type/parent"), flag_nodes));
    flag_nodes.erase(flag_nodes.begin(), flag_nodes.end());
    flag_nodes.insert(0);
    flag_nodes.insert(1);
    node_sets.insert(make_pair(string("boundary/reflective"), flag_nodes));
    flag_nodes.erase(flag_nodes.begin(), flag_nodes.end());
    flag_nodes.insert(2);
    flag_nodes.insert(3);
    node_sets.insert(make_pair(string("boundary/vacuum"), flag_nodes));
    flag_nodes.erase(flag_nodes.begin(), flag_nodes.end());
    flag_nodes.insert(0);
    flag_nodes.insert(1);
    flag_nodes.insert(2);
    node_sets.insert(make_pair(string("source/no_source"), flag_nodes));
    flag_nodes.erase(flag_nodes.begin(), flag_nodes.end());
    flag_nodes.insert(3);
    node_sets.insert(make_pair(string("source/rad_source"), flag_nodes));
    flag_nodes.erase(flag_nodes.begin(), flag_nodes.end());
    // load the element (i.e., sides + cell) sets
    flag_elements.insert(1);
    flag_elements.insert(2);
    flag_elements.insert(3);
    element_sets.insert(make_pair(string("boundary/reflective"), flag_elements));
    flag_elements.erase(flag_elements.begin(), flag_elements.end());
    flag_elements.insert(0);
    element_sets.insert(make_pair(string("boundary/vacuum"), flag_elements));
    flag_elements.erase(flag_elements.begin(), flag_elements.end());
    flag_elements.insert(4);
    element_sets.insert(make_pair(string("material/control_rod"), flag_elements));
    flag_elements.erase(flag_elements.begin(), flag_elements.end());
    element_sets.insert(make_pair(string("material/shield"), flag_elements));
    flag_elements.erase(flag_elements.begin(), flag_elements.end());
    element_sets.insert(make_pair(string("rad_source/src_name1"), flag_elements));
    flag_elements.insert(4);
    element_sets.insert(make_pair(string("rad_source/src_name2"), flag_elements));
    flag_elements.erase(flag_elements.begin(), flag_elements.end());
    // set the mesh title
    title = "RTT_format mesh file definition, version 7.";
    break;

  default:
    FAILMSG("Invalid mesh type encountered.");
    return false;
  }
  // Check node coords
  if (node_coords != mesh.get_node_coords()) {
    FAILMSG("Node coordinates not obtained.");
  }
  // Check coordinate units.
  if (node_coord_units != mesh.get_node_coord_units()) {
    FAILMSG("Coordinate units not obtained.");
  }
  if (element_nodes != mesh.get_element_nodes()) {
    FAILMSG("Element nodes not obtained.");
  }
  // Check dimension
  if (mesh.get_dims_ndim() != 3) {
    FAILMSG("Expected dimension == 3.");
  }

  // Check Element Types.
  if (element_types != mesh.get_element_types()) {
    FAILMSG("Element Types not obtained.");
  }
  // Check Element Defs.
  if (mesh.get_element_defs().size() != 1) {
    FAILMSG("Element Defs not obtained.");
  }
  // Check Unique Element Types.
  if (unique_element_types != mesh.get_unique_element_types()) {
    FAILMSG("Unique Element Types not obtained.");
  }
  // Check node sets.
  if (node_sets != mesh.get_node_sets()) {
    FAILMSG("Node sets not obtained.");
  }
  // Check Element sets.
  if (element_sets != mesh.get_element_sets()) {
    FAILMSG("Element sets not obtained.");
  }
  // Check title.
  if (title != mesh.get_title()) {
    FAILMSG("Title not obtained.");
  }
  // Check invariant.
  if (!mesh.invariant()) {
    FAILMSG("Invariant not satisfied.");
  }

  if (passed) {
    PASSMSG("Got all virtual accessors.");
    return true;
  }

  FAILMSG("Errors in some virtual accessors.");
  return false;
}

//------------------------------------------------------------------------------------------------//
// Main
//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  try {
    ScalarUnitTest ut(argc, argv, release);
    runTest(ut);
  } catch (rtt_dsxx::assertion &err) {
    std::string msg = err.what();
    if (msg != std::string("Success")) {
      cout << "ERROR: While testing " << argv[0] << ", " << err.what() << endl;
      return 1;
    }
    return 0;
  } catch (exception &err) {
    cout << "ERROR: While testing " << argv[0] << ", " << err.what() << endl;
    return 1;
  } catch (...) {
    cout << "ERROR: While testing " << argv[0] << ", An unknown exception was thrown" << endl;
    return 1;
  }
  return 0;
}

//------------------------------------------------------------------------------------------------//
// end of TestRTTMeshReader.cc
//------------------------------------------------------------------------------------------------//
