//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   meshReaders/test/TestHexMeshReader.cc
 * \author Thomas M. Evans
 * \date   Tue Mar 26 16:37:01 2002
 * \brief  Hex Mesh Reader test.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "TestHexMeshReader.hh"
#include "ds++/Release.hh"
#include "meshReaders/Hex_Mesh_Reader.hh"
#include <cmath>
#include <sstream>

using namespace std;
using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void runTest(UnitTest &ut) {
  using rtt_mesh_element::Element_Definition;
  using rtt_meshReaders::Hex_Mesh_Reader;

  cout << "\n******* CIC-19 Hex Mesh Reader Tests *******" << std::endl;

  // Read and test a 1D mesh.g
  string const inpPath = ut.getTestSourcePath();

  std::string filename = inpPath + "slab.mesh";
  std::cout << "Creating mesh from file: " << filename << std::endl;
  Hex_Mesh_Reader mesh_1D(filename);
  {
    ostringstream msg;
    msg << "Read a 1D mesh without coreing in or firing an assertion.\n";
    PASSMSG(msg.str());
  }

  rtt_meshReaders_test::check_mesh(ut, mesh_1D, "slab");
  vector<std::shared_ptr<Element_Definition>> element_defs = mesh_1D.get_element_defs();
  if (element_defs.size() > 0)
    FAILMSG("element defs is NOT empty for slab");

  // Read and test a 2D mesh.
  filename = inpPath + "quad.mesh";
  std::cout << "Creating mesh from file: " << filename << std::endl;
  Hex_Mesh_Reader mesh_2D(filename);
  {
    ostringstream message;
    message << "Read a 2D mesh without coreing in or firing an assertion." << std::endl;
    PASSMSG(message.str());
  }
  rtt_meshReaders_test::check_mesh(ut, mesh_2D, "quad");
  element_defs = mesh_2D.get_element_defs();
  if (element_defs.size() > 0)
    FAILMSG("element defs is NOT empty for quad");

  // Read and test a 3D mesh.
  filename = inpPath + "cube.mesh";
  std::cout << "Creating mesh from file: " << filename << std::endl;
  Hex_Mesh_Reader mesh_3D(filename);
  {
    ostringstream message;
    message << "Read a 3D mesh without coreing in or firing an assertion." << std::endl;
    PASSMSG(message.str());
  }
  rtt_meshReaders_test::check_mesh(ut, mesh_3D, "cube");
  element_defs = mesh_3D.get_element_defs();
  if (element_defs.size() > 0)
    FAILMSG("element defs is NOT empty for cube");
  cout << endl;

  // Report results of test.
  if (ut.numFails == 0 && ut.numPasses > 0)
    ut.passes("All tests passed.");
  else
    ut.failure("Some tests failed.");
  return;
}

//------------------------------------------------------------------------------------------------//
namespace rtt_meshReaders_test {

bool check_mesh(UnitTest &ut, const rtt_meshReaders::Hex_Mesh_Reader &mesh,
                const std::string &testid) {

  // Exercize the accessor functions for this mesh and spot check
  // the results.

  bool pass_nc = check_nodes(ut, mesh, testid);
  bool pass_cu = check_node_units(ut, mesh);
  bool pass_ns = check_node_sets(ut, mesh, testid);
  bool pass_ti = check_title(ut, mesh);
  bool pass_gdn = check_get_dims_ndim(ut, mesh, testid);
  bool pass_en = check_element_nodes(ut, mesh, testid);
  bool pass_in = check_invariant(ut, mesh);
  bool pass_es = check_element_sets(ut, mesh, testid);
  bool pass_et = check_element_types(ut, mesh, testid);
  bool pass_ut = check_unique_element_types(ut, mesh, testid);

  bool all_passed(pass_nc && pass_cu && pass_ns && pass_ti && pass_gdn && pass_en && pass_in &&
                  pass_es && pass_et && pass_ut);

  if (!all_passed)
    FAILMSG("check_mesh failed");

  return all_passed;
}

bool check_nodes(UnitTest &ut, const rtt_meshReaders::Hex_Mesh_Reader &mesh,
                 const std::string &testid) {
  // Check node coords -- Need to do a "fuzzy" check here, these are doubles!
  bool pass_nc = true;
  std::vector<std::vector<double>> point_c = mesh.get_node_coords();
  if (testid == "slab")
    pass_nc = compare_double(point_c[0][0], 0.) && compare_double(point_c[10][0], 2.5) &&
              compare_double(point_c[100][0], 25.);
  else if (testid == "quad")
    pass_nc = compare_double(point_c[0][0], 0.) && compare_double(point_c[10][0], 12.5) &&
              compare_double(point_c[440][0], 25.) && compare_double(point_c[0][1], 0.) &&
              compare_double(point_c[10][1], 0.) && compare_double(point_c[440][1], 25.);
  else if (testid == "cube") {
    bool pass1 = compare_double(point_c[0][0], 0.) && compare_double(point_c[0][1], 0.) &&
                 compare_double(point_c[0][2], 0.);
    bool pass2 = compare_double(point_c[10][0], 0.8) && compare_double(point_c[10][1], 0.2) &&
                 compare_double(point_c[10][2], 0.);
    bool pass3 = compare_double(point_c[215][0], 1.) && compare_double(point_c[215][1], 1.) &&
                 compare_double(point_c[215][2], 1.);
    pass_nc = pass1 && pass2 && pass3;
  } else
    Insist(false, "Unrecognized test id string!");
  if (pass_nc) {
    ostringstream message;
    message << "Got node coordinates." << std::endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in node coordinates." << std::endl;
    FAILMSG(message.str());
  }
  return pass_nc;
}

bool check_node_units(UnitTest &ut, const rtt_meshReaders::Hex_Mesh_Reader &mesh) {
  // Check coordinate units.
  std::string punits = mesh.get_node_coord_units();
  bool pass_cu = punits == "unknown";
  if (pass_cu) {
    ostringstream message;
    message << "Got coordinate units." << std::endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in coordinate units." << std::endl;
    FAILMSG(message.str());
  }
  return pass_cu;
}

bool check_node_sets(UnitTest &ut, const rtt_meshReaders::Hex_Mesh_Reader &mesh,
                     const std::string &testid) {
  // Check node sets.
  std::map<std::string, std::set<unsigned>> ndsets = mesh.get_node_sets();
  bool pass_ns = ndsets.size() == 1;
  if (testid == "slab") {
    pass_ns = pass_ns && check_map(ndsets, "Interior", 0, 101);
  } else if (testid == "quad") {
    pass_ns = pass_ns && check_map(ndsets, "Interior", 0, 441);
  } else if (testid == "cube") {
    pass_ns = pass_ns && check_map(ndsets, "Interior", 0, 216);
  } else
    Insist(false, "Unrecognized test id string!");
  if (pass_ns) {
    ostringstream message;
    message << "Got node sets." << std::endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in node sets." << std::endl;
    FAILMSG(message.str());
  }
  return pass_ns;
}

bool check_title(UnitTest &ut, const rtt_meshReaders::Hex_Mesh_Reader &mesh) {
  // Check title.
  std::string title = mesh.get_title();
  bool pass_ti = title == "Untitled -- CIC-19 Hex Mesh";
  if (pass_ti) {
    ostringstream message;
    message << "Got title." << std::endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in title." << std::endl;
    FAILMSG(message.str());
  }
  return pass_ti;
}

bool check_get_dims_ndim(UnitTest &ut, rtt_meshReaders::Hex_Mesh_Reader const &mesh,
                         std::string const &testid) {
  // Check dimensionality.
  int expectedDim(0);
  if (testid == "slab")
    expectedDim = 1;
  if (testid == "quad")
    expectedDim = 2;
  if (testid == "cube")
    expectedDim = 3;

  bool pass_gdn(mesh.get_dims_ndim() == static_cast<size_t>(expectedDim));
  if (pass_gdn) {
    ostringstream message;
    message << "For " << testid << " mesh, mesh.get_dims_ndim() reported " << expectedDim << ".";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "For " << testid << " mesh, mesh.get_dims_ndim() incorrectly reported "
            << expectedDim << ".";
    FAILMSG(message.str());
  }

  return pass_gdn;
}

bool check_element_nodes(UnitTest &ut, const rtt_meshReaders::Hex_Mesh_Reader &mesh,
                         const std::string &testid) {
  // Check element nodes.
  bool pass_en = true;
  std::vector<std::vector<unsigned>> enodes = mesh.get_element_nodes();
  if (testid == "slab")
    pass_en = enodes[0][0] == 0 && enodes[0][1] == 1 && enodes[10][0] == 10 &&
              enodes[10][1] == 11 && enodes[99][0] == 99 && enodes[99][1] == 100 &&
              enodes[100][0] == 0 && enodes[101][0] == 100;
  else if (testid == "quad")
    pass_en = enodes[0][0] == 0 && enodes[0][1] == 1 && enodes[0][2] == 22 && enodes[0][3] == 21 &&
              enodes[10][0] == 10 && enodes[10][1] == 11 && enodes[10][2] == 32 &&
              enodes[10][3] == 31 && enodes[399][0] == 418 && enodes[399][1] == 419 &&
              enodes[399][2] == 440 && enodes[399][3] == 439 && enodes[400][0] == 421 &&
              enodes[400][1] == 420 && enodes[419][0] == 440 && enodes[419][1] == 439 &&
              enodes[420][0] == 21 && enodes[420][1] == 0 && enodes[479][0] == 19 &&
              enodes[479][1] == 20;
  else if (testid == "cube") {
    bool pass1 = enodes[0][0] == 0 && enodes[0][1] == 1 && enodes[0][2] == 7 && enodes[0][3] == 6 &&
                 enodes[0][4] == 36 && enodes[0][5] == 37 && enodes[0][6] == 43 &&
                 enodes[0][7] == 42;
    bool pass2 = enodes[10][0] == 12 && enodes[10][1] == 13 && enodes[10][2] == 19 &&
                 enodes[10][3] == 18 && enodes[10][4] == 48 && enodes[10][5] == 49 &&
                 enodes[10][6] == 55 && enodes[10][7] == 54;
    bool pass3 = enodes[124][0] == 172 && enodes[124][1] == 173 && enodes[124][2] == 179 &&
                 enodes[124][3] == 178 && enodes[124][4] == 208 && enodes[124][5] == 209 &&
                 enodes[124][6] == 215 && enodes[124][7] == 214;
    bool pass4 =
        enodes[125][0] == 5 && enodes[125][1] == 11 && enodes[125][2] == 47 && enodes[125][3] == 41;
    bool pass5 = enodes[249][0] == 208 && enodes[249][1] == 209 && enodes[249][2] == 215 &&
                 enodes[249][3] == 214;
    bool pass6 =
        enodes[250][0] == 0 && enodes[250][1] == 36 && enodes[250][2] == 42 && enodes[250][3] == 6;
    bool pass7 = enodes[274][0] == 168 && enodes[274][1] == 204 && enodes[274][2] == 210 &&
                 enodes[274][3] == 174;
    pass_en = pass1 && pass2 && pass3 && pass4 && pass5 && pass6 && pass7;
  } else
    Insist(false, "Unrecognized test id string!");
  if (pass_en) {
    ostringstream message;
    message << "Got element nodes." << std::endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in element nodes." << std::endl;
    FAILMSG(message.str());
  }
  return pass_en;
}

bool check_invariant(UnitTest &ut, const rtt_meshReaders::Hex_Mesh_Reader &mesh) {
  // Check invariant.
  bool invr = mesh.invariant();
  if (invr) {
    ostringstream message;
    message << "Invoked invariant." << std::endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Invariant not satisfied." << std::endl;
    FAILMSG(message.str());
  }
  return invr;
}

bool check_element_sets(UnitTest &ut, const rtt_meshReaders::Hex_Mesh_Reader &mesh,
                        const std::string &testid) {
  using mt = std::map<std::string, std::set<unsigned>>;
  bool pass_es = true;
  const mt elmsets = mesh.get_element_sets();
  if (testid == "slab") {
    pass_es = pass_es && elmsets.size() == 4;
    pass_es = pass_es && check_map(elmsets, "Interior", 0, 100);
    pass_es = pass_es && check_map(elmsets, "Interior_Region_1", 0, 100);
    pass_es = pass_es && check_map(elmsets, "Vacuum_Boundary", 100, 102);
    pass_es = pass_es && check_map(elmsets, "Vacuum_Boundary_Region_1", 100, 102);
  } else if (testid == "quad") {
    pass_es = pass_es && elmsets.size() == 6;
    pass_es = pass_es && check_map(elmsets, "Interior", 0, 400);
    pass_es = pass_es && check_map(elmsets, "Interior_Region_1", 20, 400);
    pass_es = pass_es && check_map(elmsets, "Interior_Region_2", 0, 20);
    pass_es = pass_es && check_map(elmsets, "Vacuum_Boundary", 400, 420);
    pass_es = pass_es && check_map(elmsets, "Vacuum_Boundary_Region_1", 400, 420);
    pass_es = pass_es && check_map(elmsets, "Reflective_Boundary", 420, 480);

  } else if (testid == "cube") {
    bool pass0 = elmsets.size() == 10;
    bool pass1 = check_map(elmsets, "Interior", 0, 125);
    bool pass2 = check_map(elmsets, "Interior_Region_10", 0, 2);
    bool pass3 = check_map(elmsets, "Interior_Region_14", 2, 3);
    bool pass4 = check_map(elmsets, "Interior_Region_1", 3, 122);
    bool pass5 = check_map(elmsets, "Interior_Region_2", 122, 125);
    bool pass6 = check_map(elmsets, "Vacuum_Boundary", 125, 250);
    bool pass7 = check_map(elmsets, "Vacuum_Boundary_Region_3", 125, 130);
    bool pass8 = check_map(elmsets, "Vacuum_Boundary_Region_4", 130, 131);
    bool pass9 = check_map(elmsets, "Vacuum_Boundary_Region_1", 131, 250);
    bool pass10 = check_map(elmsets, "Reflective_Boundary", 250, 275);
    pass_es = pass0 && pass1 && pass2 && pass3 && pass4 && pass5 && pass6 && pass7 && pass8 &&
              pass9 && pass10;
  } else
    Insist(false, "Unrecognized test id string!");
  if (pass_es) {
    ostringstream message;
    message << "Got element sets." << std::endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in element sets." << std::endl;
    FAILMSG(message.str());
  }

  return pass_es;
}

bool check_element_types(UnitTest &ut, const rtt_meshReaders::Hex_Mesh_Reader &mesh,
                         const std::string &testid) {
  // Check Element Types.
  using et = rtt_mesh_element::Element_Definition;
  bool pass_et = true;
  std::vector<rtt_mesh_element::Element_Definition::Element_Type> etypes = mesh.get_element_types();
  if (testid == "slab") {
    for (int i = 0; i < 100; ++i)
      pass_et = pass_et && (etypes[i] == et::BAR_2);
    for (int i = 100; i < 102; ++i)
      etypes[i] = et::NODE;
  } else if (testid == "quad") {
    for (int i = 0; i < 400; ++i)
      pass_et = pass_et && (etypes[i] == et::QUAD_4);
    for (int i = 400; i < 480; ++i)
      etypes[i] = et::BAR_2;
  } else if (testid == "cube") {
    for (int i = 0; i < 125; ++i)
      pass_et = pass_et && (etypes[i] == et::HEXA_8);
    for (int i = 125; i < 275; ++i)
      etypes[i] = et::QUAD_4;
  } else
    Insist(false, "Unrecognized test id string!");
  if (pass_et) {
    ostringstream message;
    message << "Read Element Types." << std::endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in Element Types." << std::endl;
    FAILMSG(message.str());
  }
  return pass_et;
}

bool check_unique_element_types(UnitTest &ut, const rtt_meshReaders::Hex_Mesh_Reader &mesh,
                                const std::string &testid) {
  // Check Unique Element Types.
  using et = rtt_mesh_element::Element_Definition;
  bool pass_et = true;
  std::vector<rtt_mesh_element::Element_Definition::Element_Type> etypes =
      mesh.get_unique_element_types();
  if (testid == "slab") {
    pass_et = pass_et && (etypes[0] == et::NODE && etypes[1] == et::BAR_2);
  } else if (testid == "quad") {
    pass_et =
        pass_et && (etypes[0] == et::NODE && etypes[1] == et::BAR_2 && etypes[2] == et::QUAD_4);
  } else if (testid == "cube") {
    pass_et = pass_et && (etypes[0] == et::NODE && etypes[1] == et::BAR_2 &&
                          etypes[2] == et::QUAD_4 && etypes[3] == et::HEXA_8);
  } else
    Insist(false, "Unrecognized test id string!");
  if (pass_et) {
    ostringstream message;
    message << "Read Unique Element Types." << std::endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in Unique Element Types." << std::endl;
    FAILMSG(message.str());
  }
  return pass_et;
}

bool compare_double(const double &lhs, const double &rhs) {
  // Note that this is only good for doubles close to one.
  return std::fabs(lhs - rhs) <= 0.00001;
}

bool check_map(const std::map<std::string, std::set<unsigned>> &elmsets, const std::string &name,
               const unsigned &begin, const unsigned &end) {
  bool pass = true;
  auto const iter = elmsets.find(name);
  if (iter != elmsets.end()) {
    const std::set<unsigned> &elem_subset = (*iter).second;
    if (elem_subset.size() == static_cast<size_t>(end - begin)) {
      for (unsigned i = begin; i < end; ++i) {
        auto const siter = elem_subset.find(i);
        pass = pass && (siter != elem_subset.end());
      }
    } else
      pass = false;
  } else
    pass = false;
  return pass;
}
} // namespace rtt_meshReaders_test

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    runTest(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of TestHexMeshReader.cc
//------------------------------------------------------------------------------------------------//
