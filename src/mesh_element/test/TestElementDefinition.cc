//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh_element/test/TestElementDefinition.cc
 * \author Thomas M. Evans
 * \date   Tue Mar 26 16:06:55 2002
 * \brief  Test Element Definitions.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "mesh_element/Element_Definition.hh"
#include <array>
#include <list>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

namespace rtt_mesh_element_test {

bool test_node(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_bar_2(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_bar_3(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_tri_3(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_tri_6(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_quad_4(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_quad_5(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_quad_6(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_quad_6a(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_quad_6o(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_quad_7(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_quad_8(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_quad_9(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_tetra_4(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_tetra_10(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_pyra_5(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_pyra_14(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_penta_6(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_penta_15(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_penta_18(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_hexa_8(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_hexa_20(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);
bool test_hexa_27(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def);

} // end namespace rtt_mesh_element_test

//------------------------------------------------------------------------------------------------//
template <std::size_t N> std::vector<unsigned> makeuv(std::array<int, N> s) {
  return std::vector<unsigned>(s.begin(), s.end());
}

//------------------------------------------------------------------------------------------------//
void runTest(rtt_dsxx::UnitTest &ut) {
  using rtt_mesh_element::Element_Definition;

  vector<Element_Definition::Element_Type> type_list;
  type_list.push_back(Element_Definition::NODE);
  type_list.push_back(Element_Definition::BAR_2);
  type_list.push_back(Element_Definition::BAR_3);
  type_list.push_back(Element_Definition::TRI_3);
  type_list.push_back(Element_Definition::TRI_6);
  type_list.push_back(Element_Definition::QUAD_4);
  type_list.push_back(Element_Definition::QUAD_5);
  type_list.push_back(
      Element_Definition::QUAD_6); // Call all three QUAD_6 elements in a single test
  type_list.push_back(Element_Definition::QUAD_6a);
  type_list.push_back(Element_Definition::QUAD_6o);
  type_list.push_back(Element_Definition::QUAD_7);
  type_list.push_back(Element_Definition::QUAD_8);
  type_list.push_back(Element_Definition::QUAD_9);
  type_list.push_back(Element_Definition::TETRA_4);
  type_list.push_back(Element_Definition::TETRA_10);
  type_list.push_back(Element_Definition::PYRA_5);
  type_list.push_back(Element_Definition::PYRA_14);
  type_list.push_back(Element_Definition::PENTA_6);
  type_list.push_back(Element_Definition::PENTA_15);
  type_list.push_back(Element_Definition::PENTA_18);
  type_list.push_back(Element_Definition::HEXA_8);
  type_list.push_back(Element_Definition::HEXA_20);
  type_list.push_back(Element_Definition::HEXA_27);

  //----------------------------------------------------------------------------------------------//
  cout << endl << "Building Elements for Test ---" << endl << endl;
  list<Element_Definition> elem_defs;
  for (size_t i = 0; i < type_list.size(); i++) {
    elem_defs.emplace_back(type_list[i]);
    cout << elem_defs.back();
    if (!elem_defs.back().invariant_satisfied()) {
      ostringstream msg;
      msg << "invariant_satisfied() failed for element i=" << i
          << ", whose type is = " << elem_defs.back().get_name() << std::endl;
      FAILMSG(msg.str());
    }
  }
  //----------------------------------------------------------------------------------------------//
  cout << "\nChecking Elements ---\n" << endl;

  // CAREFUL HERE -- the order of the function calls must match the type_list order
  rtt_mesh_element_test::test_node(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_bar_2(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_bar_3(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_tri_3(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_tri_6(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_quad_4(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_quad_5(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_quad_6(
      ut, elem_defs.front()); // Call all three QUAD_6 elements in a single test
  elem_defs.pop_front();
  rtt_mesh_element_test::test_quad_6(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_quad_6(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_quad_7(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_quad_8(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_quad_9(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_tetra_4(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_tetra_10(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_pyra_5(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_pyra_14(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_penta_6(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_penta_15(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_penta_18(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_hexa_8(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_hexa_20(ut, elem_defs.front());
  elem_defs.pop_front();
  rtt_mesh_element_test::test_hexa_27(ut, elem_defs.front());
  elem_defs.pop_front();

  //----------------------------------------------------------------------------------------------//
  // Test the POLYGON element.
  //----------------------------------------------------------------------------------------------//

  cout << "\nBuilding POLYGON Element ---\n" << endl;

  vector<Element_Definition> polyg_elem_defs;
  polyg_elem_defs.emplace_back(Element_Definition::BAR_2);

  vector<unsigned> polyg_side_type(8,  // number of sides
                                   0); // index into elem_defs

  vector<vector<unsigned>> polyg_side_nodes(8);
  for (unsigned side = 0; side < 8; ++side) {
    polyg_side_nodes[side].push_back(side);
    polyg_side_nodes[side].push_back((side + 1) % 8);
  }

  Element_Definition polyg("OCT_8", // name
                           2,       // dimension
                           8,       // number_of_nodes
                           8,       // number_of_sides
                           polyg_elem_defs, polyg_side_type, polyg_side_nodes);

  //----------------------------------------------------------------------------------------------//
  // Test the POLYHEDRON element.
  //----------------------------------------------------------------------------------------------//

  cout << "\nBuilding POLYHEDRON Element ---\n" << endl;

  vector<Element_Definition> polyh_elem_defs;
  vector<unsigned> polyh_side_type;
  vector<vector<unsigned>> polyh_side_nodes;

  // First side is QUAD_4
  polyh_elem_defs.emplace_back(Element_Definition::QUAD_4);
  polyh_side_type.push_back(0);
  vector<unsigned> side0_nodes = {0, 1, 5, 4};
  polyh_side_nodes.push_back(side0_nodes);

  // Next four sides are QUAD_5
  polyh_elem_defs.emplace_back(Element_Definition::QUAD_5);

  polyh_side_type.push_back(1);
  vector<unsigned> side1_nodes = {1, 0, 3, 2, 11};
  polyh_side_nodes.push_back(side1_nodes);

  polyh_side_type.push_back(1);
  vector<unsigned> side2_nodes = {5, 1, 2, 6, 10};
  polyh_side_nodes.push_back(side2_nodes);

  polyh_side_type.push_back(1);
  vector<unsigned> side3_nodes = {4, 5, 6, 7, 9};
  polyh_side_nodes.push_back(side3_nodes);

  polyh_side_type.push_back(1);
  vector<unsigned> side4_nodes = {0, 4, 7, 3, 8};
  polyh_side_nodes.push_back(side4_nodes);

  // Last (sixth) side is QUAD_9
  polyh_elem_defs.emplace_back(Element_Definition::QUAD_9);
  polyh_side_type.push_back(2);
  vector<unsigned> side5_nodes = {3, 7, 6, 2, 8, 9, 10, 11, 12};
  polyh_side_nodes.push_back(side5_nodes);

  Element_Definition polyh("polyhedron1", // name
                           3,             // dimension
                           13,            // number_of_nodes
                           6,             // number_of_sides
                           polyh_elem_defs, polyh_side_type, polyh_side_nodes);

  //----------------------------------------------------------------------------------------------//
  // Merely attempting construction, with DBC active, will invoke a slew of precondition,
  // postcondition, and consistency checks.  We perform no other explicit checks here.

  if (ut.numFails == 0)
    PASSMSG("All tests passed.");
  else
    FAILMSG("Some tests failed.");
  return;
}

//------------------------------------------------------------------------------------------------//

namespace rtt_mesh_element_test {

//------------------------------------------------------------------------------------------------//
bool test_node(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the NODE element.
  using rtt_mesh_element::Element_Definition;
  string ename = "NODE";
  bool ldum = elem_def.get_name() == ename;
  ldum = ldum && elem_def.get_type() == Element_Definition::NODE;
  ldum = ldum && elem_def.get_number_of_nodes() == 1;
  ldum = ldum && elem_def.get_dimension() == 0;
  ldum = ldum && elem_def.get_number_of_face_nodes().size() == 0;
  if (ldum) {
    ostringstream message;
    message << ename << " Element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in " << ename << " Element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}

//------------------------------------------------------------------------------------------------//
bool test_bar_2(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the BAR_2 element.
  using rtt_mesh_element::Element_Definition;
  string ename = "BAR_2";
  bool ldum = elem_def.get_name() == ename;
  ldum = ldum && elem_def.get_type() == Element_Definition::BAR_2;
  ldum = ldum && elem_def.get_number_of_nodes() == 2;
  ldum = ldum && elem_def.get_dimension() == 1;
  ldum = ldum && elem_def.get_number_of_sides() == 2;
  ldum = ldum && elem_def.get_number_of_face_nodes().size() == 2;
  ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 1;
  ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 1;
  constexpr int size = 1;
  std::array<int, size> s0 = {0};
  std::array<int, size> s1 = {1};
  ldum = ldum && (elem_def.get_side_nodes(0) == makeuv(s0));
  ldum = ldum && (elem_def.get_side_nodes(1) == makeuv(s1));
  if (ldum) {
    ostringstream message;
    message << ename << " Element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in " << ename << " Element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}

//------------------------------------------------------------------------------------------------//
bool test_bar_3(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the BAR_3 element.
  using rtt_mesh_element::Element_Definition;
  string ename = "BAR_3";
  bool ldum = elem_def.get_name() == ename;
  ldum = ldum && elem_def.get_type() == Element_Definition::BAR_3;
  ldum = ldum && elem_def.get_number_of_nodes() == 3;
  ldum = ldum && elem_def.get_dimension() == 1;
  ldum = ldum && elem_def.get_number_of_sides() == 2;
  ldum = ldum && elem_def.get_number_of_face_nodes().size() == 2;
  ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 1;
  ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 1;
  constexpr int size = 1;
  std::array<int, size> s0 = {0};
  std::array<int, size> s1 = {1};
  ldum = ldum && (elem_def.get_side_nodes(0) == makeuv(s0));
  ldum = ldum && (elem_def.get_side_nodes(1) == makeuv(s1));
  if (ldum) {
    ostringstream message;
    message << ename << " Element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in " << ename << " Element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}

//------------------------------------------------------------------------------------------------//
bool test_tri_3(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the TRI_3 element.
  using rtt_mesh_element::Element_Definition;
  string ename = "TRI_3";
  bool ldum = elem_def.get_name() == ename;
  ldum = ldum && elem_def.get_type() == Element_Definition::TRI_3;
  ldum = ldum && elem_def.get_number_of_nodes() == 3;
  ldum = ldum && elem_def.get_dimension() == 2;
  ldum = ldum && elem_def.get_number_of_sides() == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes().size() == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 2;
  ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 2;
  ldum = ldum && elem_def.get_number_of_face_nodes()[2] == 2;
  constexpr int size = 2;
  std::array<int, size> s0 = {0, 1};
  std::array<int, size> s1 = {1, 2};
  std::array<int, size> s2 = {2, 0};
  ldum = ldum && elem_def.get_side_nodes(0) == makeuv(s0);
  ldum = ldum && elem_def.get_side_nodes(1) == makeuv(s1);
  ldum = ldum && elem_def.get_side_nodes(2) == makeuv(s2);
  if (ldum) {
    ostringstream message;
    message << ename << " Element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in " << ename << " Element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}

//------------------------------------------------------------------------------------------------//
bool test_tri_6(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the element.
  using rtt_mesh_element::Element_Definition;
  string ename = "TRI_6";
  bool ldum = elem_def.get_name() == ename;
  ldum = ldum && elem_def.get_type() == Element_Definition::TRI_6;
  ldum = ldum && elem_def.get_number_of_nodes() == 6;
  ldum = ldum && elem_def.get_dimension() == 2;
  ldum = ldum && elem_def.get_number_of_sides() == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes().size() == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes()[2] == 3;
  constexpr int size = 3;
  std::array<int, size> s0 = {0, 1, 3};
  std::array<int, size> s1 = {1, 2, 4};
  std::array<int, size> s2 = {2, 0, 5};
  ldum = ldum && elem_def.get_side_nodes(0) == makeuv(s0);
  ldum = ldum && elem_def.get_side_nodes(1) == makeuv(s1);
  ldum = ldum && elem_def.get_side_nodes(2) == makeuv(s2);
  if (ldum) {
    ostringstream message;
    message << ename << " Element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in " << ename << " Element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}

//------------------------------------------------------------------------------------------------//
bool test_quad_4(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the QUAD_4 element.
  using rtt_mesh_element::Element_Definition;
  string ename = "QUAD_4";
  bool ldum = elem_def.get_name() == ename;
  ldum = ldum && elem_def.get_type() == Element_Definition::QUAD_4;
  ldum = ldum && elem_def.get_number_of_nodes() == 4;
  ldum = ldum && elem_def.get_dimension() == 2;
  ldum = ldum && elem_def.get_number_of_sides() == 4;
  ldum = ldum && elem_def.get_number_of_face_nodes().size() == 4;
  ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 2;
  ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 2;
  ldum = ldum && elem_def.get_number_of_face_nodes()[2] == 2;
  ldum = ldum && elem_def.get_number_of_face_nodes()[3] == 2;
  constexpr int size = 2;
  std::array<int, size> s0 = {0, 1};
  std::array<int, size> s1 = {1, 2};
  std::array<int, size> s2 = {2, 3};
  std::array<int, size> s3 = {3, 0};
  ldum = ldum && elem_def.get_side_nodes(0) == makeuv(s0);
  ldum = ldum && elem_def.get_side_nodes(1) == makeuv(s1);
  ldum = ldum && elem_def.get_side_nodes(2) == makeuv(s2);
  ldum = ldum && elem_def.get_side_nodes(3) == makeuv(s3);
  if (ldum) {
    ostringstream message;
    message << ename << " Element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in " << ename << " Element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}

//------------------------------------------------------------------------------------------------//
bool test_quad_5(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the QUAD_5 element.
  using rtt_mesh_element::Element_Definition;
  string ename = "QUAD_5";
  bool ldum = elem_def.get_name() == ename;

  ldum = ldum && elem_def.get_type() == Element_Definition::QUAD_5;
  ldum = ldum && elem_def.get_number_of_nodes() == 5;
  ldum = ldum && elem_def.get_dimension() == 2;
  ldum = ldum && elem_def.get_number_of_sides() == 4;

  ldum = ldum && elem_def.get_number_of_face_nodes().size() == 4;

  ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 2;
  ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 2;
  ldum = ldum && elem_def.get_number_of_face_nodes()[2] == 2;
  ldum = ldum && elem_def.get_number_of_face_nodes()[3] == 3;

  std::array<int, 2> s0 = {0, 1};
  std::array<int, 2> s1 = {1, 2};
  std::array<int, 2> s2 = {2, 3};
  std::array<int, 3> s3 = {3, 0, 4};
  ldum = ldum && elem_def.get_side_nodes(0) == makeuv(s0);
  ldum = ldum && elem_def.get_side_nodes(1) == makeuv(s1);
  ldum = ldum && elem_def.get_side_nodes(2) == makeuv(s2);
  ldum = ldum && elem_def.get_side_nodes(3) == makeuv(s3);

  if (ldum) {
    ostringstream message;
    message << ename << " Element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in " << ename << " Element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}

//------------------------------------------------------------------------------------------------//
bool test_quad_6(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the QUAD_6 element.
  using rtt_mesh_element::Element_Definition;

  bool ldum = (elem_def.get_type() == Element_Definition::QUAD_6 ||
               elem_def.get_type() == Element_Definition::QUAD_6a ||
               elem_def.get_type() == Element_Definition::QUAD_6o);

  if (elem_def.get_type() == Element_Definition::QUAD_6 ||
      elem_def.get_type() == Element_Definition::QUAD_6a) {
    string ename = "QUAD_6";
    ldum = ldum && elem_def.get_name() == ename;

    ldum = ldum && elem_def.get_number_of_nodes() == 6;
    ldum = ldum && elem_def.get_dimension() == 2;
    ldum = ldum && elem_def.get_number_of_sides() == 4;

    ldum = ldum && elem_def.get_number_of_face_nodes().size() == 4;

    ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 2;
    ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 2;
    ldum = ldum && elem_def.get_number_of_face_nodes()[2] == 3;
    ldum = ldum && elem_def.get_number_of_face_nodes()[3] == 3;

    std::array<int, 2> s0 = {0, 1};
    std::array<int, 2> s1 = {1, 2};
    std::array<int, 3> s2 = {2, 3, 4};
    std::array<int, 3> s3 = {3, 0, 5};
    ldum = ldum && elem_def.get_side_nodes(0) == makeuv(s0);
    ldum = ldum && elem_def.get_side_nodes(1) == makeuv(s1);
    ldum = ldum && elem_def.get_side_nodes(2) == makeuv(s2);
    ldum = ldum && elem_def.get_side_nodes(3) == makeuv(s3);
  } else {
    string ename = "QUAD_6o";
    ldum = ldum && elem_def.get_name() == ename;

    ldum = ldum && elem_def.get_number_of_nodes() == 6;
    ldum = ldum && elem_def.get_dimension() == 2;
    ldum = ldum && elem_def.get_number_of_sides() == 4;

    ldum = ldum && elem_def.get_number_of_face_nodes().size() == 4;

    ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 2;
    ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 3;
    ldum = ldum && elem_def.get_number_of_face_nodes()[2] == 2;
    ldum = ldum && elem_def.get_number_of_face_nodes()[3] == 3;

    std::array<int, 2> s0 = {0, 1};
    std::array<int, 3> s1 = {1, 2, 4};
    std::array<int, 2> s2 = {2, 3};
    std::array<int, 3> s3 = {3, 0, 5};
    ldum = ldum && elem_def.get_side_nodes(0) == makeuv(s0);
    ldum = ldum && elem_def.get_side_nodes(1) == makeuv(s1);
    ldum = ldum && elem_def.get_side_nodes(2) == makeuv(s2);
    ldum = ldum && elem_def.get_side_nodes(3) == makeuv(s3);
  }

  if (ldum) {
    ostringstream message;
    message << " QUAD_6 element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in QUAD_6 element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}

//------------------------------------------------------------------------------------------------//
bool test_quad_7(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the QUAD_7 element.
  using rtt_mesh_element::Element_Definition;
  string ename = "QUAD_7";
  bool ldum = elem_def.get_name() == ename;

  ldum = ldum && elem_def.get_type() == Element_Definition::QUAD_7;
  ldum = ldum && elem_def.get_number_of_nodes() == 7;
  ldum = ldum && elem_def.get_dimension() == 2;
  ldum = ldum && elem_def.get_number_of_sides() == 4;

  ldum = ldum && elem_def.get_number_of_face_nodes().size() == 4;

  ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 2;
  ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes()[2] == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes()[3] == 3;

  std::array<int, 2> s0 = {0, 1};
  std::array<int, 3> s1 = {1, 2, 4};
  std::array<int, 3> s2 = {2, 3, 5};
  std::array<int, 3> s3 = {3, 0, 6};
  ldum = ldum && elem_def.get_side_nodes(0) == makeuv(s0);
  ldum = ldum && elem_def.get_side_nodes(1) == makeuv(s1);
  ldum = ldum && elem_def.get_side_nodes(2) == makeuv(s2);
  ldum = ldum && elem_def.get_side_nodes(3) == makeuv(s3);

  if (ldum) {
    ostringstream message;
    message << ename << " Element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in " << ename << " Element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}

//------------------------------------------------------------------------------------------------//
bool test_quad_8(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the QUAD_8 element.
  using rtt_mesh_element::Element_Definition;
  string ename = "QUAD_8";
  bool ldum = elem_def.get_name() == ename;
  ldum = ldum && elem_def.get_type() == Element_Definition::QUAD_8;
  ldum = ldum && elem_def.get_number_of_nodes() == 8;
  ldum = ldum && elem_def.get_dimension() == 2;
  ldum = ldum && elem_def.get_number_of_sides() == 4;
  ldum = ldum && elem_def.get_number_of_face_nodes().size() == 4;
  ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes()[2] == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes()[3] == 3;
  constexpr int size = 3;
  std::array<int, size> s0 = {0, 1, 4};
  std::array<int, size> s1 = {1, 2, 5};
  std::array<int, size> s2 = {2, 3, 6};
  std::array<int, size> s3 = {3, 0, 7};
  ldum = ldum && elem_def.get_side_nodes(0) == makeuv(s0);
  ldum = ldum && elem_def.get_side_nodes(1) == makeuv(s1);
  ldum = ldum && elem_def.get_side_nodes(2) == makeuv(s2);
  ldum = ldum && elem_def.get_side_nodes(3) == makeuv(s3);
  if (ldum) {
    ostringstream message;
    message << ename << " Element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in " << ename << " Element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}

//------------------------------------------------------------------------------------------------//
bool test_quad_9(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the QUAD_9 element.
  using rtt_mesh_element::Element_Definition;
  string ename = "QUAD_9";
  bool ldum = elem_def.get_name() == ename;
  ldum = ldum && elem_def.get_type() == Element_Definition::QUAD_9;
  ldum = ldum && elem_def.get_number_of_nodes() == 9;
  ldum = ldum && elem_def.get_dimension() == 2;
  ldum = ldum && elem_def.get_number_of_sides() == 4;
  ldum = ldum && elem_def.get_number_of_face_nodes().size() == 4;
  ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes()[2] == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes()[3] == 3;
  constexpr int size = 3;
  array<int, size> s0 = {0, 1, 4};
  array<int, size> s1 = {1, 2, 5};
  array<int, size> s2 = {2, 3, 6};
  array<int, size> s3 = {3, 0, 7};
  ldum = ldum && elem_def.get_side_nodes(0) == makeuv(s0);
  ldum = ldum && elem_def.get_side_nodes(1) == makeuv(s1);
  ldum = ldum && elem_def.get_side_nodes(2) == makeuv(s2);
  ldum = ldum && elem_def.get_side_nodes(3) == makeuv(s3);
  if (ldum) {
    ostringstream message;
    message << ename << " Element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in " << ename << " Element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}

//------------------------------------------------------------------------------------------------//
bool test_tetra_4(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the TETRA_4 element.
  using rtt_mesh_element::Element_Definition;
  string ename = "TETRA_4";
  bool ldum = elem_def.get_name() == ename;
  ldum = ldum && elem_def.get_type() == Element_Definition::TETRA_4;
  ldum = ldum && elem_def.get_number_of_nodes() == 4;
  ldum = ldum && elem_def.get_dimension() == 3;
  ldum = ldum && elem_def.get_number_of_sides() == 4;
  ldum = ldum && elem_def.get_number_of_face_nodes().size() == 4;
  ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes()[2] == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes()[3] == 3;
  constexpr int size = 3;
  array<int, size> s0 = {0, 2, 1};
  array<int, size> s1 = {0, 1, 3};
  array<int, size> s2 = {1, 2, 3};
  array<int, size> s3 = {2, 0, 3};
  ldum = ldum && elem_def.get_side_nodes(0) == makeuv(s0);
  ldum = ldum && elem_def.get_side_nodes(1) == makeuv(s1);
  ldum = ldum && elem_def.get_side_nodes(2) == makeuv(s2);
  ldum = ldum && elem_def.get_side_nodes(3) == makeuv(s3);
  if (ldum) {
    ostringstream message;
    message << ename << " Element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in " << ename << " Element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}

//------------------------------------------------------------------------------------------------//
bool test_tetra_10(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the TETRA_10 element.
  using rtt_mesh_element::Element_Definition;
  string ename = "TETRA_10";
  bool ldum = elem_def.get_name() == ename;
  ldum = ldum && elem_def.get_type() == Element_Definition::TETRA_10;
  ldum = ldum && elem_def.get_number_of_nodes() == 10;
  ldum = ldum && elem_def.get_dimension() == 3;
  ldum = ldum && elem_def.get_number_of_sides() == 4;
  ldum = ldum && elem_def.get_number_of_face_nodes().size() == 4;
  ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 6;
  ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 6;
  ldum = ldum && elem_def.get_number_of_face_nodes()[2] == 6;
  ldum = ldum && elem_def.get_number_of_face_nodes()[3] == 6;
  constexpr int size = 6;
  array<int, size> s0 = {0, 2, 1, 6, 5, 4};
  array<int, size> s1 = {0, 1, 3, 4, 8, 7};
  array<int, size> s2 = {1, 2, 3, 5, 9, 8};
  array<int, size> s3 = {2, 0, 3, 6, 7, 9};
  ldum = ldum && elem_def.get_side_nodes(0) == makeuv(s0);
  ldum = ldum && elem_def.get_side_nodes(1) == makeuv(s1);
  ldum = ldum && elem_def.get_side_nodes(2) == makeuv(s2);
  ldum = ldum && elem_def.get_side_nodes(3) == makeuv(s3);
  if (ldum) {
    ostringstream message;
    message << ename << " Element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in " << ename << " Element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}

//------------------------------------------------------------------------------------------------//
bool test_pyra_5(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the PYRA_5 element.
  using rtt_mesh_element::Element_Definition;
  string ename = "PYRA_5";
  bool ldum = elem_def.get_name() == ename;
  ldum = ldum && elem_def.get_type() == Element_Definition::PYRA_5;
  ldum = ldum && elem_def.get_number_of_nodes() == 5;
  ldum = ldum && elem_def.get_dimension() == 3;
  ldum = ldum && elem_def.get_number_of_sides() == 5;
  ldum = ldum && elem_def.get_side_type(0).get_type() == Element_Definition::QUAD_4;
  ldum = ldum && elem_def.get_number_of_face_nodes().size() == 5;
  ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 4;
  ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes()[2] == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes()[3] == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes()[4] == 3;
  constexpr int sizeq = 4;
  array<int, sizeq> s0 = {0, 3, 2, 1};
  constexpr int sizet = 3;
  array<int, sizet> s1 = {0, 1, 4};
  array<int, sizet> s2 = {1, 2, 4};
  array<int, sizet> s3 = {2, 3, 4};
  array<int, sizet> s4 = {3, 0, 4};
  ldum = ldum && elem_def.get_side_nodes(0) == makeuv(s0);
  ldum = ldum && elem_def.get_side_nodes(1) == makeuv(s1);
  ldum = ldum && elem_def.get_side_nodes(2) == makeuv(s2);
  ldum = ldum && elem_def.get_side_nodes(3) == makeuv(s3);
  ldum = ldum && elem_def.get_side_nodes(4) == makeuv(s4);
  if (ldum) {
    ostringstream message;
    message << ename << " Element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in " << ename << " Element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}

//------------------------------------------------------------------------------------------------//
bool test_pyra_14(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the PYRA_14 element.
  using rtt_mesh_element::Element_Definition;
  string ename = "PYRA_14";
  bool ldum = elem_def.get_name() == ename;
  ldum = ldum && elem_def.get_type() == Element_Definition::PYRA_14;
  ldum = ldum && elem_def.get_number_of_nodes() == 14;
  ldum = ldum && elem_def.get_dimension() == 3;
  ldum = ldum && elem_def.get_number_of_sides() == 5;
  ldum = ldum && elem_def.get_side_type(0).get_type() == Element_Definition::QUAD_8;
  constexpr int sizeq = 8;
  array<int, sizeq> s0 = {0, 3, 2, 1, 8, 7, 6, 5};
  constexpr int sizet = 6;
  array<int, sizet> s1 = {0, 1, 4, 5, 10, 9};
  array<int, sizet> s2 = {1, 2, 4, 6, 11, 10};
  array<int, sizet> s3 = {2, 3, 4, 7, 12, 11};
  array<int, sizet> s4 = {3, 0, 4, 8, 9, 12};
  ldum = ldum && elem_def.get_side_nodes(0) == makeuv(s0);
  ldum = ldum && elem_def.get_side_nodes(1) == makeuv(s1);
  ldum = ldum && elem_def.get_side_nodes(2) == makeuv(s2);
  ldum = ldum && elem_def.get_side_nodes(3) == makeuv(s3);
  ldum = ldum && elem_def.get_side_nodes(4) == makeuv(s4);
  ldum = ldum && elem_def.get_number_of_face_nodes().size() == 5;
  ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 8;
  ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 6;
  ldum = ldum && elem_def.get_number_of_face_nodes()[2] == 6;
  ldum = ldum && elem_def.get_number_of_face_nodes()[3] == 6;
  ldum = ldum && elem_def.get_number_of_face_nodes()[4] == 6;
  if (ldum) {
    ostringstream message;
    message << ename << " Element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in " << ename << " Element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}

//------------------------------------------------------------------------------------------------//
bool test_penta_6(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the PENTA_6 element.
  using rtt_mesh_element::Element_Definition;
  string ename = "PENTA_6";
  bool ldum = elem_def.get_name() == ename;
  ldum = ldum && elem_def.get_type() == Element_Definition::PENTA_6;
  ldum = ldum && elem_def.get_number_of_nodes() == 6;
  ldum = ldum && elem_def.get_dimension() == 3;
  ldum = ldum && elem_def.get_number_of_sides() == 5;
  for (int j = 0; j < 3; ++j)
    ldum = ldum && elem_def.get_side_type(j).get_type() == Element_Definition::QUAD_4;
  for (int j = 3; j < 5; ++j)
    ldum = ldum && elem_def.get_side_type(j).get_type() == Element_Definition::TRI_3;
  constexpr int sizeq = 4;
  array<int, sizeq> s0 = {0, 1, 4, 3};
  array<int, sizeq> s1 = {1, 2, 5, 4};
  array<int, sizeq> s2 = {2, 0, 3, 5};
  constexpr int sizet = 3;
  array<int, sizet> s3 = {0, 2, 1};
  array<int, sizet> s4 = {3, 4, 5};
  ldum = ldum && elem_def.get_side_nodes(0) == makeuv(s0);
  ldum = ldum && elem_def.get_side_nodes(1) == makeuv(s1);
  ldum = ldum && elem_def.get_side_nodes(2) == makeuv(s2);
  ldum = ldum && elem_def.get_side_nodes(3) == makeuv(s3);
  ldum = ldum && elem_def.get_side_nodes(4) == makeuv(s4);
  ldum = ldum && elem_def.get_number_of_face_nodes().size() == 5;
  ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 4;
  ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 4;
  ldum = ldum && elem_def.get_number_of_face_nodes()[2] == 4;
  ldum = ldum && elem_def.get_number_of_face_nodes()[3] == 3;
  ldum = ldum && elem_def.get_number_of_face_nodes()[4] == 3;
  if (ldum) {
    ostringstream message;
    message << ename << " Element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in " << ename << " Element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}

//------------------------------------------------------------------------------------------------//
bool test_penta_15(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the PENTA_15 element.
  using rtt_mesh_element::Element_Definition;
  string ename = "PENTA_15";
  bool ldum = elem_def.get_name() == ename;
  ldum = ldum && elem_def.get_type() == Element_Definition::PENTA_15;
  ldum = ldum && elem_def.get_number_of_nodes() == 15;
  ldum = ldum && elem_def.get_dimension() == 3;
  ldum = ldum && elem_def.get_number_of_sides() == 5;
  for (int j = 0; j < 3; ++j)
    ldum = ldum && elem_def.get_side_type(j).get_type() == Element_Definition::QUAD_8;
  for (int j = 3; j < 5; ++j)
    ldum = ldum && elem_def.get_side_type(j).get_type() == Element_Definition::TRI_6;
  constexpr int sizeq = 8;
  array<int, sizeq> s0 = {0, 1, 4, 3, 6, 10, 12, 9};
  array<int, sizeq> s1 = {1, 2, 5, 4, 7, 11, 13, 10};
  array<int, sizeq> s2 = {2, 0, 3, 5, 8, 9, 14, 11};
  constexpr int sizet = 6;
  array<int, sizet> s3 = {0, 2, 1, 8, 7, 6};
  array<int, sizet> s4 = {3, 4, 5, 12, 13, 14};
  ldum = ldum && elem_def.get_side_nodes(0) == makeuv(s0);
  ldum = ldum && elem_def.get_side_nodes(1) == makeuv(s1);
  ldum = ldum && elem_def.get_side_nodes(2) == makeuv(s2);
  ldum = ldum && elem_def.get_side_nodes(3) == makeuv(s3);
  ldum = ldum && elem_def.get_side_nodes(4) == makeuv(s4);
  ldum = ldum && elem_def.get_number_of_face_nodes().size() == 5;
  ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 8;
  ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 8;
  ldum = ldum && elem_def.get_number_of_face_nodes()[2] == 8;
  ldum = ldum && elem_def.get_number_of_face_nodes()[3] == 6;
  ldum = ldum && elem_def.get_number_of_face_nodes()[4] == 6;
  if (ldum) {
    ostringstream message;
    message << ename << " Element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in " << ename << " Element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}

//------------------------------------------------------------------------------------------------//
bool test_penta_18(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the PENTA_18 element.
  using rtt_mesh_element::Element_Definition;
  string ename = "PENTA_18";
  bool ldum = elem_def.get_name() == ename;
  ldum = ldum && elem_def.get_type() == Element_Definition::PENTA_18;
  ldum = ldum && elem_def.get_number_of_nodes() == 18;
  ldum = ldum && elem_def.get_dimension() == 3;
  ldum = ldum && elem_def.get_number_of_sides() == 5;
  for (int j = 0; j < 3; ++j)
    ldum = ldum && elem_def.get_side_type(j).get_type() == Element_Definition::QUAD_9;
  for (int j = 3; j < 5; ++j)
    ldum = ldum && elem_def.get_side_type(j).get_type() == Element_Definition::TRI_6;
  constexpr int sizeq = 9;
  array<int, sizeq> s0 = {0, 1, 4, 3, 6, 10, 12, 9, 15};
  array<int, sizeq> s1 = {1, 2, 5, 4, 7, 11, 13, 10, 16};
  array<int, sizeq> s2 = {2, 0, 3, 5, 8, 9, 14, 11, 17};
  constexpr int sizet = 6;
  array<int, sizet> s3 = {0, 2, 1, 8, 7, 6};
  array<int, sizet> s4 = {3, 4, 5, 12, 13, 14};
  ldum = ldum && elem_def.get_side_nodes(0) == makeuv(s0);
  ldum = ldum && elem_def.get_side_nodes(1) == makeuv(s1);
  ldum = ldum && elem_def.get_side_nodes(2) == makeuv(s2);
  ldum = ldum && elem_def.get_side_nodes(3) == makeuv(s3);
  ldum = ldum && elem_def.get_side_nodes(4) == makeuv(s4);
  ldum = ldum && elem_def.get_number_of_face_nodes().size() == 5;
  ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 9;
  ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 9;
  ldum = ldum && elem_def.get_number_of_face_nodes()[2] == 9;
  ldum = ldum && elem_def.get_number_of_face_nodes()[3] == 6;
  ldum = ldum && elem_def.get_number_of_face_nodes()[4] == 6;

  if (ldum) {
    ostringstream message;
    message << ename << " Element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in " << ename << " Element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}

//------------------------------------------------------------------------------------------------//
bool test_hexa_8(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the HEXA_8 element.
  using rtt_mesh_element::Element_Definition;
  string ename = "HEXA_8";
  bool ldum = elem_def.get_name() == ename;
  ldum = ldum && elem_def.get_type() == Element_Definition::HEXA_8;
  ldum = ldum && elem_def.get_number_of_nodes() == 8;
  ldum = ldum && elem_def.get_dimension() == 3;
  ldum = ldum && elem_def.get_number_of_sides() == 6;
  for (int j = 0; j < 6; ++j)
    ldum = ldum && elem_def.get_side_type(j).get_type() == Element_Definition::QUAD_4;
  constexpr int size = 4;
  array<int, size> s0 = {0, 3, 2, 1};
  array<int, size> s1 = {0, 4, 7, 3};
  array<int, size> s2 = {2, 3, 7, 6};
  array<int, size> s3 = {1, 2, 6, 5};
  array<int, size> s4 = {0, 1, 5, 4};
  array<int, size> s5 = {4, 5, 6, 7};
  ldum = ldum && elem_def.get_side_nodes(0) == makeuv(s0);
  ldum = ldum && elem_def.get_side_nodes(1) == makeuv(s1);
  ldum = ldum && elem_def.get_side_nodes(2) == makeuv(s2);
  ldum = ldum && elem_def.get_side_nodes(3) == makeuv(s3);
  ldum = ldum && elem_def.get_side_nodes(4) == makeuv(s4);
  ldum = ldum && elem_def.get_side_nodes(5) == makeuv(s5);
  ldum = ldum && elem_def.get_number_of_face_nodes().size() == 6;
  ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 4;
  ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 4;
  ldum = ldum && elem_def.get_number_of_face_nodes()[2] == 4;
  ldum = ldum && elem_def.get_number_of_face_nodes()[3] == 4;
  ldum = ldum && elem_def.get_number_of_face_nodes()[4] == 4;
  ldum = ldum && elem_def.get_number_of_face_nodes()[5] == 4;
  if (ldum) {
    ostringstream message;
    message << ename << " Element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in " << ename << " Element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}

//------------------------------------------------------------------------------------------------//
bool test_hexa_20(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the HEXA_20 element.
  using rtt_mesh_element::Element_Definition;
  string ename = "HEXA_20";
  bool ldum = elem_def.get_name() == ename;
  ldum = ldum && elem_def.get_type() == Element_Definition::HEXA_20;
  ldum = ldum && elem_def.get_number_of_nodes() == 20;
  ldum = ldum && elem_def.get_dimension() == 3;
  ldum = ldum && elem_def.get_number_of_sides() == 6;
  for (int j = 0; j < 6; ++j)
    ldum = ldum && elem_def.get_side_type(j).get_type() == Element_Definition::QUAD_8;
  constexpr int size = 8;
  array<int, size> s0 = {0, 3, 2, 1, 11, 10, 9, 8};
  array<int, size> s1 = {0, 4, 7, 3, 12, 19, 15, 11};
  array<int, size> s2 = {2, 3, 7, 6, 10, 15, 18, 14};
  array<int, size> s3 = {1, 2, 6, 5, 9, 14, 17, 13};
  array<int, size> s4 = {0, 1, 5, 4, 8, 13, 16, 12};
  array<int, size> s5 = {4, 5, 6, 7, 16, 17, 18, 19};
  ldum = ldum && elem_def.get_side_nodes(0) == makeuv(s0);
  ldum = ldum && elem_def.get_side_nodes(1) == makeuv(s1);
  ldum = ldum && elem_def.get_side_nodes(2) == makeuv(s2);
  ldum = ldum && elem_def.get_side_nodes(3) == makeuv(s3);
  ldum = ldum && elem_def.get_side_nodes(4) == makeuv(s4);
  ldum = ldum && elem_def.get_side_nodes(5) == makeuv(s5);
  ldum = ldum && elem_def.get_number_of_face_nodes().size() == 6;
  ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 8;
  ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 8;
  ldum = ldum && elem_def.get_number_of_face_nodes()[2] == 8;
  ldum = ldum && elem_def.get_number_of_face_nodes()[3] == 8;
  ldum = ldum && elem_def.get_number_of_face_nodes()[4] == 8;
  ldum = ldum && elem_def.get_number_of_face_nodes()[5] == 8;

  if (ldum) {
    ostringstream message;
    message << ename << " Element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in " << ename << " Element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}

//------------------------------------------------------------------------------------------------//
bool test_hexa_27(rtt_dsxx::UnitTest &ut, const rtt_mesh_element::Element_Definition &elem_def) {
  // Test the HEXA_27 element.
  using rtt_mesh_element::Element_Definition;
  string ename = "HEXA_27";
  bool ldum = elem_def.get_name() == ename;
  ldum = ldum && elem_def.get_type() == Element_Definition::HEXA_27;
  ldum = ldum && elem_def.get_number_of_nodes() == 27;
  ldum = ldum && elem_def.get_dimension() == 3;
  ldum = ldum && elem_def.get_number_of_sides() == 6;
  for (int j = 0; j < 6; ++j)
    ldum = ldum && elem_def.get_side_type(j).get_type() == Element_Definition::QUAD_9;
  constexpr int size = 9;
  array<int, size> s0 = {0, 3, 2, 1, 11, 10, 9, 8, 20};
  array<int, size> s1 = {0, 4, 7, 3, 12, 19, 15, 11, 21};
  array<int, size> s2 = {2, 3, 7, 6, 10, 15, 18, 14, 22};
  array<int, size> s3 = {1, 2, 6, 5, 9, 14, 17, 13, 23};
  array<int, size> s4 = {0, 1, 5, 4, 8, 13, 16, 12, 24};
  array<int, size> s5 = {4, 5, 6, 7, 16, 17, 18, 19, 25};
  ldum = ldum && elem_def.get_side_nodes(0) == makeuv(s0);
  ldum = ldum && elem_def.get_side_nodes(1) == makeuv(s1);
  ldum = ldum && elem_def.get_side_nodes(2) == makeuv(s2);
  ldum = ldum && elem_def.get_side_nodes(3) == makeuv(s3);
  ldum = ldum && elem_def.get_side_nodes(4) == makeuv(s4);
  ldum = ldum && elem_def.get_side_nodes(5) == makeuv(s5);
  ldum = ldum && elem_def.get_number_of_face_nodes().size() == 6;
  ldum = ldum && elem_def.get_number_of_face_nodes()[0] == 9;
  ldum = ldum && elem_def.get_number_of_face_nodes()[1] == 9;
  ldum = ldum && elem_def.get_number_of_face_nodes()[2] == 9;
  ldum = ldum && elem_def.get_number_of_face_nodes()[3] == 9;
  ldum = ldum && elem_def.get_number_of_face_nodes()[4] == 9;
  ldum = ldum && elem_def.get_number_of_face_nodes()[5] == 9;
  FAIL_IF_NOT(ldum);

#ifndef __clang_analyzer__
  // clang-analyze sees the throw and warns that 'ldum=false' is 'deadcode'.
  try {
    elem_def.get_side_type(6);
    ldum = false;
  } catch (...) {
  }

  try {
    elem_def.get_side_nodes(6);
    ldum = false;
  } catch (...) {
  }
#endif

  if (ldum) {
    ostringstream message;
    message << ename << " Element OK." << endl;
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "Error in " << ename << " Element." << endl;
    FAILMSG(message.str());
  }
  return ldum;
}
} // namespace rtt_mesh_element_test

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    runTest(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of TestElementDefinition.cc
//------------------------------------------------------------------------------------------------//
