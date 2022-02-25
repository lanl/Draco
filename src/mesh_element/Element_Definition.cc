//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh_element/Element_Definition.cc
 * \author John McGhee
 * \date   Fri Feb 25 10:03:18 2000
 * \brief  Provides some descriptive information for the standard mesh elements.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Element_Definition.hh"
#include <algorithm>
#include <iostream>

namespace rtt_mesh_element {

//------------------------------------------------------------------------------------------------//
Element_Definition::Element_Definition(Element_Type const &type_)
    : name(), type(type_), dimension(0), number_of_nodes(0), number_of_sides(0), elem_defs(),
      side_type(), side_nodes() {
  switch (type) {

  case NODE:
    construct_node();
    break;

  case BAR_2:
  case BAR_3:
    construct_bar();
    Ensure(invariant_satisfied());
    break;

  case TRI_3:
  case TRI_6:
    construct_tri();
    Ensure(invariant_satisfied());
    break;

  case QUAD_4:
  case QUAD_5:
  // Note: see Element_Definition.hh for the distinction between 6/6a and 6o
  // elements.
  case QUAD_6:
  case QUAD_6a:
  case QUAD_6o:
  case QUAD_7:
  case QUAD_8:
  case QUAD_9:
    construct_quad();
    Ensure(invariant_satisfied());
    break;

  case TETRA_4:
  case TETRA_10:
    construct_tetra();
    Ensure(invariant_satisfied());
    break;

  case PYRA_5:
  case PYRA_14:
    construct_pyra();
    Ensure(invariant_satisfied());
    break;

  case PENTA_6:
  case PENTA_15:
  case PENTA_18:
    construct_penta();
    Ensure(invariant_satisfied());
    break;

  case HEXA_8:
  case HEXA_20:
  case HEXA_27:
    construct_hexa();
    Ensure(invariant_satisfied());
    break;

  case POLYHEDRON:
    dimension = 3;
    break;

  case POLYGON:
    dimension = 2;
    break;

  default:
    Insist(false, "Unrecognized Element-Type Flag");
  }
}

//------------------------------------------------------------------------------------------------//
Element_Definition::Element_Definition(std::string const &name_, unsigned dimension_,
                                       unsigned number_of_nodes_, unsigned number_of_sides_,
                                       std::vector<Element_Definition> const &elem_defs_,
                                       std::vector<unsigned> const &side_type_,
                                       std::vector<std::vector<unsigned>> const &side_nodes_)
    : name(name_), type((dimension_ == 2) ? POLYGON : POLYHEDRON), dimension(dimension_),
      number_of_nodes(number_of_nodes_), number_of_sides(number_of_sides_), elem_defs(elem_defs_),
      side_type(side_type_), side_nodes(side_nodes_) {
  //----------------------------------------------------------------------------------------------//
  // Check input first, before any modifications
  Require(number_of_nodes_ > 0);

#ifdef REQUIRE_ON
  for (const auto &elem_def : elem_defs_)
    Require(elem_def.get_dimension() + 1 == dimension_);
#endif
  Require(side_type_.size() == number_of_sides_);
  for (unsigned i = 0; i < number_of_sides_; ++i) {
    Require(static_cast<unsigned>(side_type_[i]) < elem_defs_.size());
  }
  Require(side_nodes_.size() == number_of_sides_);
  for (unsigned i = 0; i < number_of_sides_; ++i) {
    Require(side_nodes_[i].size() == elem_defs_[side_type_[i]].get_number_of_nodes());
#ifdef REQUIRE_ON
    for (unsigned const side_node : side_nodes_[i])
      Require(side_node < number_of_nodes_);
#endif
  }

  // Only time this constructor should be called
  Ensure(get_type() == Element_Definition::POLYGON || get_type() == Element_Definition::POLYHEDRON);

  Ensure(get_name() == name_);
  Ensure(get_dimension() == dimension_);
  Ensure(get_number_of_nodes() == number_of_nodes_);
  Ensure(get_number_of_sides() == number_of_sides_);
  for (unsigned i = 0; i < number_of_sides; ++i) {
    Ensure(get_side_nodes(i) == side_nodes_[i]);
  }

  Ensure(invariant_satisfied());
}

//------------------------------------------------------------------------------------------------//
bool Element_Definition::invariant_satisfied() const {
  bool ldum = (name.empty() == false);

  if (type == NODE) {
    ldum = ldum && (dimension == 0);
    ldum = ldum && (number_of_nodes == 1);
    ldum = ldum && (number_of_sides == 0);
    ldum = ldum && (elem_defs.size() == 0);
  } else {
    ldum = ldum && (dimension > 0);
    ldum = ldum && (dimension < 4);
    ldum = ldum && (number_of_nodes > dimension);
    ldum = ldum && (number_of_sides <= number_of_nodes);
    ldum = ldum && (number_of_sides > dimension);
    ldum = ldum && (elem_defs.size() > 0);
  }

  ldum = ldum && (side_type.size() == number_of_sides);
  ldum = ldum && (side_nodes.size() == number_of_sides);

  for (const auto &elem_def : elem_defs)
    ldum = ldum && (elem_def.dimension == dimension - 1);

  for (size_t i = 0; i < side_nodes.size(); i++) {
    ldum = ldum && (side_nodes[i].size() > 0);
    ldum = ldum && (side_nodes[i].size() == elem_defs[side_type[i]].number_of_nodes);
    for (unsigned int side_node : side_nodes[i]) {
      ldum = ldum && (side_node < number_of_nodes);
    }
  }

  return ldum;
}

//------------------------------------------------------------------------------------------------//
void Element_Definition::construct_node() {
  name = "NODE";
  dimension = 0;
  number_of_sides = 0;
  number_of_nodes = 1;
}

//------------------------------------------------------------------------------------------------//
void Element_Definition::construct_bar() {
  dimension = 1;
  number_of_sides = 2;
  side_nodes = {{0}, {1}};
  switch (type) {
  case BAR_2:
    name = "BAR_2";
    number_of_nodes = 2;
    break;
  case BAR_3:
    name = "BAR_3";
    number_of_nodes = 3;
    break;
  default:
    Insist(false, "#2 Unrecognized Element-Type Flag");
  }
  elem_defs.emplace_back(NODE);
  side_type = std::vector<unsigned>(number_of_sides, 0);
}

//------------------------------------------------------------------------------------------------//
void Element_Definition::construct_tri() {
  dimension = 2;
  number_of_sides = 3;
  switch (type) {
  case TRI_3:
    name = "TRI_3";
    elem_defs.emplace_back(BAR_2);
    number_of_nodes = 3;
    side_nodes = {{0, 1}, {1, 2}, {2, 0}};
    break;
  case TRI_6:
    name = "TRI_6";
    elem_defs.emplace_back(BAR_3);
    number_of_nodes = 6;
    side_nodes = {{0, 1, 3}, {1, 2, 4}, {2, 0, 5}};
    break;
  default:
    Insist(false, "#3 Unrecognized Element-Type Flag");
  }
  side_type = std::vector<unsigned>(number_of_sides, 0);
}

//------------------------------------------------------------------------------------------------//
void Element_Definition::construct_quad() {
  dimension = 2;
  number_of_sides = 4;
  side_nodes = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
  switch (type) {
  case QUAD_4:
    name = "QUAD_4";
    number_of_nodes = 4;
    elem_defs.emplace_back(BAR_2);
    side_type = std::vector<unsigned>(number_of_sides, 0);
    break;

  case QUAD_5:
    name = "QUAD_5";
    number_of_nodes = 5;
    elem_defs.emplace_back(BAR_2);
    elem_defs.emplace_back(BAR_3);

    // Three BAR_2 sides
    side_type = std::vector<unsigned>(3, 0);

    // One BAR_3 side (add midpoint nodes)
    for (unsigned i = 3; i < number_of_sides; i++) {
      side_type.push_back(1);
      side_nodes[i].push_back(i - 3 + number_of_sides);
    }
    break;

  case QUAD_6:
  case QUAD_6a:
    name = "QUAD_6";
    number_of_nodes = 6;
    elem_defs.emplace_back(BAR_2);
    elem_defs.emplace_back(BAR_3);

    // Nodes on adjacent sides

    // Two BAR_2 sides
    side_type = std::vector<unsigned>(2, 0);

    // Two BAR_3 sides (add midpoint nodes)
    for (unsigned i = 2; i < number_of_sides; i++) {
      side_type.push_back(1);
      side_nodes[i].push_back(i - 2 + number_of_sides);
    }
    break;

  case QUAD_6o:
    name = "QUAD_6o";
    number_of_nodes = 6;
    elem_defs.emplace_back(BAR_2);
    elem_defs.emplace_back(BAR_3);

    // Nodes on opposite sides

    // First BAR_2 side
    side_type.push_back(0);

    // First BAR_3 side (add midpoint node)
    side_type.push_back(1);
    side_nodes[1].push_back(4);

    // Second BAR_2 side
    side_type.push_back(0);

    // Second BAR_3 side (add midpoint node)
    side_type.push_back(1);
    side_nodes[3].push_back(5);
    break;

  case QUAD_7:
    name = "QUAD_7";
    number_of_nodes = 7;
    elem_defs.emplace_back(BAR_2);
    elem_defs.emplace_back(BAR_3);

    // One BAR_2 side
    side_type = {0};

    // Three BAR_3 sides (add midpoint nodes)
    for (unsigned i = 1; i < number_of_sides; i++) {
      side_type.push_back(1);
      side_nodes[i].push_back(i - 1 + number_of_sides);
    }
    break;

  case QUAD_8:
    name = "QUAD_8";
    number_of_nodes = 8;
    elem_defs.emplace_back(BAR_3);
    side_type = std::vector<unsigned>(number_of_sides, 0);
    for (unsigned i = 0; i < 4; ++i)
      side_nodes[i].push_back(i + 4);
    break;

  case QUAD_9:
    name = "QUAD_9";
    number_of_nodes = 9;
    elem_defs.emplace_back(BAR_3);
    side_type = std::vector<unsigned>(number_of_sides, 0);
    for (unsigned i = 0; i < 4; ++i)
      side_nodes[i].push_back(i + 4);
    break;

  default:
    Insist(false, "#5 Unrecognized Element-Type Flag");
  }
}

//------------------------------------------------------------------------------------------------//
void Element_Definition::construct_tetra() {
  dimension = 3;
  number_of_sides = 4;
  side_nodes = {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};

  switch (type) {
  case TETRA_4:
    name = "TETRA_4";
    number_of_nodes = 4;
    elem_defs.emplace_back(TRI_3);
    break;
  case TETRA_10:
    name = "TETRA_10";
    number_of_nodes = 10;
    elem_defs.emplace_back(TRI_6);
    side_nodes[0].insert(side_nodes[0].end(), {6, 5, 4});
    side_nodes[1].insert(side_nodes[1].end(), {4, 8, 7});
    side_nodes[2].insert(side_nodes[2].end(), {5, 9, 8});
    side_nodes[3].insert(side_nodes[3].end(), {6, 7, 9});
    break;
  default:
    Insist(false, "#6 Unrecognized Element-Type Flag");
  }
  side_type = std::vector<unsigned>(number_of_sides, 0);
}

//------------------------------------------------------------------------------------------------//
void Element_Definition::construct_pyra() {
  dimension = 3;
  number_of_sides = 5;
  side_nodes = {{0, 3, 2, 1}, {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}};

  switch (type) {
  case PYRA_5:
    name = "PYRA_5";
    number_of_nodes = 5;
    elem_defs.emplace_back(QUAD_4);
    elem_defs.emplace_back(TRI_3);
    break;
  case PYRA_14:
    name = "PYRA_14";
    number_of_nodes = 14;
    elem_defs.emplace_back(QUAD_8);
    elem_defs.emplace_back(TRI_6);
    side_nodes[0].insert(side_nodes[0].end(), {8, 7, 6, 5});
    side_nodes[1].insert(side_nodes[1].end(), {5, 10, 9});
    side_nodes[2].insert(side_nodes[2].end(), {6, 11, 10});
    side_nodes[3].insert(side_nodes[3].end(), {7, 12, 11});
    side_nodes[4].insert(side_nodes[4].end(), {8, 9, 12});
    break;
  default:
    Insist(false, "#7 Unrecognized Element-Type Flag");
  }

  side_type.push_back(0);
  for (size_t i = 1; i < number_of_sides; i++)
    side_type.push_back(1);
}

//------------------------------------------------------------------------------------------------//
void Element_Definition::construct_penta() {
  dimension = 3;
  number_of_sides = 5;
  side_nodes = {{0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}, {0, 2, 1}, {3, 4, 5}};

  switch (type) {
  case PENTA_6:
    name = "PENTA_6";
    number_of_nodes = 6;
    elem_defs.emplace_back(QUAD_4);
    elem_defs.emplace_back(TRI_3);
    break;

  case PENTA_15:
  case PENTA_18:
    side_nodes[0].insert(side_nodes[0].end(), {6, 10, 12, 9});
    side_nodes[1].insert(side_nodes[1].end(), {7, 11, 13, 10});
    side_nodes[2].insert(side_nodes[2].end(), {8, 9, 14, 11});
    side_nodes[3].insert(side_nodes[3].end(), {8, 7, 6});
    side_nodes[4].insert(side_nodes[4].end(), {12, 13, 14});
    switch (type) {
    case PENTA_15:
      name = "PENTA_15";
      number_of_nodes = 15;
      elem_defs.emplace_back(QUAD_8);
      elem_defs.emplace_back(TRI_6);
      break;
    case PENTA_18:
      name = "PENTA_18";
      number_of_nodes = 18;
      elem_defs.emplace_back(QUAD_9);
      elem_defs.emplace_back(TRI_6);
      side_nodes[0].push_back(15);
      side_nodes[1].push_back(16);
      side_nodes[2].push_back(17);
      break;
    default:
      Insist(false, "#8 Unrecognized Element-Type Flag");
    }
    break;

  default:
    Insist(false, "#9 Unrecognized Element-Type Flag");
  }

  side_type = {0, 0, 0};
  for (size_t i = 3; i < number_of_sides; i++)
    side_type.push_back(1);
}

//------------------------------------------------------------------------------------------------//
void Element_Definition::construct_hexa() {
  dimension = 3;
  number_of_sides = 6;
  side_nodes = {{0, 3, 2, 1}, {0, 4, 7, 3}, {2, 3, 7, 6}, {1, 2, 6, 5}, {0, 1, 5, 4}, {4, 5, 6, 7}};

  switch (type) {
  case HEXA_8:
    name = "HEXA_8";
    number_of_nodes = 8;
    elem_defs.emplace_back(QUAD_4);
    break;
  case HEXA_20:
  case HEXA_27:
    side_nodes[0].insert(side_nodes[0].end(), {11, 10, 9, 8});
    side_nodes[1].insert(side_nodes[1].end(), {12, 19, 15, 11});
    side_nodes[2].insert(side_nodes[2].end(), {10, 15, 18, 14});
    side_nodes[3].insert(side_nodes[3].end(), {9, 14, 17, 13});
    side_nodes[4].insert(side_nodes[4].end(), {8, 13, 16, 12});
    side_nodes[5].insert(side_nodes[5].end(), {16, 17, 18, 19});
    switch (type) {
    case HEXA_20:
      name = "HEXA_20";
      number_of_nodes = 20;
      elem_defs.emplace_back(QUAD_8);
      break;
    case HEXA_27:
      name = "HEXA_27";
      number_of_nodes = 27;
      elem_defs.emplace_back(QUAD_9);
      for (unsigned i = 0; i < number_of_sides; i++)
        side_nodes[i].push_back(i + 20);
      break;
    default:
      Insist(false, "#10 Unrecognized Element-Type Flag");
    }
    break;
  default:
    Insist(false, "#11 Unrecognized Element-Type Flag");
  }
  side_type = std::vector<unsigned>(number_of_sides, 0);
}

//------------------------------------------------------------------------------------------------//
std::ostream &Element_Definition::print(std::ostream &os_out) const {
  os_out << "Element Type   : " << get_type() << "\n"
         << "Element Name   : " << get_name() << "\n"
         << "Number of Nodes: " << get_number_of_nodes() << "\n"
         << "Dimension      : " << get_dimension() << "\n"
         << "Number of Sides: " << get_number_of_sides() << "\n"
         << "\n";
  if (get_number_of_sides() != 0) {
    os_out << "Side Types     : ";
    for (unsigned j = 0; j < get_number_of_sides(); j++)
      os_out << get_side_type(j).get_name() << " ";
    os_out << "\n";

    std::vector<unsigned> tmp;
    os_out << "Side Nodes     : "
           << "\n";
    for (unsigned j = 0; j < get_number_of_sides(); j++) {
      tmp = get_side_nodes(j);
      os_out << "  side# " << j << " -    ";
      for (unsigned k : tmp)
        os_out << k << " ";
      os_out << "\n";
    }
  }
  {
    std::vector<unsigned> num_face_nodes = get_number_of_face_nodes();
    std::vector<std::vector<unsigned>> face_nodes = get_face_nodes();
    os_out << "Face Nodes: " << num_face_nodes.size() << "\n";
    for (size_t j = 0; j < num_face_nodes.size(); ++j) {
      os_out << "  Face " << j << ": " << num_face_nodes[j] << " nodes : ";
      for (size_t k = 0; k < num_face_nodes[j]; ++k) {
        os_out << face_nodes[j][k] << " ";
      }
      os_out << "\n";
    }
  }
  os_out << std::endl;
  return os_out;
}

} // end namespace rtt_mesh_element

//------------------------------------------------------------------------------------------------//
// end of Element_Definition.cc
//------------------------------------------------------------------------------------------------//
