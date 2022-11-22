//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/CellDefs.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/CellDefs class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "CellDefs.hh"

namespace rtt_RTT_Format_Reader {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Parses the cell_defs (cell definitions) data block from the mesh file via calls to private
 *        member functions.
 * \param meshfile Mesh file name.
 */
void CellDefs::readCellDefs(ifstream &meshfile) {
  readKeyword(meshfile);
  readDefs(meshfile);
  readEndKeyword(meshfile);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the cell_defs block (cell definitions) keyword.
 * \param meshfile Mesh file name.
 */
void CellDefs::readKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "cell_defs", "Invalid mesh file: cell_defs block missing");
  std::getline(meshfile, dummyString);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the cell_defs (cell definitions) block data.
 * \param meshfile Mesh file name.
 */
void CellDefs::readDefs(ifstream &meshfile) {
  int cellDefNum;
  string dummyString;

  for (size_t i = 0; i < dims.get_ncell_defs(); ++i) {
    meshfile >> cellDefNum >> dummyString;
    Insist(static_cast<size_t>(cellDefNum) == i + 1, "Invalid mesh file: cell def out of order");
    // Ignore plurals in cell definitions
    if (dummyString[dummyString.size() - 1] == 's')
      dummyString.resize(dummyString.size() - 1);
    Check(i < defs.size());
    defs[i] = std::make_shared<CellDef>(*this, dummyString);
    std::getline(meshfile, dummyString);
    defs[i]->readDef(meshfile);
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Reads and validates the end_cell_defs block keyword.
 * \param meshfile Mesh file name.
 */
void CellDefs::readEndKeyword(ifstream &meshfile) {
  string dummyString;

  meshfile >> dummyString;
  Insist(dummyString == "end_cell_defs", "Invalid mesh file: cell_defs block missing end");
  std::getline(meshfile, dummyString); // read and discard blank line.
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Changes the cell definitions specified in the RTT_Format file to an alternative
 *        coordinate-system independent cell definition (e.g., CYGNUS).
 * \param cell_side_types New side types for each of the existing cell definitions.
 * \param cell_ordered_sides New ordered sides for each of the existing cell definitions.
 */
void CellDefs::redefineCellDefs(vector_vector_uint const &cell_side_types,
                                std::vector<vector_vector_uint> const &cell_ordered_sides) {
  Insist(cell_side_types.size() == dims.get_ncell_defs(),
         "Error in supplied cell redefinition side types data.");
  Insist(cell_ordered_sides.size() == dims.get_ncell_defs(),
         "Error in supplied cell redefinition ordered side data.");

  redefined = true;

  for (size_t cd = 0; cd < dims.get_ncell_defs(); cd++) {
    Check(cd < defs.size());
    Check(cd < cell_side_types.size());
    Check(cd < cell_ordered_sides.size());
    defs[cd]->redefineCellDef(cell_side_types[cd], cell_ordered_sides[cd], dims.get_ndim());
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Used by the CellDefs class objects to parse the number of nodes and sides per cell, the
 *        side type indices, and the nodes for each side.
 * \param meshfile Mesh file name.
 */
void CellDef::readDef(ifstream &meshfile) {
  string dummyString;

  meshfile >> nnodes >> nsides;
  side_types.resize(nsides);
  sides.resize(nsides);
  ordered_sides.resize(nsides);
  std::getline(meshfile, dummyString);

  for (unsigned i = 0; i < nsides; ++i) {
    Check(i < side_types.size());
    meshfile >> side_types[i];
    --side_types[i];
  }
  if (nsides > 0)
    std::getline(meshfile, dummyString);

  // note that this implementation does not preserve the "right hand rule" of the cell definitions
  // due to the use of a set container (which is sorted).  It is slicker than snail snot when it
  // comes time to implement the connectivity, however. The ordered_sides vector was added to allow
  // the original ordered data to be retained.
  int side;
  for (unsigned i = 0; i < nsides; ++i) {
    Check(i < side_types.size());
    size_t numb_nodes = cellDefs.get_cell_def(side_types[i]).get_nnodes();
    Check(i < ordered_sides.size());
    ordered_sides[i].resize(numb_nodes);
    Check(i < sides.size());
    for (size_t j = 0; j < numb_nodes; ++j) {
      meshfile >> side;
      --side;
      sides[i].push_back(side);
      Check(j < ordered_sides[i].size());
      ordered_sides[i][j] = side;
    }
    if (sides[i].size() > 0)
      std::getline(meshfile, dummyString);
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Changes the cell definitions specified in the RTT_Format file to an alternative
 *        coordinate-system independent cell definition (e.g., CYGNUS).
 * \param new_side_types New cell side types.
 * \param new_ordered_sides New cell ordered sides.
 * \param ndim Topological dimension of the cells in the mesh.
 */
void CellDef::redefineCellDef(vector_uint const &new_side_types,
                              vector_vector_uint const &new_ordered_sides, size_t const ndim) {
  Insist(new_side_types.size() == nsides, "New side types input error");
  Insist(new_ordered_sides.size() == nsides, "New ordered sides input error");

  Require(ndim <= 3);

  node_map.resize(nnodes);

  if (name == "point") {
    node_map[0] = 0;
  } else if (name == "line" || name == "bar2") {
    Check(new_ordered_sides[0][0] < INT_MAX);
    node_map[ordered_sides[0][0]] = static_cast<int>(new_ordered_sides[0][0]);
    Check(new_ordered_sides[1][0] < INT_MAX);
    node_map[ordered_sides[1][0]] = static_cast<int>(new_ordered_sides[1][0]);
  } else if (name == "line_qdr" || name == "bar3") {
    Check(new_ordered_sides[0][0] < INT_MAX);
    node_map[ordered_sides[0][0]] = static_cast<int>(new_ordered_sides[0][0]);
    Check(new_ordered_sides[1][0] < INT_MAX);
    node_map[ordered_sides[1][0]] = static_cast<int>(new_ordered_sides[1][0]);
    // kgb (060307): I'm guessing BTA never thought about how to map internal nodes, so it's not
    // clear how to proceed here. My best guess is we find the node not in each map.
    unsigned old_node;
    for (old_node = 0; old_node < 3; ++old_node) {
      if (ordered_sides[0][0] != old_node && ordered_sides[1][0] != old_node) {
        break;
      }
    }
    unsigned new_node;
    for (new_node = 0; new_node < 3; ++new_node) {
      if (new_ordered_sides[0][0] != new_node && new_ordered_sides[1][0] != new_node) {
        break;
      }
    }
    node_map[old_node] = new_node;
  } else if (name == "triangle" || name == "tri3" || name == "quad" || name == "quad4") {
    // Arbitrarily assign the first node in the old and the new cell definitions to be the
    // same. This assumption is necessary because the cell definitions do not assume a specific
    // orientation relative to any coordinate system.  The transformed cell may be rotated about
    // it's outward normal relative to the input cell definition.
    node_map[0] = 0;
    // The right hand rule has to apply, so only the ordering of the nodes (edges) can change for a
    // two-dimensional cell.
    size_t old_node = 0;
    size_t new_node = 0;
    for (size_t n = 0; n < static_cast<size_t>(nnodes - 1); n++) {
      // Find the new side that starts with this node.
      size_t new_side = 0;
      while (static_cast<size_t>(new_ordered_sides[new_side][0]) != new_node) {
        ++new_side;
        Insist(new_side < nsides, "Edge error for new two dimensional cell definition.");
      }
      new_node = new_ordered_sides[new_side][1];
      // Find the old side that starts with this node.
      size_t old_side = 0;
      while (ordered_sides[old_side][0] != old_node) {
        ++old_side;
        Insist(old_side < nsides, "Edge error for old two dimensional cell definition.");
      }
      old_node = ordered_sides[old_side][1];
      Check(new_node < INT_MAX);
      node_map[old_node] = static_cast<int>(new_node);
    }
  } else if (name == "triangle_qdr" || name == "tri6") {
    // Arbitrarily assign the first node in the old and the new cell definitions to be the
    // same. This assumption is necessary because the cell definitions do not assume a specific
    // orientation relative to any coordinate system. The transformed cell may be rotated about it's
    // outward normal relative to the input cell definition.
    node_map[0] = 0;
    // The right hand rule has to apply, so only the ordering of the nodes (edges) can change for a
    // two-dimensional cell.
    size_t old_node = 0;
    size_t new_node = 0;
    for (size_t n = 0; n < nnodes - nsides; n++) {
      // Find the new side that starts with this node.
      size_t new_side = 0;
      while (new_ordered_sides[new_side][0] != new_node) {
        ++new_side;
        Insist(new_side < nsides, "Edge error for new two dimensional cell definition.");
      }
      new_node = new_ordered_sides[new_side][1];
      // Find the old side that starts with this node.
      size_t old_side = 0;
      while (ordered_sides[old_side][0] != old_node) {
        ++old_side;
        Insist(old_side < nsides, "Edge error for old two dimensional cell definition.");
      }
      old_node = ordered_sides[old_side][1];
      Check(new_node < INT_MAX);
      node_map[old_node] = static_cast<int>(new_node);
      Check(new_ordered_sides[new_side][2] < INT_MAX);
      int new_mid_node = static_cast<int>(new_ordered_sides[new_side][2]);
      Check(new_ordered_sides[old_side][2] < INT_MAX);
      int old_mid_node = static_cast<int>(ordered_sides[old_side][2]);
      node_map[old_mid_node] = new_mid_node;
    }
  } else if (name == "tetrahedron") {
    // Arbitrarily assign the first node in the old and the new cell definitions to be the
    // same. This assumption is necessary because the cell definitions do not assume a specific
    // orientation relative to any coordinate system. The transformed cell may be rotated about the
    // outward normal of the opposite face relative to the input cell definition.
    node_map[0] = 0;
    // Find the one side definition that does not contain the first node.
    size_t new_side = 0;
    while (std::count(new_ordered_sides[new_side].begin(), new_ordered_sides[new_side].end(),
                      node_map[0]) > 0) {
      ++new_side;
      Insist(new_side < nsides, "Side error for new tetrahedron cell definition.");
    }
    // Find the one side definition that does not contain the first node.
    size_t old_side = 0;
    while (std::count(ordered_sides[old_side].begin(), ordered_sides[old_side].end(), node_map[0]) >
           0) {
      ++old_side;
      Insist(old_side < nsides, "Side error for old tetrahedron cell definition.");
    }
    // Now just apply the right-hand rule.
    for (size_t n = 0; n < ordered_sides[old_side].size(); n++) {
      Check(new_ordered_sides[new_side][n] < INT_MAX);
      node_map[ordered_sides[old_side][n]] = static_cast<int>(new_ordered_sides[new_side][n]);
    }
  } else if (name == "quad_pyr") {
    // Find the side that is the quad. The transformed cell may be rotated about the outward normal
    // of this face relative to the input cell definition.
    size_t new_side = 0;
    while (new_ordered_sides[new_side].size() != 4) {
      ++new_side;
      Insist(new_side < nsides, "Quad side error for new quad pyramid cell definition.");
    }
    size_t old_side = 0;
    while (ordered_sides[old_side].size() != 4) {
      ++old_side;
      Insist(old_side < nsides, "Quad side error for old quad pyramid cell definition.");
    }
    // Find the single node that is not included in the quad side definition and assign this to the
    // node map.
    int new_node = 0;
    int old_node = 0;
    for (size_t n = 0; n < nnodes; n++) {
      if (std::count(new_ordered_sides[new_side].begin(), new_ordered_sides[new_side].end(), n) ==
          0) {
        Check(n < INT_MAX);
        new_node = static_cast<int>(n);
      }
      if (std::count(ordered_sides[old_side].begin(), ordered_sides[old_side].end(), n) == 0) {
        Check(n < INT_MAX);
        old_node = static_cast<int>(n);
      }
    }
    node_map[old_node] = new_node;
    // Now just apply the right-hand rule to the quad side.
    for (size_t n = 0; n < ordered_sides[old_side].size(); n++) {
      Check(new_ordered_sides[new_side][n] < INT_MAX);
      node_map[ordered_sides[old_side][n]] = static_cast<int>(new_ordered_sides[new_side][n]);
    }
  } else if (name == "tri_prism") {
    // Find the one quad side definition that does not contain the first node.  The transformed cell
    // may be rotated about the outward normal of this face relative to the input cell definition.
    size_t new_quad = 0;
    while (new_ordered_sides[new_quad].size() != 4 ||
           std::count(new_ordered_sides[new_quad].begin(), new_ordered_sides[new_quad].end(), 0U) >
               0) {
      ++new_quad;
      Insist(new_quad < nsides, "Quad side error for new tri-prism cell definition.");
    }
    // Find the one quad side definition that does not contain the first node.
    size_t old_quad = 0;
    while (ordered_sides[old_quad].size() != 4 ||
           std::count(ordered_sides[old_quad].begin(), ordered_sides[old_quad].end(), 0U) > 0) {
      ++old_quad;
      Insist(old_quad < nsides, "Quad side error for old tri-prism cell definition.");
    }
    // Apply the right-hand rule to this quad.
    for (size_t n = 0; n < ordered_sides[old_quad].size(); n++) {
      Check(new_ordered_sides[new_quad][n] < INT_MAX);
      node_map[ordered_sides[old_quad][n]] = static_cast<int>(new_ordered_sides[new_quad][n]);
    }
    // Equate the two remaining triangle nodes. Find the first node.
    size_t old_tri = 0;
    while (ordered_sides[old_tri].size() != 3 ||
           std::count(ordered_sides[old_tri].begin(), ordered_sides[old_tri].end(), 0U) > 0) {
      ++old_tri;
      Insist(old_tri < nsides, "Triangle side error for old tri-prism cell definition.");
    }
    // Find out which of the two triangles this is by identifying one of the nodes common to both
    // the triangle and the previous quad.
    size_t old_node = 0;
    while (std::count(ordered_sides[old_quad].begin(), ordered_sides[old_quad].end(),
                      ordered_sides[old_tri][old_node]) == 0) {
      ++old_node;
      Insist(old_node < ordered_sides[old_tri].size(),
             "Node error for old tri-prism cell definition.");
    }
    size_t new_tri = 0;
    while (new_ordered_sides[new_tri].size() != 3 ||
           std::count(new_ordered_sides[new_tri].begin(), new_ordered_sides[new_tri].end(),
                      node_map[ordered_sides[old_tri][old_node]]) == 0) {
      ++new_tri;
      Insist(new_tri < nsides, "Triangle side error for new tri-prism cell definition.");
    }
    --old_node;
    size_t new_node = 0;
    while (std::count(new_ordered_sides[new_quad].begin(), new_ordered_sides[new_quad].end(),
                      new_ordered_sides[new_tri][new_node]) != 0) {
      ++new_node;
      Insist(new_node < new_ordered_sides[new_tri].size(),
             "Node error for new tri-prism cell definition.");
    }
    Check(new_ordered_sides[new_tri][new_node] < INT_MAX);
    node_map[ordered_sides[old_tri][old_node]] =
        static_cast<int>(new_ordered_sides[new_tri][new_node]);
    // The node that is neither in the previous quad or triangle is all that is left.
    for (size_t n = 0; n < nnodes; n++) {
      if (std::count(new_ordered_sides[new_quad].begin(), new_ordered_sides[new_quad].end(), n) ==
              0 &&
          std::count(new_ordered_sides[new_tri].begin(), new_ordered_sides[new_tri].end(), n) == 0)
        new_node = n;

      if (std::count(ordered_sides[old_quad].begin(), ordered_sides[old_quad].end(), n) == 0 &&
          std::count(ordered_sides[old_tri].begin(), ordered_sides[old_tri].end(), n) == 0)
        old_node = n;
    }
    Check(new_node < INT_MAX);
    node_map[old_node] = static_cast<int>(new_node);
  } else if (name == "hexahedron") {
    // Arbitrarily assign the first quad and the associated nodes in the old and the new cell
    // definitions to be the same. This assumption is necessary because the cell definitions do not
    // assume a specific orientation relative to any coordinate system. The transformed cell may be
    // rotated about it's coordinate system relative to the input cell definition.
    int quad = 0;
    vector_int new_node_count(nnodes, 0);
    vector_int old_node_count(nnodes, 0);
    for (size_t n = 0; n < ordered_sides[quad].size(); n++) {
      unsigned new_node = new_ordered_sides[quad][n];
      size_t old_node = ordered_sides[quad][n];
      node_map[old_node] = new_node;
      for (size_t s = 0; s < nsides; s++) {
        if (std::count(new_ordered_sides[s].begin(), new_ordered_sides[s].end(), new_node) > 0)
          for (uint64_t side : new_ordered_sides[s])
            ++new_node_count[side];

        if (std::count(ordered_sides[s].begin(), ordered_sides[s].end(), old_node) > 0)
          for (uint64_t side : ordered_sides[s])
            ++old_node_count[side];
      }
      // The node located diagonally across the hexahedron relative to the first node will have a
      // count of zero from the previous loop.
      for (size_t c = 0; c < nnodes; c++) {
        if (new_node_count[c] == 0) {
          Check(c < INT_MAX);
          new_node = static_cast<int>(c);
        }
        if (old_node_count[c] == 0)
          old_node = c;
      }
      node_map[old_node] = new_node;
      std::fill(new_node_count.begin(), new_node_count.end(), 0);
      std::fill(old_node_count.begin(), old_node_count.end(), 0);
    }
  } else {
    if (ndim == 2) // POLYGON
    {
      // Arbitrarily assign the first node in the old and the new cell definitions to be the
      // same. This assumption is necessary because the cell definitions do not assume a specific
      // orientation relative to any coordinate system. The transformed cell may be rotated about
      // it's outward normal relative to the input cell definition.
      node_map[0] = 0;
      // The right hand rule has to apply, so only the ordering of the nodes (edges) can change for
      // a two-dimensional cell.
      size_t old_node = 0;
      size_t new_node = 0;
      for (size_t n = 0; n < nnodes - 1; n++) {
        // Find the new side that starts with this node.
        size_t new_side = 0;
        while (new_ordered_sides[new_side][0] != new_node) {
          ++new_side;
          Insist(new_side < nsides, "Edge error for new two dimensional cell definition.");
        }
        new_node = new_ordered_sides[new_side][1];
        // Find the old side that starts with this node.
        size_t old_side = 0;
        while (ordered_sides[old_side][0] != old_node) {
          ++old_side;
          Insist(old_side < nsides, "Edge error for old two dimensional cell definition.");
        }
        old_node = ordered_sides[old_side][1];
        Check(new_node < INT_MAX);
        node_map[old_node] = static_cast<int>(new_node);
      }
    } else if (ndim == 3) // POLYHEDRON
    {
      for (unsigned i = 0; i < nnodes; ++i)
        node_map[i] = i;
    }
  }

  // Assign the new side types, sides, and ordered sides to this cell definition.
  side_types = new_side_types;
  ordered_sides = new_ordered_sides;
  for (size_t i = 0; i < ordered_sides.size(); ++i) {
    sides[i].erase(sides[i].begin(), sides[i].end());
    for (size_t n = 0; n < ordered_sides[i].size(); ++n)
      sides[i].push_back(ordered_sides[i][n]);
  }
}

} // end namespace rtt_RTT_Format_Reader

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/CellDefs.cc
//------------------------------------------------------------------------------------------------//
