//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   meshReaders/Hex_Mesh_Reader.cc
 * \author John McGhee
 * \date   Tue Mar  7 08:38:04 2000
 * \brief  Implements a CIC-19 Hex Mesh Format mesh reader.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Hex_Mesh_Reader.hh"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

namespace rtt_meshReaders {

using rtt_mesh_element::Element_Definition;

//------------------------------------------------------------------------------------------------//
Hex_Mesh_Reader::Hex_Mesh_Reader(std::string const &filename)
    : meshfile_name(filename), version("unknown"), npoints(0), ncells(0), nvrtx(0), nvrpf(0),
      ndim(0), nvb_faces(0), nrb_faces(0), nmat(0), point_coords(), ipar(), imat_index(),
      irgn_vb_index(), ipar_vb(), ipar_rb(), node_sets() {

  // Open the mesh file for read.
  std::ifstream meshfile(filename.c_str(), std::ios::in);
  if (!meshfile)
    Insist(false, "Could not open mesh-file!");

  // Check to make sure this is a CIC-19 Hex format file.
  std::string chdum;
  meshfile >> chdum;
  if (chdum != keyword())
    Insist(false, "Not a CIC-19 Hex Mesh File!");

  // Read in the dimensions of the problem.
  meshfile >> npoints >> ncells >> nvrtx >> nvrpf >> ndim >> nvb_faces >> nrb_faces >> nmat;
  Insist(check_dims(), "Error in Mesh Dimension data!");

  // Read the point coordinates data.
  point_coords.resize(npoints);
  for (unsigned i = 0; i < npoints; i++) {
    point_coords[i].resize(ndim);
    if (ndim == 1)
      meshfile >> point_coords[i][0];
    else if (ndim == 2)
      meshfile >> point_coords[i][0] >> point_coords[i][1];
    else if (ndim == 3)
      meshfile >> point_coords[i][0] >> point_coords[i][1] >> point_coords[i][2];
    else
      Insist(false, "Dimension index out of range!");
  }

  // Read in the mesh connectivity.
  ipar.resize(ncells);
  for (unsigned i = 0; i < ncells; i++) {
    ipar[i].resize(nvrtx);
    if (ndim == 1)
      meshfile >> ipar[i][0] >> ipar[i][1];
    else if (ndim == 2)
      meshfile >> ipar[i][0] >> ipar[i][1] >> ipar[i][2] >> ipar[i][3];
    else if (ndim == 3)
      meshfile >> ipar[i][0] >> ipar[i][1] >> ipar[i][2] >> ipar[i][3] >> ipar[i][4] >>
          ipar[i][5] >> ipar[i][6] >> ipar[i][7];
    else
      Insist(false, "Dimension index out of range!");
    for (size_t j = 0; j < nvrtx; j++)
      ipar[i][j] = ipar[i][j] - 1;
  }

  // Read in the mesh interior-region data.
  imat_index.resize(ncells);
  for (unsigned i = 0; i < ncells; i++)
    meshfile >> imat_index[i];

  // Read in the mesh vacuum boundary data.
  ipar_vb.resize(nvb_faces);
  irgn_vb_index.resize(nvb_faces);
  for (size_t i = 0; i < nvb_faces; i++) {
    ipar_vb[i].resize(nvrpf);
    if (ndim == 1)
      meshfile >> ipar_vb[i][0] >> irgn_vb_index[i];
    else if (ndim == 2)
      meshfile >> ipar_vb[i][0] >> ipar_vb[i][1] >> irgn_vb_index[i];
    else if (ndim == 3)
      meshfile >> ipar_vb[i][0] >> ipar_vb[i][1] >> ipar_vb[i][2] >> ipar_vb[i][3] >>
          irgn_vb_index[i];
    else
      Insist(false, "Dimension index out of range!");
    for (size_t j = 0; j < nvrpf; j++)
      ipar_vb[i][j] = ipar_vb[i][j] - 1;
  }

  // Read in the mesh reflective boundary data.
  ipar_rb.resize(nrb_faces);
  for (size_t i = 0; i < nrb_faces; i++) {
    ipar_rb[i].resize(nvrpf);
    if (ndim == 1)
      meshfile >> ipar_rb[i][0];
    else if (ndim == 2)
      meshfile >> ipar_rb[i][0] >> ipar_rb[i][1];
    else if (ndim == 3)
      meshfile >> ipar_rb[i][0] >> ipar_rb[i][1] >> ipar_rb[i][2] >> ipar_rb[i][3];
    else
      Insist(false, "Dimension index out of range!");
    for (size_t j = 0; j < nvrpf; j++)
      ipar_rb[i][j] = ipar_rb[i][j] - 1;
  }

  // Load a default node set which contains all the nodes in the mesh.

  std::set<unsigned> stmp;
  for (unsigned i = 0; i < npoints; ++i)
    stmp.insert(i);
  using resultT = std::map<std::string, std::set<unsigned>>;
  node_sets.insert(resultT::value_type("Interior", stmp));
}

//------------------------------------------------------------------------------------------------//
/*!
 * Returns all the ndim-dimensional interior elements as well as the (ndim-1) dimensional vacuum and
 * reflective boundary face elements
 */
std::vector<std::vector<unsigned>> Hex_Mesh_Reader::get_element_nodes() const {
  // Collate the interior, vacuum, and reflective mesh elements into one vector. Note that the order
  // is important as we will rely on it later to output the element set data.

  // Alternatively, the private data of the class could be changed so that the work done here is
  // done in the constructor.  This would be more efficient if this is going to be used
  // repetitively.
  std::vector<std::vector<unsigned>> result;
  for (unsigned i = 0; i < ncells; i++)
    result.push_back(ipar[i]);
  for (size_t i = 0; i < nvb_faces; i++)
    result.push_back(ipar_vb[i]);
  for (size_t i = 0; i < nrb_faces; i++)
    result.push_back(ipar_rb[i]);
  return result;
}

//------------------------------------------------------------------------------------------------//
/*!
 * Returns an element type for each element in the mesh. Will always be one of
 * - rtt_mesh_element::Element_Definition::NODE,
 * - rtt_mesh_element::Element_Definition::BAR_2,
 * - rtt_mesh_element::Element_Definition::QUAD_4, or
 * - rtt_mesh_element::Element_Definition::HEXA_8.
 */
std::vector<Element_Definition::Element_Type> Hex_Mesh_Reader::get_element_types() const {
  Element_Definition::Element_Type d1, d2;
  switch (ndim) {
  case (1):
    d1 = Element_Definition::BAR_2;
    d2 = Element_Definition::NODE;
    break;
  case (2):
    d1 = Element_Definition::QUAD_4;
    d2 = Element_Definition::BAR_2;
    break;
  case (3):
    d1 = Element_Definition::HEXA_8;
    d2 = Element_Definition::QUAD_4;
    break;
  default:
    Insist(false, "Dimension index out of range!");
  }
  std::vector<Element_Definition::Element_Type> tmp;
  for (unsigned i = 0; i < ncells; i++)
    tmp.push_back(d1);
  for (size_t i = 0; i < nvb_faces + nrb_faces; i++)
    tmp.push_back(d2);
  return tmp;
}

//------------------------------------------------------------------------------------------------//
/*!
 * Returns the unique element types defined in the mesh. Will always be one of
 * - rtt_mesh_element::Element_Definition::NODE,
 * - rtt_mesh_element::Element_Definition::BAR_2,
 * - rtt_mesh_element::Element_Definition::QUAD_4, or
 * - rtt_mesh_element::Element_Definition::HEXA_8.
 */
std::vector<Element_Definition::Element_Type> Hex_Mesh_Reader::get_unique_element_types() const {
  std::vector<Element_Definition::Element_Type> tmp;
  tmp.push_back(Element_Definition::NODE);
  tmp.push_back(Element_Definition::BAR_2);
  switch (ndim) {
  case (1):
    break;
  case (2):
    tmp.push_back(Element_Definition::QUAD_4);
    break;
  case (3):
    tmp.push_back(Element_Definition::QUAD_4);
    tmp.push_back(Element_Definition::HEXA_8);
    break;
  default:
    Insist(false, "Dimension index out of range!");
  }
  return tmp;
}

//------------------------------------------------------------------------------------------------//
/*!
 * There is no provision for naming element sets in the Hex format. The following default names are
 * provided for the sets found on the mesh file:
 * - "Interior" -- All the ndim-dimensional cells in the problem.
 * - "Interior_Region_x" - Interior cells with integer flag "x".
 * - "Vacumm_Boundary" -- All the (ndim-1)dimensional vacuum boundary faces.
 * - "Vacuum_Boundary_Region_x" -- Vacuum boundary faces with flag integer "x"
 * - "Reflective_Boundary" -- All the (ndim-1) dimensional reflective boundary faces.
 */
std::map<std::string, std::set<unsigned>> Hex_Mesh_Reader::get_element_sets() const {
  // Alternatively, the private data of the class could be changed so that the work done here is
  // done in the constructor. This would be more efficient if this is going to be used repetively.
  using resultT = std::map<std::string, std::set<unsigned>>;
  resultT result;
  std::vector<int> tmp;
  std::set<int> rgn_index;
  std::set<unsigned> stmp;

  // Create a set that flags all interior cells.
  stmp.clear();
  for (unsigned i = 0; i < ncells; i++)
    stmp.insert(i);
  result.insert(resultT::value_type("Interior", stmp));

  // Create sets for all the interior mesh sub-regions.  This loops over the whole mesh
  // number_of_mesh_regions times. Could be made to do it more efficiently in one loop? Note that
  // this depends on the elements being stored in a specific order.
  rgn_index = std::set<int>(imat_index.begin(), imat_index.end());
  for (auto const &iregion : rgn_index) {
    std::ostringstream os_chdum("");
    // os_chdum << "Interior_Region_" << *i;
    os_chdum << "Interior_Region_" << iregion;
    stmp.clear();
    for (unsigned j = 0; j < ncells; j++)
      if (imat_index[j] == iregion)
        stmp.insert(j);
    result.insert(resultT::value_type(os_chdum.str(), stmp));
  }

  if (nvb_faces > 0) {
    // Create a vacuum boundary set. Note that this depends on the elements being stored in a
    // specific order.
    stmp.clear();
    for (unsigned i = ncells; i < ncells + nvb_faces; i++)
      stmp.insert(i);
    result.insert(resultT::value_type("Vacuum_Boundary", stmp));

    // Create sets for all the vacuum boundary regions.  This loops over the whole mesh
    // number_of_vb_regions times. Could be made to do it more efficiently in one loop? Note that
    // this depends on the elements being stored in a specific order.
    rgn_index = std::set<int>(irgn_vb_index.begin(), irgn_vb_index.end());
    for (auto const &iregion : rgn_index) {
      std::ostringstream os_chdum("");
      os_chdum << "Vacuum_Boundary_Region_" << iregion;
      stmp.clear();
      Check(nvb_faces < UINT_MAX);
      for (size_t j = 0; j < nvb_faces; j++) {
        if (irgn_vb_index[j] == iregion)
          stmp.insert(static_cast<unsigned>(j) + ncells);
      }
      if (stmp.size() != 0)
        result.insert(resultT::value_type(os_chdum.str(), stmp));
    }
  }

  // Create a reflective boundary set. Note that this depends on the elements being stored in a
  // specific order.
  if (nrb_faces > 0) {
    stmp.clear();
    for (unsigned i = ncells + nvb_faces; i < ncells + nvb_faces + nrb_faces; i++)
      stmp.insert(i);
    result.insert(resultT::value_type("Reflective_Boundary", stmp));
  }

  return result;
}

//------------------------------------------------------------------------------------------------//
//! Checks the internal consistancy of the Hex_Mesh_Reader private data.
bool Hex_Mesh_Reader::invariant() const {
  bool ldum = check_dims() && (point_coords.size() == npoints) && (ipar.size() == ncells) &&
              (imat_index.size() == ncells) && (irgn_vb_index.size() == nvb_faces) &&
              (ipar_vb.size() == nvb_faces) && (ipar_rb.size() == nrb_faces);

  for (unsigned i = 0; i < ncells; ++i) {
    ldum = ldum && ipar[i].size() == nvrtx;
    for (size_t j = 0; j < nvrtx; ++j)
      ldum = ldum && ipar[i][j] < npoints;
  }
  for (size_t i = 0; i < nvb_faces; ++i) {
    ldum = ldum && ipar_vb[i].size() == nvrpf;
    for (size_t j = 0; j < nvrpf; ++j)
      ldum = ldum && ipar_vb[i][j] < npoints;
  }
  for (size_t i = 0; i < nrb_faces; ++i) {
    ldum = ldum && ipar_rb[i].size() == nvrpf;
    for (size_t j = 0; j < nvrpf; ++j)
      ldum = ldum && ipar_rb[i][j] < npoints;
  }
  return ldum;
}

bool Hex_Mesh_Reader::check_dims() const {
  bool ldum = (npoints > 0) && (ncells > 0) &&
              ((ndim == 3 && nvrtx == 8 && nvrpf == 4) || (ndim == 2 && nvrtx == 4 && nvrpf == 2) ||
               (ndim == 1 && nvrtx == 2 && nvrpf == 1)) &&
              nmat > 0;
  return ldum;
}

} // end namespace rtt_meshReaders

//------------------------------------------------------------------------------------------------//
// end of Hex_Mesh_Reader.cc
//------------------------------------------------------------------------------------------------//
