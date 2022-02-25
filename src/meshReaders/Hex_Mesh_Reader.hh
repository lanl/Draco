//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   meshReaders/Hex_Mesh_Reader.hh
 * \author John McGhee
 * \date   Tue Mar  7 08:38:04 2000
 * \brief  Header file for CIC-19 Hex format mesh reader.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_meshReaders_Hex_Mesh_Reader_hh
#define rtt_meshReaders_Hex_Mesh_Reader_hh

#include "Mesh_Reader.hh"

namespace rtt_meshReaders {

//================================================================================================//
/*!
 * \class Hex_Mesh_Reader
 *
 * \brief Reads a CIC-19 Hex Mesh format mesh data file.
 *
 * \sa The rtt_mesh_element::Element_Definition class provides information on the hex, quad, and
 *     line elements used in this class. The \ref overview_rtt_format_reader page provides an
 *     overview of the other utilities in the rtt_meshReaders namespace. the \ref
 *     rtt_meshreaders_hexformat page provides a description of the Hex file format.
 */
//================================================================================================//
class Hex_Mesh_Reader : public rtt_meshReaders::Mesh_Reader {

  // NESTED CLASSES AND TYPEDEFS

  // DATA

  std::string meshfile_name;
  static std::string keyword() { return "cic19_hex_mesh"; }
  std::string version;
  unsigned npoints;
  unsigned ncells;
  unsigned nvrtx;
  unsigned nvrpf;
  unsigned ndim;
  unsigned nvb_faces;
  unsigned nrb_faces;
  unsigned nmat;
  std::vector<std::vector<double>> point_coords;
  std::vector<std::vector<unsigned>> ipar;
  std::vector<int> imat_index;
  std::vector<int> irgn_vb_index;
  std::vector<std::vector<unsigned>> ipar_vb;
  std::vector<std::vector<unsigned>> ipar_rb;

  std::map<std::string, std::set<unsigned>> node_sets;

public:
  // CREATORS

  explicit Hex_Mesh_Reader(std::string const &filename);

  // MANIPULATORS

  // ACCESSORS

  //! Returns the point coordinates.
  std::vector<std::vector<double>> get_node_coords() const override { return point_coords; }

  /*!
   * The Hex mesh format has no provision for labeling coordinate units Consequently, this method
   * always returns the default string: "unknown".
   */
  std::string get_node_coord_units() const override { return "unknown"; }

  /*!
   * The Hex mesh format has no provision for flagging nodes.  This method therefore always returns
   * a map with one entry which contains all the nodes.
   */
  std::map<std::string, std::set<unsigned>> get_node_sets() const override { return node_sets; }

  /*!
   * There is no provision in the Hex format for naming a mesh.  This function always returns the
   * default string: "Untitled -- CIC-19 Hex Mesh"
   */
  std::string get_title() const override { return "Untitled -- CIC-19 Hex Mesh"; }

  std::vector<std::vector<unsigned>> get_element_nodes() const override;

  bool invariant() const override;

  std::map<std::string, std::set<unsigned>> get_element_sets() const override;

  std::vector<rtt_mesh_element::Element_Definition::Element_Type>
  get_element_types() const override;

  std::vector<rtt_mesh_element::Element_Definition::Element_Type>
  get_unique_element_types() const override;

  size_t get_dims_ndim() const override { return ndim; }

private:
  bool check_dims() const;

  // IMPLEMENTATION
};

} // end namespace rtt_meshReaders

#endif // rtt_meshReaders_Hex_Mesh_Reader_hh

/*!
 * \page rtt_meshreaders_hexformat The CIC-19 Hex Mesh File Format
 *
 * **Introduction**

 * The CIC-19 Hex Format was developed primarily for small-scale testing and development
 * purposes. It's chief virtue is its simplicity.  A Hex mesh format file is usually only a few
 * hundred lines long, and is in ASCII text format that can easily be directly modified by the
 * developer with any text editor. Moreover, the format does not requires a sophisticated mesh
 * generator. These characteristics make the Hex format very useful for the creation and
 * manipulation of small, simple test problems in support of an initial debugging and development
 * effort.
 *
 * The format is restricted to three element types: Line elements in 1D, quadrilateral elements in
 * 2D, and hexahedra in 3D. One flag field is provided for interior elements, and one for
 * non-reflective boundary elements. Reflective boundary elements are listed separately.
 *
 * Support for reading this file format is provided by the rtt_meshReaders::Hex_Mesh_Reader class.

 * **Format Details**
 *
 * Node and element numbering is one-based on the
 * file. The format is as follows:
 * - Line 1 -- "cic19_hex_mesh" The keyword to ID the file.
 * - Line 2 -- "npoints, ncells, nvrtx, nvrpf, ndim, nvb_faces, nrb_faces, nmat" Eight integers
 *   which describe the dimensions of the mesh.  npoints is the number of points in the mesh. ncells
 *   is the number of cells. For example a 3x3 2D quadrilateral mesh would have npoints=16 and
 *   ncells=9. nvrtx is the number of vertexes on a cell (8 for hex, 4 for quad, 2 for lines). nvrpf
 *   is the number of vertexes on cell face (4 for hex, 2 for quad, and 1 for lines).  ndim is the
 *   number of spatial dimensions(3,2 or 1). nvb_faces is the number of vacuum boundary
 *   faces. nrb_faces is the number of reflective boundary faces. nmat is the number of unique flag
 *   values assigned to interior cells.
 * - Next npoints lines -- (x,y,z)mesh point coordinates, ndim real numbers on each line.
 * - Next ncells lines -- cell vertex numbers. These numbers refer to the point coordinates line
 *   numbers in the previous section. Vertex labeling follows the conventions detailed in
 *   rtt_meshReaders::Element_Definition::BAR_2, rtt_meshReaders::Element_Definition::QUAD_4, and
 *   rtt_meshReaders::Element_Definition::HEXA_8 There are nvrtx integers on each line.
 * - Next ncells/10 (+1?) lines -- Interior cell flag data. An integer flag to be associated with
 *   each interior cell. 10 integers on each line, with less on the last line as needed.
 * - Next nvb_faces lines -- vacuum boundary face vertex numbers followed by face flag. nvrpf+1
 *   integers per line.
 * - Next nrb_faces lines -- reflective boundary face vertex numbers.  nvrpf+1 integers per line.
 *
 * The "cube.mesh.in" file found in the Examples section provides an example CIC-19 Hex format mesh
 * file.
 */

/*!
 * \example meshReaders/test/cube.mesh
 *
 * The following provides an example of a 3D, 5x5x5 hexahedra CIC-19 Hex mesh format file. The mesh
 * has 125 cells, 125 vacuum boundary faces, 25 reflective boundary faces, and four cell flag
 * values.
 */

//------------------------------------------------------------------------------------------------//
// end of meshReaders/Hex_Mesh_Reader.hh
//------------------------------------------------------------------------------------------------//
