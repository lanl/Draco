//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/RTT_Format_Reader.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader library.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_RTT_Format_Reader_RTT_Format_Reader_hh
#define rtt_RTT_Format_Reader_RTT_Format_Reader_hh

#include "CellData.hh"
#include "CellDataIDs.hh"
#include "Header.hh"
#include "NodeData.hh"
#include "NodeDataIDs.hh"
#include "SideData.hh"
#include "SideDataIDs.hh"

namespace rtt_RTT_Format_Reader {

//================================================================================================//
/*!
 * \brief A generalized input routine to parse an RTT Format mesh file.
 *
 * \sa The RTT_Format_Reader class constructor automatically instantiates and executes the readMesh
 *     member function used to parse the mesh data.  Accessor functions are provided for all of the
 *     remaining member classes to allow data retrieval. The \ref overview_rtt_format_reader page
 *     presents a summary of the capabilities provided by the class.
 */
//================================================================================================//

class RTT_Format_Reader {

private:
  // NESTED CLASSES AND TYPEDEFS
  using ifstream = std::ifstream;
  using string = std::string;
  using vector_int = std::vector<int>;
  using vector_dbl = std::vector<double>;
  using vector_vector_dbl = std::vector<std::vector<double>>;
  using vector_str = std::vector<string>;
  using set_uint = std::set<unsigned int>;
  using vector_uint = std::vector<unsigned int>;
  using vector_vector_uint = std::vector<std::vector<unsigned int>>;

  // DATA
  Header header;
  Dims dims;
  std::shared_ptr<NodeFlags> spNodeFlags;
  std::shared_ptr<SideFlags> spSideFlags;
  std::shared_ptr<CellFlags> spCellFlags;
  std::shared_ptr<NodeDataIDs> spNodeDataIds;
  std::shared_ptr<SideDataIDs> spSideDataIds;
  std::shared_ptr<CellDataIDs> spCellDataIds;
  std::shared_ptr<CellDefs> spCellDefs;
  std::shared_ptr<Nodes> spNodes;
  std::shared_ptr<Sides> spSides;
  std::shared_ptr<Cells> spCells;
  std::shared_ptr<NodeData> spNodeData;
  std::shared_ptr<SideData> spSideData;
  std::shared_ptr<CellData> spCellData;

public:
  //! Constructor
  explicit RTT_Format_Reader(string const &RTT_File);

  //! Destructor
  ~RTT_Format_Reader() = default;
  RTT_Format_Reader(RTT_Format_Reader const &rhs) = delete;
  RTT_Format_Reader(RTT_Format_Reader &&rhs) noexcept = delete;
  RTT_Format_Reader &operator=(RTT_Format_Reader const &rhs) = delete;
  RTT_Format_Reader &operator=(RTT_Format_Reader &&rhs) noexcept = delete;

  // ACCESSORS

  /*!
   * \brief Returns the mesh file version number.
   * \return Version number.
   */
  string get_header_version() const { return header.get_version(); }

  /*!
   * \brief Returns the mesh file title.
   * \return Title.
   */
  string get_header_title() const { return header.get_title(); }

  /*!
   * \brief Returns the mesh file date.
   * \return Date the mesh file was generated.
   */
  string get_header_date() const { return header.get_date(); }

  /*!
   * \brief Returns the mesh file cycle number.
   * \return Cycle number.
   */
  size_t get_header_cycle() const { return header.get_cycle(); }

  /*!
   * \brief Returns the mesh file problem time.
   * \return Problem time.
   */
  double get_header_time() const { return header.get_time(); }

  /*!
   * \brief Returns the number of comment lines in the mesh file.
   * \return The number of comment lines.
   */
  size_t get_header_ncomments() const { return header.get_ncomments(); }

  /*!
   * \brief Returns the specified comment line from the mesh file.
   * \param i Line number of the comment to be returned.
   * \return The comment line.
   */
  string get_header_comments(size_t i) const { return header.get_comments(i); }

  // dimensions units and cell definition data access
  /*!
   * \brief Returns the problem coordinate units (e.g, cm).
   * \return  Coordinate units.
   */
  string get_dims_coor_units() const { return dims.get_coor_units(); }

  /*!
   * \brief Returns the problem time units (e.g, shakes).
   * \return Time units.
   */
  string get_dims_prob_time_units() const { return dims.get_prob_time_units(); }

  /*!
   * \brief Returns the number of unique cell type definitions.
   * \return The number of cell definitions.
   */
  size_t get_dims_ncell_defs() const { return dims.get_ncell_defs(); }

  /*!
   * \brief Returns the maximum number of nodes per cell type.
   * \return The maximum number of nodes per cell type.
   */
  size_t get_dims_nnodes_max() const { return dims.get_nnodes_max(); }

  /*!
   * \brief Returns the maximum number of sides per cell type.
   * \return The maximum number of sides per cell type.
   */
  size_t get_dims_nsides_max() const { return dims.get_nsides_max(); }

  /*!
   * \brief Returns the maximum number of nodes per cell side.
   * \return The maximum number of nodes per cell side.
   */
  size_t get_dims_nnodes_side_max() const { return dims.get_nnodes_side_max(); }

  // dimensions node data access
  /*!
   * \brief Returns the number of spatial dimensions.
   * \return The number of spatial dimensions.
   */
  unsigned get_dims_ndim() const {
    Check(dims.get_ndim() < UINT_MAX);
    return static_cast<unsigned>(dims.get_ndim());
  }

  /*!
   * \brief Returns the number of topological dimensions.
   * \return The number of topological dimensions.
   */
  size_t get_dims_ndim_topo() const { return dims.get_ndim_topo(); }

  /*!
   * \brief Returns the number of nodes.
   * \return The number of nodes.
   */
  size_t get_dims_nnodes() const { return dims.get_nnodes(); }

  /*!
   * \brief Returns the number of node flag types.
   * \return The number of node flag types.
   */
  size_t get_dims_nnode_flag_types() const { return dims.get_nnode_flag_types(); }

  /*!
   * \brief Returns the number of node flags for the specified node flag type.
   * \param i Node flag type number.
   * \return The number of node flags.
   */
  size_t get_dims_nnode_flags(size_t i) const { return dims.get_nnode_flags(i); }

  /*!
   * \brief Returns the number of node data fields.
   * \return The number of node data fields.
   */
  size_t get_dims_nnode_data() const { return dims.get_nnode_data(); }

  // dimensions side data access
  /*!
   * \brief Returns the number of sides.
   * \return The number of sides.
   */
  size_t get_dims_nsides() const { return dims.get_nsides(); }

  /*!
   * \brief Returns the number of side types that are present in the "sides" block.
   * \return The number of side types.
   */
  size_t get_dims_nside_types() const { return dims.get_nside_types(); }
  /*!
   * \brief Returns the side type index for the specified side type.
   * \param i Side type number.
   * \return The side type index.
   */
  int get_dims_side_types(size_t i) const { return dims.get_side_types(i); }
  /*!
   * \brief Returns the number of side flag types.
   * \return The number of side flag types.
   */
  size_t get_dims_nside_flag_types() const { return dims.get_nside_flag_types(); }
  /*!
   * \brief Returns the number of side flags for the specified side flag type.
   * \param i Side flag type number.
   * \return The number of side flags.
   */
  size_t get_dims_nside_flags(size_t i) const { return dims.get_nside_flags(i); }
  /*!
   * \brief Returns the number of side data fields.
   * \return The number of side data fields.
   */
  size_t get_dims_nside_data() const { return dims.get_nside_data(); }

  // dimensions cell data access
  /*!
   * \brief Returns the number of cells.
   */
  size_t get_dims_ncells() const { return dims.get_ncells(); }
  /*!
   * \brief Returns the number of cell types that are present in the "cells" block.
   * \return The number of cell types.
   */
  size_t get_dims_ncell_types() const { return dims.get_ncell_types(); }
  /*!
   * \brief Returns the cell type index for the specified cell type.
   * \param i Cell type number.
   * \return The cell type index.
   */
  int get_dims_cell_types(size_t i) const { return dims.get_cell_types(i); }
  /*!
   * \brief Returns the number of cell flag types.
   * \return The number of cell flag types.
   */
  size_t get_dims_ncell_flag_types() const { return dims.get_ncell_flag_types(); }

  /*!
   * \brief Returns the number of cell flags for the specified cell flag type.
   * \param i Cell flag type number.
   * \return The number of cell flags.
   */
  size_t get_dims_ncell_flags(size_t i) const { return dims.get_ncell_flags(i); }

  /*!
   * \brief Returns the number of cell data fields.
   * \return The number of cell data fields.
   */
  size_t get_dims_ncell_data() const { return dims.get_ncell_data(); }

  // node flags access
  /*!
   * \brief Returns the name of specified node flag type.
   * \param flagtype Node flag type number.
   * \return The node flag type name.
   */
  string get_node_flags_flag_type(size_t flagtype) const {
    return spNodeFlags->get_flag_type(flagtype);
  }

  /*!
   * \brief Returns the index to the node flag type that contains the specified string.
   * \param desired_flag_type Flag type.
   * \return The node flag type index.
   */
  int get_node_flags_flag_type_index(string &desired_flag_type) const {
    return spNodeFlags->get_flag_type_index(desired_flag_type);
  }

  /*!
   * \brief Returns the node flag number associated with the specified node flag type and node flag
   *        index.
   * \param flagtype Node flag type number.
   * \param flag_index Node flag index.
   * \return The node flag number.
   */
  int get_node_flags_flag_number(size_t flagtype, size_t flag_index) const {
    return spNodeFlags->get_flag_number(flagtype, flag_index);
  }

  /*!
   * \brief Returns the number of node flags for the specified node flag type.
   * \param flagtype Node flag type number.
   * \return The number of node flags.
   */
  size_t get_node_flags_flag_size(size_t flagtype) const {
    return spNodeFlags->get_flag_size(flagtype);
  }

  /*!
   * \brief Returns the node flag name associated with the specified node flag type and node flag
   *        type index.
   * \param flagtype Node flag type number.
   * \param flag_index Node flag index.
   * \return The node flag name.
   */
  string get_node_flags_flag_name(size_t flagtype, size_t flag_index) const {
    return spNodeFlags->get_flag_name(flagtype, flag_index);
  }

  // side flags access
  /*!
   * \brief Returns the name of specified side flag type
   * \param flagtype Side flag type number.
   * \return The side flag type name.
   */
  string get_side_flags_flag_type(size_t flagtype) const {
    return spSideFlags->get_flag_type(flagtype);
  }

  /*!
   * \brief Returns the index to the side flag type that contains the specified string.
   * \param desired_flag_type Flag type.
   * \return The side flag type index.
   */
  int get_side_flags_flag_type_index(string &desired_flag_type) const {
    return spSideFlags->get_flag_type_index(desired_flag_type);
  }

  /*!
   * \brief Returns the side flag number associated with the specified side flag type and side flag
   *        index.
   * \param flagtype Side flag index.
   * \param flag_index Side flag index.
   * \return The side flag number.
   */
  int get_side_flags_flag_number(size_t flagtype, size_t flag_index) const {
    return spSideFlags->get_flag_number(flagtype, flag_index);
  }

  /*!
   * \brief Returns the number of side flags for the specified side flag type.
   * \param flagtype Side flag type number.
   * \return The number of side flags.
   */
  size_t get_side_flags_flag_size(size_t flagtype) const {
    return spSideFlags->get_flag_size(flagtype);
  }

  /*!
   * \brief Returns the side flag name associated with the specified side flag index and side flag
   *        type.
   * \param flagtype Side flag index.
   * \param flag_index Side flag index.
   * \return The side flag name.
   */
  string get_side_flags_flag_name(size_t flagtype, size_t flag_index) const {
    return spSideFlags->get_flag_name(flagtype, flag_index);
  }

  // cell flags access
  /*!
   * \brief Returns the name of specified cell flag type
   * \param flagtype Cell flag type number.
   * \return The cell flag type name.
   */
  string get_cell_flags_flag_type(size_t flagtype) const {
    return spCellFlags->get_flag_type(flagtype);
  }

  /*!
   * \brief Returns the index to the cell flag type that contains the specified string.
   * \param desired_flag_type Flag type.
   * \return The cell flag type index.
   */
  int get_cell_flags_flag_type_index(string &desired_flag_type) const {
    return spCellFlags->get_flag_type_index(desired_flag_type);
  }

  /*!
   * \brief Returns the cell flag number associated with the specified cell flag type and cell flag
   *        index.
   * \param flagtype Cell flag type number.
   * \param flag_index Cell flag index.
   * \return The cell flag number.
   */
  int get_cell_flags_flag_number(size_t flagtype, size_t flag_index) const {
    return spCellFlags->get_flag_number(flagtype, flag_index);
  }

  /*!
   * \brief Returns the number of cell flags for the specified cell flag type.
   * \param flagtype Cell flag type number.
   * \return The number of cell flags.
   */
  size_t get_cell_flags_flag_size(size_t flagtype) const {
    return spCellFlags->get_flag_size(flagtype);
  }

  /*!
   * \brief Returns the cell flag name associated with the specified cell flag type and cell flag
   *        index.
   * \param flagtype Cell flag type number.
   * \param flag_index Cell flag index.
   * \return The cell flag name.
   */
  string get_cell_flags_flag_name(size_t flagtype, size_t flag_index) const {
    return spCellFlags->get_flag_name(flagtype, flag_index);
  }

  // node data ids access
  /*!
   * \brief Returns the specified node_data_id name.
   * \param id_numb node_data_id index number.
   * \return The node_data_id name.
   */
  string get_node_data_id_name(size_t id_numb) const {
    return spNodeDataIds->get_data_id_name(id_numb);
  }

  /*!
   * \brief Returns the units associated with the specified node_data_id.
   * \param id_numb node_data_id index number.
   * \return The node_data_id units.
   */
  string get_node_data_id_units(size_t id_numb) const {
    return spNodeDataIds->get_data_id_units(id_numb);
  }

  // side data ids access
  /*!
   * \brief Returns the specified side_data_id name.
   * \param id_numb side_data_id index number.
   * \return The side_data_id name.
   */
  string get_side_data_id_name(size_t id_numb) const {
    return spSideDataIds->get_data_id_name(id_numb);
  }

  /*!
   * \brief Returns the units associated with the specified side_data_id.
   * \param id_numb side_data_id index number.
   * \return The side_data_id units.
   */
  string get_side_data_id_units(size_t id_numb) const {
    return spSideDataIds->get_data_id_units(id_numb);
  }

  // cell data ids access
  /*!
   * \brief Returns the specified cell_data_id name.
   * \param id_numb cell_data_id index number.
   * \return The cell_data_id name.
   */
  string get_cell_data_id_name(size_t id_numb) const {
    return spCellDataIds->get_data_id_name(id_numb);
  }

  /*!
   * \brief Returns the units associated with the specified cell_data_id.
   * \param id_numb cell_data_id index number.
   * \return The cell_data_id units.
   */
  string get_cell_data_id_units(size_t id_numb) const {
    return spCellDataIds->get_data_id_units(id_numb);
  }

  // cell definitions access
  /*!
   * \brief Returns the name of the specified cell definition.
   * \param i Cell definition index number.
   * \return The cell definition name.
   */
  string get_cell_defs_name(size_t i) const { return spCellDefs->get_name(i); }

  /*!
   * \brief Returns the specified cell definition.
   * \param i Cell definition index number.
   * \return The cell definition.
   */
  const CellDef &get_cell_defs_cell_def(int i) const { return spCellDefs->get_cell_def(i); }

  std::shared_ptr<CellDef> get_cell_defs_def(int i) const { return spCellDefs->get_def(i); }
  /*!
   * \brief Returns the number of nodes associated with the specified cell definition.
   * \param i Cell definition index number.
   * \return The number of nodes comprising the cell definition.
   */
  size_t get_cell_defs_nnodes(size_t i) const { return spCellDefs->get_nnodes(i); }

  /*!
   * \brief Returns the number of sides associated with the specified cell definition.
   * \param i Cell definition index number.
   * \return The number of sides comprising the cell definition.
   */
  size_t get_cell_defs_nsides(size_t i) const { return spCellDefs->get_nsides(i); }

  /*!
   * \brief Returns the side type number associated with the specified side index and cell
   *        definition.
   * \param i Cell definition index number.
   * \param s Side index number.
   * \return The side type number.
   */
  int get_cell_defs_side_types(size_t i, size_t s) const {
    return spCellDefs->get_side_types(i, s);
  }

  /*!
   * \brief Returns the side definition associated with the specified cell definition and side index
   *        with the returned cell-node indexes in sorted order.
   * \param i Cell definition index number.
   * \param s Side index number.
   * \return The side definition (i.e., the cell-node indexes that comprise the side).
   */
  vector_uint const &get_cell_defs_side(size_t i, size_t s) const {
    return spCellDefs->get_side(i, s);
  }

  /*!
   * \brief Returns the side definition associated with the specified cell definition and side index
   *        with the returned cell-node indexes ordered to preserve the right hand rule for the
   *        outward-directed normal.
   * \param i Cell definition index number.
   * \param s Side index number.
   * \return The side definition (i.e., the cell-node indexes that comprise the side).
   */
  vector_uint const &get_cell_defs_ordered_side(size_t i, size_t s) const {
    return spCellDefs->get_ordered_side(i, s);
  }

  /*!
   * \brief Returns the status of the flag indicating that the cell definitions have been redefined.
   * \return The status of the redefined flag.
   */
  bool get_cell_defs_redefined() const { return spCellDefs->get_redefined(); }

  /*!
   * \brief Returns the new node map for the specified cell definition when redefinition has been
   *        performed.
   * \param cell_def Cell definition index.
   * \return New cell definition node map.
   */
  const vector_uint &get_cell_defs_node_map(int cell_def) const {
    return spCellDefs->get_node_map(cell_def);
  }

  /*!
   * \brief Returns the coordinate values for each of the nodes.
   * \return The coordinate values for the nodes.
   */
  vector_vector_dbl get_nodes_coords() const { return spNodes->get_coords(); }

  /*!
   * \brief Returns all of the coordinate values for the specified node.
   * \param node_numb Node number.
   * \return The node coordinate values.
   */
  vector_dbl get_nodes_coords(size_t node_numb) const { return spNodes->get_coords(node_numb); }

  /*!
   * \brief Returns the coordinate value for the specified node and direction (i.e., x, y, and z).
   * \param node_numb Node number.
   * \param coord_index Coordinate index number (x = 0, y = 1, z = 2).
   * \return The node coordinate value.
   */
  double get_nodes_coords(size_t node_numb, size_t coord_index) const {
    return spNodes->get_coords(node_numb, coord_index);
  }

  /*!
   * \brief Returns the node parent for the specified node.
   * \param node_numb Node number.
   * \return The node parent.
   */
  int get_nodes_parents(size_t node_numb) const { return spNodes->get_parents(node_numb); }

  /*!
   * \brief Returns the node flag for the specified node and flag.
   * \param node_numb Node number.
   * \param flag_numb Node flag number.
   * \return The node flag.
   */
  int get_nodes_flags(size_t node_numb, size_t flag_numb) const {
    return spNodes->get_flags(node_numb, flag_numb);
  }

  // sides access
  /*!
   * \brief Returns the side type associated with the specified side.
   * \param side_numb Side number.
   * \return The side type.
   */
  int get_sides_type(size_t side_numb) const { return spSides->get_type(side_numb); }

  /*!
   * \brief Returns the node numbers associated with each side.
   * \return The node numbers for all of the sides.
   */
  vector_vector_uint get_sides_nodes() const { return spSides->get_nodes(); }

  /*!
   * \brief Returns the node numbers associated with the specified side.
   * \param side_numb Side number.
   * \return The side node numbers.
   */
  vector_uint get_sides_nodes(size_t side_numb) const { return spSides->get_nodes(side_numb); }

  /*!
   * \brief Returns the node number associated with the specified side and side-node index.
   * \param side_numb Side number.
   * \param node_numb Side-node index number.
   * \return The side node number.
   */
  unsigned get_sides_nodes(size_t side_numb, size_t node_numb) const {
    return spSides->get_nodes(side_numb, node_numb);
  }

  /*!
   * \brief Returns the side flag for the specified side and flag.
   * \param side_numb Side number.
   * \param flag_numb Side flag number.
   * \return The side flag.
   */
  unsigned get_sides_flags(size_t side_numb, size_t flag_numb) const {
    return spSides->get_flags(side_numb, flag_numb);
  }

  // cells access
  /*!
   * \brief Returns the cell type associated with the specified cell.
   * \param cell_numb Cell number.
   * \return The cell type.
   */
  int get_cells_type(size_t cell_numb) const { return spCells->get_type(cell_numb); }

  /*!
   * \brief Returns all of the node numbers for each of the cells.
   * \return The node numbers for all cells.
   */
  vector_vector_uint get_cells_nodes() const { return spCells->get_nodes(); }

  /*!
   * \brief Returns all of the node numbers associated with the specified cell.
   * \param cell_numb Cell number.
   * \return The cell node numbers.
   */
  vector_uint get_cells_nodes(size_t cell_numb) const { return spCells->get_nodes(cell_numb); }

  /*!
   * \brief Returns the node number associated with the specified cell and cell-node index.
   * \param cell_numb Cell number.
   * \param node_numb Cell-node index number.
   * \return The node number.
   */
  unsigned get_cells_nodes(size_t cell_numb, size_t node_numb) const {
    return spCells->get_nodes(cell_numb, node_numb);
  }

  /*!
   * \brief Returns the cell flag for the specified cell and flag.
   * \param cell_numb Cell number.
   * \param flag_numb Cell flag number.
   * \return The cell flag.
   */
  int get_cells_flags(size_t cell_numb, size_t flag_numb) const {
    return spCells->get_flags(cell_numb, flag_numb);
  }

  // node_data access
  /*!
   * \brief Returns all of the data field values for each of the nodes.
   * \return The data field values for each of the nodes.
   */
  vector_vector_dbl get_node_data() const { return spNodeData->get_data(); }

  /*!
   * \brief Returns all of the data field values for the specified node.
   * \param node_numb Node number.
   * \return The node data field values.
   */
  vector_dbl get_node_data(size_t node_numb) const { return spNodeData->get_data(node_numb); }

  /*!
   * \brief Returns the specified data field value for the specified node.
   * \param node_numb Node number.
   * \param data_index Data field.
   * \return The node data field value.
   */
  double get_node_data(size_t node_numb, size_t data_index) const {
    return spNodeData->get_data(node_numb, data_index);
  }

  // side_data access
  /*!
   * \brief Returns all of the data field values for each of the sides.
   * \return The data field values for each of the sides.
   */
  vector_vector_dbl get_side_data() const { return spSideData->get_data(); }

  /*!
   * \brief Returns all of the data field values for the specified side.
   * \param side_numb Side number.
   * \return The side data field values.
   */
  vector_dbl get_side_data(size_t side_numb) const { return spSideData->get_data(side_numb); }

  /*!
   * \brief Returns the specified data field value for the specified side.
   * \param side_numb Side number.
   * \param data_index Data field.
   * \return The side data field value.
   */
  double get_side_data(size_t side_numb, size_t data_index) const {
    return spSideData->get_data(side_numb, data_index);
  }

  // cell_data access
  /*!
   * \brief Returns all of the data field values for each of the cells.
   * \return The data field values for each of the cells.
   */
  vector_vector_dbl get_cell_data() const { return spCellData->get_data(); }

  /*!
   * \brief Returns all of the data field values for the specified cell.
   * \param cell_numb Cell number.
   * \return The cell data field values.
   */
  vector_dbl get_cell_data(size_t cell_numb) const { return spCellData->get_data(cell_numb); }

  /*!
   * \brief Returns the specified data field value for the specified cell.
   * \param cell_numb Cell number.
   * \param data_index Data field.
   * \return The cell data field value.
   */
  double get_cell_data(size_t cell_numb, size_t data_index) const {
    return spCellData->get_data(cell_numb, data_index);
  }

  // IMPLEMENTATION

  void reformatData(vector_vector_uint const &cell_side_types_,
                    std::vector<vector_vector_uint> const &cell_ordered_sides_);

private:
  void readMesh(const string &RTT_file);
  void readKeyword(ifstream &meshfile);
  void createMembers();
  void readFlagBlocks(ifstream &meshfile);
  void readDataIDs(ifstream &meshfile);
  void readEndKeyword(ifstream &meshfile);
};

} // end namespace rtt_RTT_Format_Reader

#endif // rtt_RTT_Format_Reader_RTT_Format_Reader_hh

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/RTT_Format_Reader.hh
//------------------------------------------------------------------------------------------------//
