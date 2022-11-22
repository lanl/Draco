//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/Dims.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader/Dims class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_RTT_Format_Reader_Dims_hh
#define rtt_RTT_Format_Reader_Dims_hh

#include "ds++/Assert.hh"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

namespace rtt_RTT_Format_Reader {

//================================================================================================//
/*!
 * \class Dims
 * \brief Controls parsing, storing, and accessing the data contained in the dimensions block of the
 *        mesh file.
 */
//================================================================================================//
class Dims {

  // typedefs
  using ifstream = std::ifstream;
  using string = std::string;
  using vector_int = std::vector<int>;

  string coor_units;
  string prob_time_units;

  size_t ncell_defs{0};
  size_t nnodes_max{0};
  size_t nsides_max{0};
  size_t nnodes_side_max{0};

  size_t ndim{0};
  size_t ndim_topo{0};

  size_t nnodes{0};
  size_t nnode_flag_types{0};
  vector_int nnode_flags;
  size_t nnode_data{0};

  size_t nsides{0};
  size_t nside_types{0};
  vector_int side_types;
  size_t nside_flag_types{0};
  vector_int nside_flags;
  size_t nside_data{0};

  size_t ncells{0};
  size_t ncell_types{0};
  vector_int cell_types;
  size_t ncell_flag_types{0};
  vector_int ncell_flags;
  size_t ncell_data{0};

public:
  Dims()
      : coor_units(std::string()), prob_time_units(std::string()), nnode_flags(std::vector<int>()),
        side_types(std::vector<int>()), nside_flags(std::vector<int>()),
        cell_types(std::vector<int>()), ncell_flags(std::vector<int>()) { /* empty */
  }
  ~Dims() = default;
  Dims(Dims const &rhs) = delete;
  Dims(Dims &&rhs) noexcept = delete;
  Dims &operator=(Dims const &rhs) = delete;
  Dims &operator=(Dims &&rhs) noexcept = delete;

  void readDims(ifstream &meshfile);

  // units and cell definition data access
  /*!
   * \brief Returns the problem coordinate units (e.g, cm).
   * \return Coordinate units.
   */
  string get_coor_units() const { return coor_units; }
  /*!
   * \brief Returns the problem time units (e.g, shakes).
   * \return Time units.
   */
  string get_prob_time_units() const { return prob_time_units; }
  /*!
   * \brief Returns the number of unique cell type definitions.
   * \return The number of cell definitions.
   */
  size_t get_ncell_defs() const { return ncell_defs; }
  /*!
   * \brief Returns the maximum number of nodes per cell type.
   * \return The maximum number of nodes per cell type.
   */
  size_t get_nnodes_max() const { return nnodes_max; }
  /*!
   * \brief Returns the maximum number of sides per cell type.
   * \return The maximum number of sides per cell type.
   */
  size_t get_nsides_max() const { return nsides_max; }
  /*!
   * \brief Returns the maximum number of nodes per cell side.
   * \return The maximum number of nodes per cell side.
   */
  size_t get_nnodes_side_max() const { return nnodes_side_max; }

  // dimensions node data access
  /*!
   * \brief Returns the number of spatial dimensions.
   * \return The number of spatial dimensions.
   */
  size_t get_ndim() const { return ndim; }
  /*!
   * \brief Returns the number of topological dimensions.
   * \return The number of topological dimensions.
   */
  size_t get_ndim_topo() const { return ndim_topo; }
  /*!
   * \brief Returns the number of nodes.
   * \return The number of nodes.
   */
  size_t get_nnodes() const { return nnodes; }
  /*!
   * \brief Returns the number of node flag types.
   * \return The number of node flag types.
   */
  size_t get_nnode_flag_types() const { return nnode_flag_types; }
  /*!
   * \brief Returns the number of node flags for the specified node flag type.
   * \param i Node flag type number.
   * \return The number of node flags.
   */
  size_t get_nnode_flags(size_t i) const { return nnode_flags[i]; }
  /*!
   * \brief Returns the number of node data fields.
   * \return The number of node data fields.
   */
  size_t get_nnode_data() const { return nnode_data; }

  // dimensions side data access
  /*!
   * \brief Returns the number of sides read.
   * \return The number of sides.
   */
  size_t get_nsides() const { return nsides; }
  /*!
   * \brief Returns the number of side types that are present in the "sides"
   *        block.
   * \return The number of side types.
   */
  size_t get_nside_types() const { return nside_types; }
  /*!
   * \brief Returns the side type index for the specified side type.
   * \param i Side type number.
   * \return The side type index.
   */
  int get_side_types(size_t i) const { return side_types[i]; }
  /*!
   * \brief Returns the number of side flag types.
   * \return The number of side flag types.
   */
  size_t get_nside_flag_types() const { return nside_flag_types; }
  /*!
   * \brief Returns the number of side flags for the specified side flag type.
   * \param i Side flag type number.
   * \return The number of side flags.
   */
  size_t get_nside_flags(size_t i) const { return nside_flags[i]; }
  /*!
   * \brief Returns the number of side data fields.
   * \return The number of side data fields.
   */
  size_t get_nside_data() const { return nside_data; }

  // cell data access
  /*!
   * \brief Returns the number of cells.
   * \return The number of cells.
   */
  size_t get_ncells() const { return ncells; }
  /*!
   * \brief Returns the number of cell types that are present in the "cells"
   *        block.
   * \return The number of cell types.
   */
  size_t get_ncell_types() const { return ncell_types; }
  /*!
   * \brief Returns the cell type index for the specified cell type.
   * \param i Cell type number.
   * \return The cell type index.
   */
  int get_cell_types(size_t i) const { return cell_types[i]; }
  /*!
   * \brief Returns the number of cell flag types.
   * \return The number of cell flag types.
   */
  size_t get_ncell_flag_types() const { return ncell_flag_types; }
  /*!
   * \brief Returns the number of cell flags for the specified cell flag type.
   * \param i Cell flag type number.
   * \return The number of cell flags.
   */
  size_t get_ncell_flags(size_t i) const { return ncell_flags[i]; }
  /*!
   * \brief Returns the number of cell data fields.
   * \return The number of cell data fields.
   */
  size_t get_ncell_data() const { return ncell_data; }

  /*!
   * \brief Validates the specified side type.
   * \param sidetype Side type number.
   * \return The existence of the side type.
   */
  bool allowed_side_type(size_t sidetype) const {
    return side_types.end() !=
           std::find(side_types.begin(), side_types.end(), static_cast<int>(sidetype));
  }
  /*!
   * \brief Validates the specified cell type.
   * \param celltype Cell type number.
   * \return The existence of the cell type.
   */
  bool allowed_cell_type(size_t celltype) const {
    return cell_types.end() !=
           std::find(cell_types.begin(), cell_types.end(), static_cast<int>(celltype));
  }

private:
  void readKeyword(ifstream &meshfile);
  void readUnits(ifstream &meshfile);
  void readCellDefs(ifstream &meshfile);
  void readDimensions(ifstream &meshfile);
  void readNodes(ifstream &meshfile);
  void readSides(ifstream &meshfile);
  void readCells(ifstream &meshfile);
  void readEndKeyword(ifstream &meshfile);
};

} // end namespace rtt_RTT_Format_Reader

#endif // rtt_RTT_Format_Reader_Dims_hh

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/Dims.hh
//------------------------------------------------------------------------------------------------//
