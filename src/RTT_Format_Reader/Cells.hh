//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/Cells.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader/Cells class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_RTT_Format_Reader_Cells_hh
#define rtt_RTT_Format_Reader_Cells_hh

#include "CellDefs.hh"
#include "CellFlags.hh"
#include "Nodes.hh"

namespace rtt_RTT_Format_Reader {

//! Controls parsing, storing, and accessing the data specific to the cells block of the mesh file.
class Cells {
  // typedefs
  using ifstream = std::ifstream;
  using string = std::string;
  using vector_int = std::vector<int>;
  using vector_vector_int = std::vector<std::vector<int>>;
  using vector_uint = std::vector<unsigned int>;
  using vector_vector_uint = std::vector<std::vector<unsigned int>>;

  const CellFlags &cellFlags;
  const Dims &dims;
  const CellDefs &cellDefs;
  vector_int cellType;
  vector_vector_uint nodes;
  vector_vector_int flags;

public:
  Cells(const CellFlags &cellFlags_, const Dims &dims_, const CellDefs &cellDefs_)
      : cellFlags(cellFlags_), dims(dims_), cellDefs(cellDefs_), cellType(dims.get_ncells()),
        nodes(dims.get_ncells()),
        flags(dims.get_ncells(), vector_int(dims.get_ncell_flag_types())) { /* empty */
  }
  ~Cells() = default;
  Cells(Cells const &rhs) = delete;
  Cells(Cells &&rhs) noexcept = delete;
  Cells &operator=(Cells const &rhs) = delete;
  Cells &operator=(Cells &&rhs) noexcept = delete;

  void readCells(ifstream &meshfile);
  void redefineCells();

  /*!
   * \brief Returns the cell type associated with the specified cell.
   * \param cell_numb Cell number.
   * \return The cell type.
   */
  int get_type(size_t cell_numb) const { return cellType[cell_numb]; }

  /*!
   * \brief Returns all of the node numbers for each of the cells.
   * \return The node numbers for all cells.
  */
  vector_vector_uint get_nodes() const { return nodes; }

  /*!
   * \brief Returns all of the node numbers associated with the specified cell.
   * \param cell_numb Cell number.
   * \return The cell node numbers.
   */
  vector_uint get_nodes(size_t cell_numb) const { return nodes[cell_numb]; }

  /*!
   * \brief Returns the node number associated with the specified cell and cell-node index.
   * \param cell_numb Cell number.
   * \param node_numb Cell-node index number.
   * \return The cell node number.
   */
  int get_nodes(size_t cell_numb, size_t node_numb) const { return nodes[cell_numb][node_numb]; }

  /*!
   * \brief Returns the cell flag for the specified cell and flag index
   * \param cell_numb Cell number.
   * \param flag_numb Cell flag index.
   * \return The cell flag.
   */
  int get_flags(size_t cell_numb, size_t flag_numb) const { return flags[cell_numb][flag_numb]; }

private:
  void readKeyword(ifstream &meshfile);
  void readData(ifstream &meshfile);
  void readEndKeyword(ifstream &meshfile);
};

} // end namespace rtt_RTT_Format_Reader

#endif // rtt_RTT_Format_Reader_Cells_hh

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/Cells.hh
//------------------------------------------------------------------------------------------------//
