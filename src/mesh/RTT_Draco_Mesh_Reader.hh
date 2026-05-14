//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/RTT_Draco_Mesh_Reader.hh
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \brief  RTT_Draco_Mesh_Reader header file.
 * \note   Copyright (C) 2018-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_mesh_RTT_Draco_Mesh_Reader_hh
#define rtt_mesh_RTT_Draco_Mesh_Reader_hh

#include "Draco_Mesh_Reader.hh"
#include "RTT_Format_Reader/RTT_Format_Reader.hh"

namespace rtt_mesh {

//================================================================================================//
/*!
 * \class RTT_Draco_Mesh_Reader
 *
 * \brief Wrap RTT file format reader and provide data to Draco_Mesh_Builder.
 */
//================================================================================================//

class RTT_Draco_Mesh_Reader : public Draco_Mesh_Reader {
private:
  // >>> DATA

  const std::string filename;

  std::shared_ptr<rtt_RTT_Format_Reader::RTT_Format_Reader> rtt_reader;

public:
  //! Constructor
  explicit RTT_Draco_Mesh_Reader(const std::string &filename_);

  // >>> SERVICES

  void read_mesh() override;

  // >>> ACCESSORS

  [[nodiscard]] unsigned get_numdim() const override { return rtt_reader->get_dims_ndim(); }
  [[nodiscard]] size_t get_numcells() const override { return rtt_reader->get_dims_ncells(); }
  [[nodiscard]] size_t get_numnodes() const override { return rtt_reader->get_dims_nnodes(); }
  [[nodiscard]] std::vector<double> get_nodecoord(size_t node) const override {
    return rtt_reader->get_nodes_coords(node);
  }
  [[nodiscard]] size_t get_numsides() const override { return rtt_reader->get_dims_nsides(); }
  [[nodiscard]] unsigned get_sideflag(size_t side) const override {
    return rtt_reader->get_sides_flags(side, 0);
  }
  [[nodiscard]] std::vector<unsigned> get_sidenodes(size_t side) const override {
    return rtt_reader->get_sides_nodes(side);
  }

  // accessors with deferred implementation

  /*!
  * \brief Get the number of nodes for a cell.
  * \param[in] cell index of cell (0-based)
  * \return number of nodes for cell
  */
  [[nodiscard]] unsigned get_celltype(size_t cell) const override;

  /*!
   * \brief Get the vector of node indices for a side.
   * \param[in] side index of side (0-based)
   * \return number of nodes for side
   */
  [[nodiscard]] size_t get_sidetype(size_t side) const override;
  [[nodiscard]] std::vector<unsigned> get_cellnodes(size_t cell) const override;
  [[nodiscard]] std::vector<unsigned> get_cellfacenodes(size_t cell, size_t face) const;
};

} // end namespace rtt_mesh

#endif // rtt_mesh_RTT_Draco_Mesh_Reader_hh

//------------------------------------------------------------------------------------------------//
// end of mesh/RTT_Draco_Mesh_Reader.hh
//------------------------------------------------------------------------------------------------//
