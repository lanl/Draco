//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/RTT_Mesh_Reader.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Mesh_Reader library.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_RTT_Format_Reader_RTT_Mesh_Reader_hh
#define rtt_RTT_Format_Reader_RTT_Mesh_Reader_hh

#include "RTT_Format_Reader.hh"
#include "meshReaders/Mesh_Reader.hh"

namespace rtt_RTT_Format_Reader {

//================================================================================================//
/*!
 * \class RTT_Mesh_Reader
 *
 * \brief An input routine to parse an RTT Format mesh file using the DRACO meshReaders standard
 *        interface.
 *
 * \sa The RTT_Mesh_Reader class is a derived-type of the Mesh_Reader abstract base class specified
 *     in the meshReaders package. Packages using the RTT_Mesh_Reader should thus include the
 *     RTT_Mesh_Reader.hh decorator file located in the meshReaders directory to resolve the
 *     namespace. RTT_Mesh_Reader contains a RTT_Format_Reader as a private data member, so none of
 *     the RTT_Format_Reader class public accessor functions are accessible.
 */
//================================================================================================//
class RTT_Mesh_Reader : public rtt_meshReaders::Mesh_Reader {

private:
  // NESTED CLASSES AND TYPEDEFS
  using string = std::string;
  using set_int = std::set<int>;
  using vector_int = std::vector<int>;
  using vector_vector_int = std::vector<std::vector<int>>;
  using vector_vector_vector_int = std::vector<std::vector<std::vector<int>>>;
  using vector_vector_dbl = std::vector<std::vector<double>>;
  using set_uint = std::set<unsigned int>;
  using vector_uint = std::vector<unsigned int>;
  using vector_vector_uint = std::vector<std::vector<unsigned int>>;
  using vector_vector_vector_uint = std::vector<std::vector<std::vector<unsigned int>>>;

  // DATA
  std::shared_ptr<RTT_Format_Reader> rttMesh;
  std::vector<std::shared_ptr<rtt_mesh_element::Element_Definition>> element_defs;
  std::vector<rtt_mesh_element::Element_Definition::Element_Type> element_types;
  std::vector<rtt_mesh_element::Element_Definition::Element_Type> unique_element_types;

public:
  /*!
   * \brief Constructs an RTT_Mesh_Reader class object.
   * \param RTT_File Mesh file name.
   */
  explicit RTT_Mesh_Reader(string const &RTT_File)
      : rttMesh(new RTT_Format_Reader(RTT_File)), element_defs(), element_types(),
        unique_element_types() {
    transform2CGNS();
  }

  //! Destroys an RTT_Mesh_Reader class object
  ~RTT_Mesh_Reader() override = default;
  RTT_Mesh_Reader(RTT_Mesh_Reader const &rhs) = delete;
  RTT_Mesh_Reader(RTT_Mesh_Reader &&rhs) noexcept = delete;
  RTT_Mesh_Reader &operator=(RTT_Mesh_Reader const &rhs) = delete;
  RTT_Mesh_Reader &operator=(RTT_Mesh_Reader &&rhs) noexcept = delete;

  // ACCESSORS

  /*!
   * \brief Returns the coordinate values for each of the nodes.
   * \return The coordinate values for the nodes.
   */
  vector_vector_dbl get_node_coords() const override { return rttMesh->get_nodes_coords(); }

  /*!
   * \brief Returns the problem coordinate units (e.g, cm).
   * \return Coordinate units.
   */
  string get_node_coord_units() const override { return rttMesh->get_dims_coor_units(); }

  /*!
   * \brief Returns the topological dimenstion (1, 2 or 3).
   * \return Topological dimension.
   */
  size_t get_dims_ndim() const override { return rttMesh->get_dims_ndim(); }

  size_t get_dims_ncells() const { return rttMesh->get_dims_ncells(); }

  size_t get_dims_nsides() const { return rttMesh->get_dims_nsides(); }

  vector_vector_uint get_element_nodes() const override;

  /*!
   * \brief Returns the element (i.e., sides and cells) types (e.g., TRI_3 and TETRA_4).
   * \return Element definitions.
   */
  std::vector<rtt_mesh_element::Element_Definition::Element_Type>
  get_element_types() const override {
    return element_types;
  }

  std::vector<std::shared_ptr<rtt_mesh_element::Element_Definition>>
  get_element_defs() const override {
    return element_defs;
  }

  /*!
   * \brief Returns the unique element types (e.g., TRI_3 and TETRA_4) that are defined in the mesh
   *        file.
   * \return Element definitions.
   */
  std::vector<rtt_mesh_element::Element_Definition::Element_Type>
  get_unique_element_types() const override {
    return unique_element_types;
  }

  std::map<string, set_uint> get_node_sets() const override;

  std::map<string, set_uint> get_element_sets() const override;

  /*!
   * \brief Returns the mesh file title.
   * \return Mesh file title.
   */
  string get_title() const override { return rttMesh->get_header_title(); }

  bool invariant() const override;

  // IMPLEMENTATION

private:
  void transform2CGNS();
};

} // end namespace rtt_RTT_Format_Reader

#endif // rtt_RTT_Format_Reader_RTT_Mesh_Reader_hh

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/RTT_Mesh_Reader.hh
//------------------------------------------------------------------------------------------------//
