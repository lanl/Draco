//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   RTT_Format_Reader/CellDefs.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for CellDefs library.
 * \note   Copyright (C) 2010-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_RTT_Format_Reader_CellDefs_hh
#define rtt_RTT_Format_Reader_CellDefs_hh

#include "Dims.hh"
#include <memory>
#include <set>

namespace rtt_RTT_Format_Reader {

// forward declaration.
class CellDefs;

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Controls parsing, storing, and accessing the data contained in the cell definition block
 *        of the mesh file.
 */
class CellDef {
  /* TYPEDEFS */

  using ifstream = std::ifstream;
  using string = std::string;
  using set_int = std::set<int>;
  using vector_int = std::vector<int>;
  using vector_vector_int = std::vector<std::vector<int>>;
  using set_uint = std::set<unsigned int>;
  using vector_uint = std::vector<unsigned int>;
  using vector_vector_uint = std::vector<std::vector<unsigned int>>;

  /* DATA */

  const CellDefs &cellDefs;
  const string name;
  unsigned nnodes;
  unsigned nsides;
  vector_uint side_types;
  vector_vector_uint sides;
  // Add the capability to maintain the sense of the outward normals.
  vector_vector_uint ordered_sides;
  // Mapping between the old and new cell definition nodes.
  vector_uint node_map;

public:
  CellDef(const CellDefs &cellDefs_, const string name_)
      : cellDefs(cellDefs_), name(std::move(name_)), nnodes(0), nsides(0),
        side_types(vector_uint()), sides(vector_vector_uint()), ordered_sides(0),
        node_map() { /* empty */
  }

  ~CellDef() = default;

  //! Copy constructor
  CellDef(CellDef const &rhs) = default;
  CellDef(CellDef &&rhs) = default;

  // Assignment operator - "explicitly defaulted assignment op is implicitly deleted because
  // CellsDefs is of reference type.
  // CellDef & operator=(CellDef const & rhs) = default;
  // CellDef & operator=(CellDef && rhs) = default;

  void readDef(ifstream &meshfile);
  void redefineCellDef(vector_uint const &new_side_types_,
                       vector_vector_uint const &new_ordered_sides, size_t const ndim);

public:
  /*!
   * \brief  Returns the cell definition name.
   * \return The cell definition name.
   */
  string get_name() const { return name; }

  /*!
   * \brief  Returns the number of nodes associated with the cell definition.
   * \return The number of nodes comprising the cell definition.
   */
  unsigned get_nnodes() const { return nnodes; }

  /*!
   * \brief  Returns the number of sides associated with the cell definition.
   * \return The number of sides comprising the cell definition.
   */
  unsigned get_nsides() const { return nsides; }

  /*!
   * \brief Returns the side type number associated with the cell definition specified side index.
   * \param s Side index number.
   * \return  The side type number.
   */
  int get_side_types(size_t s) const { return side_types[s]; }

  vector_uint get_all_side_types() const { return side_types; }
  vector_vector_uint get_all_sides() const { return sides; }
  vector_vector_uint get_all_ordered_sides() const { return ordered_sides; }

  /*!
   * \brief Returns the side definition of the specified side index of this cell definition with the
   *        returned cell-node indexes in sorted order.
   * \param[in] s Side index number.
   * \return The side definition (i.e., the cell-node indexes that comprise the side).
   */
  vector_uint const &get_side(size_t s) const { return sides[s]; }

  /*!
   * \brief Returns the side definition of the specified side index of this cell definition with the
   *        returned cell-node indexes ordered to preserve the right hand rule for the
   *        outward-directed normal.
   * \param[in] s Side index number.
   * \return The side definition (i.e., the cell-node indexes that comprise the side).
   */
  vector_uint const &get_ordered_side(size_t s) const { return ordered_sides[s]; }

  //-------------------------------------------------------------------------//
  /*!
   * \brief Returns the new nodes map when cell redefinition has been performed.
   * \return New nodes map.
   */
  const vector_uint &get_node_map() const { return node_map; }
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Controls parsing, storing, and accessing the data contained in the cell definition block
 *        of the mesh file.
 */
class CellDefs {
  // typedefs
  using ifstream = std::ifstream;
  using string = std::string;
  using set_int = std::set<int>;
  using vector_int = std::vector<int>;
  using vector_vector_int = std::vector<std::vector<int>>;
  using vector_vector_vector_int = std::vector<std::vector<std::vector<int>>>;
  using set_uint = std::set<unsigned int>;
  using vector_uint = std::vector<unsigned int>;
  using vector_vector_uint = std::vector<std::vector<unsigned int>>;

  const Dims &dims;
  std::vector<std::shared_ptr<CellDef>> defs;
  bool redefined;

public:
  CellDefs(const Dims &dims_)
      : dims(dims_), defs(dims.get_ncell_defs()), redefined(false) { /*empty*/
  }
  ~CellDefs() = default;

  void readCellDefs(ifstream &meshfile);
  void redefineCellDefs(vector_vector_uint const &cell_side_types,
                        std::vector<vector_vector_uint> const &cell_ordered_sides);

private:
  void readKeyword(ifstream &meshfile);
  void readDefs(ifstream &meshfile);
  void readEndKeyword(ifstream &meshfile);

public:
  /*!
   * \brief Returns the name of the specified cell definition.
   * \param[in] i Cell definition index number.
   * \return The cell definition name.
   */
  string get_name(size_t i) const {
    Require(i < defs.size());
    return defs[i]->get_name();
  }

  /*!
   * \brief Returns the specified cell definition.
   * \param[in] i Cell definition index number.
   * \return The cell definition.
   */
  const CellDef &get_cell_def(size_t i) const { return *(defs[i]); }
  std::shared_ptr<CellDef> get_def(size_t i) const { return defs[i]; }

  /*!
   * \brief Returns the number of nodes associated with the specified cell definition.
   * \param[in] i Cell definition index number.
   * \return The number of nodes comprising the cell definition.
   */
  size_t get_nnodes(size_t i) const { return defs[i]->get_nnodes(); }

  /*!
   * \brief Returns the number of sides associated with the specified cell definition.
   * \param[in] i Cell definition index number.
   * \return The number of sides comprising the cell definition.
   */
  size_t get_nsides(size_t i) const { return defs[i]->get_nsides(); }

  /*!
   * \brief Returns the side type number associated with the specified side index and cell
   *        definition.
   * \param[in] i Cell definition index number.
   * \param[in] s Side index number.
   * \return The side type number.
   */
  int get_side_types(size_t i, size_t s) const { return defs[i]->get_side_types(s); }

  /*!
   * \brief Returns the side definition associated with the specified cell definition and side index
   *        with the returned cell-node indexes in sorted order.
   * \param[in] i Cell definition index number.
   * \param[in] s Side index number.
   * \return The side definition (i.e., the cell-node indexes that comprise the side).
   */
  vector_uint const &get_side(size_t i, size_t s) const { return defs[i]->get_side(s); }

  /*!
   * \brief Returns the side definition associated with the specified cell definition and side index
   *        with the returned cell-node indexes ordered to preserve the right hand rule for the
   *        outward-directed normal.
   * \param[in] i Cell definition index number.
   * \param[in] s Side index number.
   * \return The side definition (i.e., the cell-node indexes that comprise the side).
   */
  vector_uint const &get_ordered_side(size_t i, size_t s) const {
    return defs[i]->get_ordered_side(s);
  }

  /*!
   * \brief Returns the status of the flag indicating that the cell definitions have been redefined.
   * \return The status of the redefined flag.
   */
  bool get_redefined() const { return redefined; }

  /*!
   * \brief Returns the new node map for the specified cell definition when redefinition has been
   *        performed.
   * \param[in] cell_def Cell definition index.
   * \return New cell definition node map.
   */
  const vector_uint &get_node_map(int cell_def) const { return defs[cell_def]->get_node_map(); }
};

} // end namespace rtt_RTT_Format_Reader

#endif // rtt_RTT_Format_Reader_CellDefs_hh

//------------------------------------------------------------------------------------------------//
// end of RTT_Format_Reader/CellDefs.hh
//------------------------------------------------------------------------------------------------//
