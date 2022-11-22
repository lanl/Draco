//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   viz/Ensight_Translator.hh
 * \author Thomas M. Evans
 * \date   Fri Jan 21 16:36:10 2000
 * \brief  Ensight_Translator header file.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_viz_Ensight_Translator_hh
#define rtt_viz_Ensight_Translator_hh

#include "Ensight_Stream.hh"
#include "Viz_Traits.hh"
#include "c4/C4_Functions.hh"
#include <memory>
#include <set>

namespace rtt_viz {

// helper for SFINAE
template <bool C, typename T> using enable_if_t = typename std::enable_if<C, T>::type;

//================================================================================================//
/*!
 * \enum Ensight_Cell_Types
 *
 * \brief Supported Ensight cell types.
 *
 * Values from this enumeration are assigned (by the user) to each cell in the mesh through the
 * rtt_viz::Ensight_Translator::ensight_dump() function input field argument iel_type.  See Section
 * 2.5 (Figure 2-22) in the Ensight 6.2 user manual for further information on Ensight cell types.
 *
 * \sa Ensight example page for examples of how this enumeration is used to generate the iel_type
 *     field.
 */
//================================================================================================//

enum Ensight_Cell_Types {
  point = 0,
  two_node_bar,
  three_node_bar,
  three_node_triangle,
  six_node_triangle,
  four_node_quadrangle,
  eight_node_quadrangle,
  four_node_tetrahedron,
  ten_node_tetrahedron,
  five_node_pyramid,
  thirteen_node_pyramid,
  eight_node_hexahedron,
  twenty_node_hexahedron,
  six_node_wedge,
  fifteen_node_wedge,
  unstructured,
  unstructured_3d
};

//================================================================================================//
/*!
 * \class Ensight_Translator
 *
 * \brief A translator for dumping problem data in EnSight format.
 *
 * This class is used to create data dumps that are readable by EnSight.  Currently, only ASCII
 * format dumps are supported.  The present incantation of this class is a C++ version of the
 * ensight_dump model from Dante.
 *
 * \anchor Ensight_Translator_strings
 *
 * A note about names: all strings containing names, be they problem names or data names, should not
 * contain the following characters:
 *
 * \arg ( ) [ ] + - \@ space ! \# * ^ \$ /
 *
 * Names in the ensight translation are often used to create directories, so they must not have
 * these characters in them.  Additionally, keep all names under 19 characters.  Finally, all names
 * of data should be unique.
 *
 * \anchor Ensight_Translator_description
 *
 * Ensight_Translator is not a templated class; however, the member functions are templated on field
 * types.  These types and there restrictions are as follows:
 *
 * \arg \b SSF (string scalar field) can be any field type that supports random access iterators and
 *         operator[] indexing
 * \arg \b ISF (unsigned integer scalar field) can be any field type that supports random access
 *         iterators and operator[] indexing
 * \arg \b IVF (integer vector field) uses the rtt_viz::Viz_Traits class to convert to i,j indexing
 *         to operator(int, int).  See the rtt_viz::Viz_Traits class to determine which 2D vector
 *         containers are specialized
 * \arg \b FVF (floating point vector field) uses the rtt_viz::Viz_Traits class to convert to i,j
 *         indexing to operator(int, int).  See the rtt_viz::Viz_Traits class to determine which 2D
 *         vector containers are specialized
 *
 * The vector-field (VF) types will most commonly be std::vector<std::vector> types.  Other types
 * will require specialization in the rtt_viz::Viz_Traits class.
 *
 * Data dumped to Ensight will be stored in the directory prefix argument given to the constructor
 * (Ensight_Translator::Ensight_Translator()).  The directory created will have the name
 * "prefix"_ensight. The directory is placed in the path specified by the gd_wpath argument to the
 * constructor.  This directory must exist.
 *
 * To launch Ensight: select the "prefix".case file that resides in the top-level ensight dump
 * directory from the "file/Data (reader)" menu.  Set the data file "Format" option in Ensight to
 * "case" and hit the "(Set) Geometry" button.  From there see the Ensight manual.
 */
/*!
 * \example viz/test/tstEnsight_Translator.cc
 *
 * \anchor ensight_test
 *
 * Example usage of the Ensight_Translator class. A two part, 3-D, nine cell mesh graphics dump is
 * performed.  Two vertex data field and two cell data fields are dumped.  After running the test
 * look in testproblem_ensight to see the output.
 */
//================================================================================================//

class Ensight_Translator {
public:
  // Ensight_Translator typedefs.
  using sf_int = std::vector<int>;
  using sf_double = std::vector<double>;
  using sf_string = std::vector<std::string>;
  using std_string = std::string;
  using sf2_int = std::vector<sf_int>;
  using sf3_int = std::vector<sf2_int>;
  using set_int = std::set<int>;
  using vec_set_int = std::vector<set_int>;

  using set_const_iterator = set_int::const_iterator;
  using vec_stream = std::vector<std::shared_ptr<Ensight_Stream>>;

private:
  // >>> DATA

  //! if true, geometry is static
  bool d_static_geom;

  //! if true, output geometry and varible data files in binary format.
  bool d_binary;

  //! Name of the directory where this will write files
  std_string d_dump_dir;

  //! Number of Ensight cell types.
  unsigned d_num_cell_types;

  //! Ensight cell names.
  sf_string d_cell_names;

  //! Number of vertices for a given Ensight cell type.
  sf_int d_vrtx_cnt;

  //! Cell types.
  sf_int d_cell_type_index;

  //! Vector of dump_times.
  sf_double d_dump_times;

  //! Ensight file prefixes.
  std_string d_prefix;

  //! Names of vertex data.
  sf_string d_vdata_names;

  //! Names of cell data.
  sf_string d_cdata_names;

  //! Name of case file.
  std_string d_case_filename;

  //! Name of geometry directory.
  std_string d_geo_dir;

  //! Names of vdata directories.
  sf_string d_vdata_dirs;

  //! Names of cdata directories.
  sf_string d_cdata_dirs;

  //! Geometry file stream.
  Ensight_Stream d_geom_out;

  //! Cell-data streams.
  vec_stream d_cell_out;

  //! Vertex-data streams.
  vec_stream d_vertex_out;

  //! Domain Decomposed flag
  const bool d_decomposed;

private:
  // >>> PRIVATE IMPLEMENTATION

  // Creates some of the file prefixes and names.
  void create_filenames(const std_string &prefix);

  // Write out case file.
  void write_case();

  // Write out geometry file.
  template <typename IVF, typename FVF, typename ISF>
  void write_geom(const uint32_t part_num, const std_string &part_name,
                  const rtt_viz::Viz_Traits<IVF> &ipar, const rtt_viz::Viz_Traits<FVF> &vrtx_data,
                  const sf2_int &cells_of_type, const sf_int &vertices, const ISF &g_vrtx_indeces,
                  const ISF &g_cell_indeces);

  // Write out vertex data.
  template <typename FVF>
  void write_vrtx_data(uint32_t part_num, const rtt_viz::Viz_Traits<FVF> &vrtx_data,
                       const sf_int &vertices);

  // Write out cell data.
  template <typename FVF>
  void write_cell_data(const uint32_t part_num, const rtt_viz::Viz_Traits<FVF> &cell_data,
                       const sf2_int &cells_of_type);

  // Initializer used by constructors
  void initialize(const bool graphics_continue);

public:
  // Constructor.
  template <typename SSF>
  Ensight_Translator(const std_string &prefix, std_string gd_wpath, SSF vdata_names,
                     SSF cdata_names, const bool overwrite = false, const bool static_geom = false,
                     const bool binary = false, const bool decomposed = false,
                     const double reset_time = -1.0);

  // Do an Ensight_Dump.
  template <typename ISF, typename IVF, typename SSF, typename FVF>
  enable_if_t<std::is_integral<typename ISF::value_type>::value &&
                  std::is_unsigned<typename ISF::value_type>::value,
              void>
  ensight_dump(int icycle, double time, double dt, const IVF &ipar, const ISF &iel_type,
               const ISF &cell_rgn_index, const FVF &pt_coor, const FVF &vrtx_data,
               const FVF &cell_data, const ISF &rgn_numbers, const SSF &rgn_name);

  // Opens the geometry and variable files.
  void open(const int icycle, const double time, const double dt);

  // Closes any open file streams.
  void close();

  // Write ensight data for a single part.
  template <typename ISF, typename IVF, typename FVF>
  enable_if_t<std::is_integral<typename ISF::value_type>::value &&
                  std::is_unsigned<typename ISF::value_type>::value,
              void>
  write_part(uint32_t part_num, const std_string &part_name, const IVF &ipar_in,
             const ISF &iel_type, const FVF &pt_coor_in, const FVF &vrtx_data_in,
             const FVF &cell_data_in, const ISF &g_vrtx_indices, const ISF &g_cell_indices);

  // >>> ACCESSORS

  //! Get the list of dump times.
  const sf_double &get_dump_times() const { return d_dump_times; }
};

} // end namespace rtt_viz

//------------------------------------------------------------------------------------------------//
// include template definitions so that template functions will be automatically instantiated in
// client code

#include "Ensight_Translator.i.hh"

#endif // rtt_viz_Ensight_Translator_hh

//------------------------------------------------------------------------------------------------//
// end of viz/Ensight_Translator.hh
//------------------------------------------------------------------------------------------------//
