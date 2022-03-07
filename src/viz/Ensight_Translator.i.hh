//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   viz/Ensight_Translator.i.hh
 * \author Thomas M. Evans
 * \date   Fri Jan 21 16:36:10 2000
 * \brief  Ensight_Translator template definitions.
 * \note   Copyright (C) 2015-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include <map>
#include <sstream>
#include <type_traits>

namespace rtt_viz {

//------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for Ensight_Translator.
 *
 * This constructor builds an Ensight_Translator.  The behavior of the (existing) ensight dump files
 * is controlled by the overwrite parameter.  If this is true, any existing ensight dumps (with the
 * same problem name) will be overwritten.  If overwrite is false then the ensight dumps are
 * appended.  The ensight case files are parsed to get the dump times if overwrite is false.
 *
 * \param prefix std_string giving the name of the problem
 * \param gd_wpath directory where dumps are stored
 * \param vdata_names string field containing vertex data names
 * \param cdata_names string field containing cell data names
 * \param overwrite bool that controls whether an existing ensight directory is to be appended to or
 *           overwritten.  If true, overwrites the existing ensight directory.  If false, and the
 *           ensight directory exists, the case file is appended to.  In either case, if the ensight
 *           directory does not exist it is created.  The default for overwrite is false.
 * \param static_geom optional input that if true, geometry is assumed the same across all calls to
 *           Ensight_Translator::ensight_dump.
 * \param binary If true, geometry and variable data files are output in binary format.
 * \param decomposed If true, geometry is decomposed overall all ranks
 * \param reset_time time after which to rewrite dumps, if overwrite=false
 *
 * \note If appending data (\a overwrite is false), then \a binary must be the same value as the
 * first ensight dump.  This class does NOT check for this potential error (yes, it's possible to
 * check and is left for a future exercise).
 */
template <typename SSF>
Ensight_Translator::Ensight_Translator(const std_string &prefix, std_string gd_wpath,
                                       SSF vdata_names, SSF cdata_names, const bool overwrite,
                                       const bool static_geom, const bool binary,
                                       const bool decomposed, const double reset_time)
    : d_static_geom(static_geom), d_binary(binary), d_dump_dir(std::move(gd_wpath)),
      d_num_cell_types(0), d_cell_names(), d_vrtx_cnt(0), d_cell_type_index(), d_dump_times(),
      d_prefix(), d_vdata_names(std::move(vdata_names)), d_cdata_names(std::move(cdata_names)),
      d_case_filename(), d_geo_dir(), d_vdata_dirs(), d_cdata_dirs(), d_geom_out(), d_cell_out(),
      d_vertex_out(), d_decomposed(decomposed) {
  Require(d_dump_times.empty());
  create_filenames(prefix);

  bool graphics_continue = false; // default behavior

  if (!overwrite && rtt_c4::node() == 0) {
    // then try to parse the case file.  Case files are always ascii.

    std::ifstream casefile(d_case_filename.c_str());

    if (casefile) {
      // then case file exists, so parse the dump times
      std_string key("number of steps:");
      std_string line;
      int num_steps = 0;

      for (;;) {
        std::getline(casefile, line);
        Insist(casefile.good(), "Error getting number of steps from case file!");
        if (line.find(key) == 0) {
          std::istringstream ss(line.substr(key.size()));
          ss >> num_steps;
          //std::cout << "FOUND " << num_steps
          //      << " STEPS " << std::endl;
          break;
        }
      }

      // read next three lines and discard
      std::getline(casefile, line);
      std::getline(casefile, line);
      std::getline(casefile, line);
      Insist(casefile.good(), "Error reading case file!");

      // read the dump_times

      d_dump_times.resize(num_steps);

      for (int i = 0; i < num_steps; ++i) {
        casefile >> d_dump_times[i];
        Insist(casefile.good(), "Error reading dump_times from case file!");
      }

      casefile.close();
      graphics_continue = true;

      // check if a valid reset time was provided
      if (reset_time > 0.0) {

        // find nearest graphics dump before reset time
        int idump = num_steps;
        for (int i = num_steps - 1; i >= 0; --i) {
          if (reset_time > d_dump_times[i]) {
            idump = i + 1;
            break;
          }
        }

        // truncate times at most recent graphics dump before restart
        if (idump < num_steps)
          d_dump_times.erase(d_dump_times.begin() + idump, d_dump_times.end());
      }
    }
  }

  initialize(graphics_continue);
}

//------------------------------------------------------------------------------------------------//
// ENSIGHT DUMP PUBLIC INTERFACES
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Do an Ensight dump to disk.
 *
 * Performs an ensight dump in the directory specified in Ensight_Translator::Ensight_Translator().
 *
 * \param icycle time cycle number associated with this dump
 * \param time elapsed problem time
 * \param dt current problem timestep
 *
 * \param[in] ipar_in IVF field of pointers to vertices.  Dimensioned [0:ncells-1,
 *           0:n_local_vertices_per_cell-1], where n_local_vertices_per_cell is the number of
 *           vertices that make up the cell.  ipar(i,j) maps the jth+1 vertex number, in the ith+1
 *           cell, to Ensight's "vertex number."  The "vertex number" is in [1:nvertices], so that
 *           for example, the corresponding x-coordinate is pt_coor(ipar(i,j)-1, 0).
 *           Alternatively, for 3D unstructured meshes, ipar_in is dimensioned [0:ncells-1,
 *           0:n_local_faces_per_cell-1, 0:n_local_vertices_per_face_per_cell-1].
 * \param[in] iel_type ISF field of Ensight_Cell_Types.  Dimensioned [0:ncells-1].  Each cell in the
 *           problem must be associated with a Ensight_Cell_Types enumeration object.
 * \param[in] cell_rgn_index ISF field of region identifiers.  Dimensioned [0:ncells-1].  This
 *           matches a region index to each cell in the problem.
 * \param[in] pt_coor_in FVF field of vertex coordinates. pt_coor is dimensioned [0:nvertices-1,
 *           0:ndim-1].  For each vertex point give the value in the appropriate dimension.
 * \param[in] vrtx_data_in FVF field of vertex data.  vrtx_data is dimensioned [0:nvertices-1,
 *           0:number of vertex data fields - 1].  The ordering of the second index must match the
 *           vdata_names field input argument to Ensight_Translator::Ensight_Translator().  The
 *           ordering of the first index must match the vertex ordering from pt_coor.
 * \param[in] cell_data_in FVF field of cell data.  cell_data is dimensioned [0:ncells-1, 0:number
 *           of cell data fields - 1].  The ordering of the second index must match the cdata_names
 *           field input argument to Ensight_Translator::Ensight_Translator().  The ordering of the
 *           first index must match the cell ordering from ipar.
 * \param[in] rgn_numbers ISF field of unique region ids.  This has dimensions of the number of
 *           unique values found in the cell_rgn_index field.
 * \param[in] rgn_name SSF field of unique region names.  This has the same dimensions and ordering
 *           as rgn_numbers.  In summary, rgn_numbers gives a list of the unique region ids in the
 *           problem and rgn_name gives a list of the names associated with each region id.
 *
 * \sa \ref Ensight_Translator_strings "Ensight_Translator class" for restrictions on name strings.
 * \sa \ref Ensight_Translator_description "Ensight_Translator class" for information on templated
 *           field types.
 * \sa Examples page for more details about how to do Ensight dumps.
 */

template <typename ISF, typename IVF, typename SSF, typename FVF>
enable_if_t<std::is_integral<typename ISF::value_type>::value &&
                std::is_unsigned<typename ISF::value_type>::value,
            void>
Ensight_Translator::ensight_dump(int icycle, double time, double dt, const IVF &ipar_in,
                                 const ISF &iel_type, const ISF &cell_rgn_index,
                                 const FVF &pt_coor_in, const FVF &vrtx_data_in,
                                 const FVF &cell_data_in, const ISF &rgn_numbers,
                                 const SSF &rgn_name) {
  using rtt_viz::Viz_Traits;
  using std::find;
  using std::string;
  using std::vector;

  // >>> PREPARE DATA TO SET ENSIGHT OUTPUT

  // load traits for vector field types
  Viz_Traits<IVF> ipar(ipar_in);
  Viz_Traits<FVF> pt_coor(pt_coor_in);
  Viz_Traits<FVF> vrtx_data(vrtx_data_in);
  Viz_Traits<FVF> cell_data(cell_data_in);

  // define sizes used throughout
  size_t ncells = ipar.nrows();
  size_t nvertices = pt_coor.nrows();
  Remember(size_t nrgn = rgn_name.size(););

  // Check sizes of all data.
  Require(iel_type.size() == ncells);
  Require(cell_rgn_index.size() == ncells);
  Require(cell_data.nrows() == ncells || cell_data.nrows() == 0);
  Require(vrtx_data.nrows() == nvertices || vrtx_data.nrows() == 0);
  Require(rgn_numbers.size() == nrgn);

  // create the parts list
  ISF parts_list;

  for (size_t i = 0; i < ncells; ++i) {
    auto const find_location = find(parts_list.begin(), parts_list.end(), cell_rgn_index[i]);

    if (find_location == parts_list.end())
      parts_list.push_back(cell_rgn_index[i]);
  }

  // store the number of parts
  size_t nparts = parts_list.size();

  // create the parts names
  vector<string> part_names;

  for (size_t i = 0; i < nparts; ++i) {
    auto const find_location_c = find(rgn_numbers.begin(), rgn_numbers.end(), parts_list[i]);

    if (find_location_c != rgn_numbers.end()) {
      auto index = find_location_c - rgn_numbers.begin();
      part_names.push_back(rgn_name[index]);
    } else if (find_location_c == rgn_numbers.end()) {
      Insist(0, "Didn't supply a region name!");
    }
  }

  Insist(parts_list.size() == part_names.size(), "Mismatch on part size!");
  Insist(rgn_name.size() == parts_list.size(), "Mismatch on region size!");

  // create the cells that make up each part

  // vertices_of_part[ipart] is the set of vertex indices that make up part ipart.
  vec_set_int vertices_of_part(nparts);

  // cells_of_type[ipart][itype][i] is the cell index of the i'th cell of type itype in part ipart.
  sf3_int cells_of_type(nparts);
  for (size_t i = 0; i < nparts; ++i)
    cells_of_type[i].resize(d_num_cell_types);

  // Initialize cells_of_type and vertices_of_part.

  for (size_t i = 0; i < ncells; ++i) {
    auto const find_location = find(parts_list.begin(), parts_list.end(), cell_rgn_index[i]);

    Check(find_location != parts_list.end());
    //Check(iel_type[i] >= 0);
    Check(static_cast<unsigned>(iel_type[i]) < d_num_cell_types);

    auto const ipart = find_location - parts_list.begin();

    Check(i < INT_MAX);
    cells_of_type[ipart][iel_type[i]].push_back(static_cast<int>(i));

    // get number of columns (vertices or faces) for this cell
    const size_t n_local_cols = ipar.ncols(i);

    // check for 0-depth in ipar, which indicates it is 2-dimensional
    if (ipar.ndpth(0, 0) == 0) {
      for (size_t iv = 0; iv < n_local_cols; ++iv)
        vertices_of_part[ipart].insert(ipar(i, iv, 0) - 1);
    } else {
      for (size_t j = 0; j < n_local_cols; ++j) {
        const size_t n_local_dpth = ipar.ndpth(i, j);
        for (size_t k = 0; k < n_local_dpth; ++k)
          vertices_of_part[ipart].insert(ipar(i, j, k) - 1);
      }
    }
  }

  // Form global cell and vertex indices.  These are the same as their local index, in this case.
  sf_int g_cell_indices(ncells);
  sf_int g_vrtx_indices(nvertices);

  int local_cell_offset = 0;
  int local_vert_offset = 0;

  if (d_decomposed) {
    local_cell_offset = static_cast<int>(ncells);
    local_vert_offset = static_cast<int>(nvertices);
    local_cell_offset = rtt_c4::prefix_sum(local_cell_offset) - local_cell_offset;
    local_vert_offset = rtt_c4::prefix_sum(local_vert_offset) - local_vert_offset;
  }

  for (size_t i = 0; i < ncells; ++i) {
    Check(i + local_cell_offset < INT_MAX);
    g_cell_indices[i] = static_cast<int>(i + local_cell_offset);
  }

  for (size_t i = 0; i < nvertices; ++i) {
    Check(i + local_vert_offset < INT_MAX);
    g_vrtx_indices[i] = static_cast<int>(i + local_vert_offset);
  }

  // >>> WRITE OUT DATA TO DIRECTORIES

  open(icycle, time, dt);

  for (uint32_t ipart = 0; ipart < part_names.size(); ipart++) {
    // Load vertices_of_part into a vector.
    set_int &v = vertices_of_part[ipart];
    sf_int vertices;

    for (auto &iv : v)
      vertices.push_back(iv);

    // write the geometry data
    Check(ipart + 1 < INT_MAX);
    write_geom(static_cast<int>(ipart + 1), part_names[ipart], ipar, pt_coor, cells_of_type[ipart],
               vertices, g_vrtx_indices, g_cell_indices);

    // write the vertex data
    write_vrtx_data(ipart + 1, vrtx_data, vertices);

    // write out the cell data
    write_cell_data(static_cast<int>(ipart + 1), cell_data, cells_of_type[ipart]);
  }

  close();
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Write ensight data for a single part.
 *
 * open() must be called before calling this function.
 *
 * Writes data for a single part (or "region," in ensight_dump parlance).
 *
 * \param part_num A part number to be used by Ensight.  Must be positive.
 * \param part_name A name for the part.
 * \param ipar_in See ensight_dump().
 * \param iel_type See ensight_dump().
 * \param pt_coor_in See ensight_dump().
 * \param vrtx_data_in See ensight_dump().
 * \param cell_data_in See ensight_dump().
 * \param g_vrtx_indices Global vertex indices.  These are used by Ensight as integer labels for
 *           each vertex.  Specifically, let i access the i'th value in \a vrtx_data.  Then \a
 *           g_vrtx_indices[i] gives the "global index" (or label index) for i.  This is referred to
 *           as "global" because \a g_vrtx_indices can be used to map each processor's local indices
 *           to global index space.
 * \param g_cell_indices Global cell indices.  Analogous to \a g_vrtx_indices, but for cell indices.
 *
 * \sa \ref Ensight_Translator_strings "Ensight_Translator class" for restrictions on name strings.
 *
 * \sa \ref Ensight_Translator_description "Ensight_Translator class" for information on templated
 *           field types.
 *
 * \sa Examples page for more details about how to do Ensight dumps.
 */
template <typename ISF, typename IVF, typename FVF>
enable_if_t<std::is_integral<typename ISF::value_type>::value &&
                std::is_unsigned<typename ISF::value_type>::value,
            void>
Ensight_Translator::write_part(uint32_t part_num, const std_string &part_name, const IVF &ipar_in,
                               const ISF &iel_type, const FVF &pt_coor_in, const FVF &vrtx_data_in,
                               const FVF &cell_data_in, const ISF &g_vrtx_indices,
                               const ISF &g_cell_indices) {
  Require(part_num > 0);

  using rtt_viz::Viz_Traits;
  using std::find;
  using std::string;
  using std::vector;

  // load traits for vector field types
  Viz_Traits<IVF> ipar(ipar_in);
  Viz_Traits<FVF> pt_coor(pt_coor_in);
  Viz_Traits<FVF> vrtx_data(vrtx_data_in);
  Viz_Traits<FVF> cell_data(cell_data_in);

  // define sizes used throughout
  size_t ncells = ipar.nrows();
  size_t nvertices = pt_coor.nrows();

  // Check sizes of all data.
  Require(iel_type.size() == ncells);
  Require(cell_data.nrows() == ncells || cell_data.nrows() == 0);
  Require(vrtx_data.nrows() == nvertices || vrtx_data.nrows() == 0);
  Require(g_vrtx_indices.size() == nvertices);
  Require(g_cell_indices.size() == ncells);

  // cells_of_type[itype][i] is the cell index of the i'th cell of type itype.
  sf2_int cells_of_type(d_num_cell_types);

  for (size_t i = 0; i < ncells; ++i) {
    // Check(iel_type[i] >= 0);
    Check(static_cast<unsigned>(iel_type[i]) < d_num_cell_types);
    Check(i < INT_MAX);
    cells_of_type[iel_type[i]].push_back(static_cast<int>(i));
  }

  // All vertices are output in this case.
  sf_int vertices(nvertices);
  for (size_t i = 0; i < nvertices; ++i) {
    Check(i < INT_MAX);
    vertices[i] = static_cast<int>(i);
  }

  // >>> WRITE OUT DATA TO DIRECTORIES

  // write the geometry data
  write_geom(part_num, part_name, ipar, pt_coor, cells_of_type, vertices, g_vrtx_indices,
             g_cell_indices);

  // write the vertex data
  write_vrtx_data(part_num, vrtx_data, vertices);

  // write out the cell data
  write_cell_data(part_num, cell_data, cells_of_type);
}

//------------------------------------------------------------------------------------------------//
// ENSIGHT DATA OUTPUT FUNCTIONS (PRIVATE)
//------------------------------------------------------------------------------------------------//

//! Write out data to ensight geometry file.
template <typename IVF, typename FVF, typename ISF>
void Ensight_Translator::write_geom(const uint32_t part_num, const std_string &part_name,
                                    const rtt_viz::Viz_Traits<IVF> &ipar,
                                    const rtt_viz::Viz_Traits<FVF> &pt_coor,
                                    const sf2_int &cells_of_type, const sf_int &vertices,
                                    const ISF &g_vrtx_indices, const ISF &g_cell_indices) {
  // Return if the geometry is static and we've already dumped the geometry.
  if (d_static_geom && d_dump_times.size() > 1) {
    return;
  }

  Insist(d_geom_out.is_open(), "Geometry file not open.  Must call open() before write_part().");

  size_t ndim = pt_coor.ncols(0);
  size_t nvertices = vertices.size();
  size_t g_nvertices = nvertices;
  if (d_decomposed)
    rtt_c4::global_sum(g_nvertices);

  // output part number and names
  if (rtt_c4::node() == 0) {
    d_geom_out << "part" << endl;
    d_geom_out << part_num << endl;
    d_geom_out << part_name << endl;
    d_geom_out << "coordinates" << endl;
    d_geom_out << int(g_nvertices) << endl; // #vertices in this part
  }
  d_geom_out.flush();

  // output the global vertex indices and form ens_vertex.  Enight demands that vertices be numbered
  // from 1 to the number of vertices *for this part* (nvertices).  Argghhh.  ens_vertex maps our
  // local vertex index to a vertex in [1,nvertices].

  // Form global cell and vertex indices.  These are the same as their local index, in this case.

  // set parallel offset for cell and vertices numbers
  int local_vert_offset = 0;
  if (d_decomposed) {
    local_vert_offset = static_cast<int>(nvertices);
    local_vert_offset = rtt_c4::prefix_sum(local_vert_offset) - local_vert_offset;
  }

  std::map<int, int> ens_vertex;
  for (size_t i = 0; i < nvertices; ++i) {
    d_geom_out << g_vrtx_indices[vertices[i]] << endl;

    // add 1 because ipar and Ensight demand indices that start at 1.
    Check(i + 1 + local_vert_offset < INT_MAX);
    ens_vertex[vertices[i] + 1] = static_cast<int>(i + 1 + local_vert_offset);
  }
  d_geom_out.flush();

  // output the coordinates
  for (size_t idim = 0; idim < ndim; idim++) {
    for (size_t i = 0; i < nvertices; ++i)
      d_geom_out << pt_coor(vertices[i], idim) << endl;
    d_geom_out.flush();
  }

  // ensight expects coordinates for three dimensions, so fill any remaining dimensions with zeroes
  double zero = 0.0;
  for (size_t idim = ndim; idim < 3; idim++) {
    for (size_t i = 0; i < nvertices; ++i)
      d_geom_out << zero << endl;
    d_geom_out.flush();
  }

  // for each cell type, dump the local vertex indices for each cell.
  for (unsigned type = 0; type < d_num_cell_types; type++) {
    const sf_int &c = cells_of_type[type];
    const size_t num_elem = c.size();
    size_t g_num_elem = num_elem;
    if (d_decomposed)
      rtt_c4::global_sum(g_num_elem);

    if (g_num_elem > 0) {
      if (rtt_c4::node() == 0) {
        d_geom_out << d_cell_names[type] << endl;
        d_geom_out << int(g_num_elem) << endl;
      }
      d_geom_out.flush();

      for (size_t i = 0; i < num_elem; ++i)
        d_geom_out << g_cell_indices[c[i]] << endl;
      d_geom_out.flush();

      // for n-sided polygons, Ensight requires number of nodes per cell
      if (d_cell_names[type] == "nsided") {
        for (size_t i = 0; i < num_elem; ++i) {
          d_geom_out << static_cast<int>(ipar.ncols(c[i])) << endl;
        }
        d_geom_out.flush();
      }

      // for n-faced polyhedra, Ensight requires:
      // 1.) number of faces per element
      // 2.) number of nodes per face per element
      // 3.) local node indices per face per element
      if (d_cell_names[type] == "nfaced") {

        // 3D unstructured cells must have non-zero depth
        Check(ipar.ndpth(c[0], 0) > 0);

        // 1.) number of faces per element
        for (size_t i = 0; i < num_elem; ++i) {
          d_geom_out << static_cast<int>(ipar.ncols(c[i])) << endl;
        }
        d_geom_out.flush();

        // 2.) number of nodes per face per element
        for (size_t i = 0; i < num_elem; ++i) {
          for (size_t j = 0; j < ipar.ncols(c[i]); j++)
            d_geom_out << static_cast<int>(ipar.ndpth(c[i], j)) << endl;
        }
        d_geom_out.flush();

        // 3.) local node indices per face per element
        for (size_t i = 0; i < num_elem; ++i) {
          for (size_t j = 0; j < ipar.ncols(c[i]); j++) {
            for (size_t k = 0; k < ipar.ndpth(c[i], j); ++k)
              d_geom_out << ens_vertex[ipar(c[i], j, k)];
            d_geom_out << endl;
          }
        }
        d_geom_out.flush();

      } else {

        // if not nface, only write the vertex list per cell
        for (size_t i = 0; i < num_elem; ++i) {
          Check(d_vrtx_cnt[type] > 0 ? static_cast<int>(ipar.ncols(c[i])) == d_vrtx_cnt[type]
                                     : true);
          for (size_t j = 0; j < ipar.ncols(c[i]); j++)
            d_geom_out << ens_vertex[ipar(c[i], j, 0)];
          d_geom_out << endl;
        }
        d_geom_out.flush();
      }
    }
  } // done looping over cell types
  d_geom_out.flush();
}

//------------------------------------------------------------------------------------------------//
//! Write out data to ensight vertex data.
template <typename FVF>
void Ensight_Translator::write_vrtx_data(const uint32_t part_num,
                                         const rtt_viz::Viz_Traits<FVF> &vrtx_data,
                                         const sf_int &vertices) {
  if (vrtx_data.nrows() == 0)
    return;

  size_t nvertices = vertices.size();
  size_t ndata = vrtx_data.ncols(0);

  std::string err = "Vertex data files not open. Must call open() before write_part().";
  Insist(d_vertex_out.size() == ndata, err.c_str());

  // loop over all vertex data fields and write out data for each field
  for (size_t nvd = 0; nvd < ndata; nvd++) {
    Ensight_Stream &vout = *d_vertex_out[nvd];

    if (rtt_c4::node() == 0) {
      vout << "part" << endl;
      vout << part_num << endl;
      vout << "coordinates" << endl;
    }
    vout.flush();

    for (size_t i = 0; i < nvertices; ++i)
      vout << vrtx_data(vertices[i], nvd) << endl;
    vout.flush();
  }
}

//------------------------------------------------------------------------------------------------//
//! Write out data to ensight cell data.
template <typename FVF>
void Ensight_Translator::write_cell_data(const uint32_t part_num,
                                         const rtt_viz::Viz_Traits<FVF> &cell_data,
                                         const sf2_int &cells_of_type) {
  if (cell_data.nrows() == 0)
    return;

  size_t const ndata = cell_data.ncols(0);

  std::string err = "Cell data files not open."
                    "  Must call open() before write_part().";

  Insist(d_cell_out.size() == static_cast<size_t>(ndata), err.c_str());

  // loop over all cell data fields and write out data for each field
  for (size_t ncd = 0; ncd < ndata; ncd++) {
    Ensight_Stream &cellout = *d_cell_out[ncd];

    if (rtt_c4::node() == 0) {
      cellout << "part" << endl;
      cellout << part_num << endl;
    }
    cellout.flush();

    // loop over ensight cell types
    for (unsigned type = 0; type < d_num_cell_types; type++) {
      const sf_int &c = cells_of_type[type];

      size_t num_elem = c.size();
      // if one processor has this type all must call the flush statements
      size_t global_num_elem = num_elem;

      if (d_decomposed)
        rtt_c4::global_sum(global_num_elem);

      // print out data if there are cells of this type
      if (global_num_elem > 0) {
        // printout cell-type name
        if (rtt_c4::node() == 0)
          cellout << d_cell_names[type] << endl;
        cellout.flush();

        // print out data
        for (size_t i = 0; i < num_elem; ++i)
          cellout << cell_data(c[i], ncd) << endl;
        cellout.flush();
      }
    }
  }
}

} // namespace rtt_viz

//------------------------------------------------------------------------------------------------//
// end of viz/Ensight_Translator.i.hh
//------------------------------------------------------------------------------------------------//
