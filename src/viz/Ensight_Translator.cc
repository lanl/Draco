//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   viz/Ensight_Translator.cc
 * \author Thomas M. Evans
 * \brief  Ensight_Translator implementation file (non-templated code).
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Ensight_Translator.hh"
#include "dsxx/Assert.hh"
#include "dsxx/SystemCall.hh"
#include "dsxx/path.hh"
#include "c4/C4_Functions.hh"
#include "viz/Ensight_Stream.hh"
#include <algorithm>
#include <cerrno>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>

namespace rtt_viz {

using std::ios;
using std::ofstream;
using std::setiosflags;
using std::setw;
using std::string;

//------------------------------------------------------------------------------------------------//
// Local Helper Functions
//------------------------------------------------------------------------------------------------//

namespace {

constexpr std::string_view forbidden_chars = R"()[]+-@!#*^$/ )";
inline bool name_allowed(const std::string &input) noexcept {
  return std::ranges::find_first_of(input, forbidden_chars) == input.end();
}

} // namespace

//------------------------------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//------------------------------------------------------------------------------------------------//

// NOLINTBEGIN(performance-avoid-endl)
// * Using '\n' instead of std::endl causes the unit tests to fail.

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Opens the geometry and variable files.
 *
 * \param[in] icycle Cycle number for this dump.
 * \param[in] time   Time value for this dump.
 * \param[in] dt     Timestep at this dump.  This parameter is only used for diagnostic and is not
 *                   placed in the Ensight dump.
 */
void Ensight_Translator::open(const int icycle, const double time, const double dt) {
  Insist(!d_geom_out.is_open(), "Attempted to open an already open geometry file!");

  using std::ostringstream;
  using std::string;

  // Increment local dump counter and add dump time
  d_dump_times.push_back(time);
  size_t const igrdump_num = d_dump_times.size();
  Check(igrdump_num < 10000);

  // create ensight postfix indicators
  ostringstream postfix_build;
  ostringstream post_number;

  if (igrdump_num < 10)
    post_number << "000" << igrdump_num;
  else if (igrdump_num < 100)
    post_number << "00" << igrdump_num;
  else if (igrdump_num < 1000)
    post_number << "0" << igrdump_num;
  else
    post_number << igrdump_num;

  postfix_build << "data." << post_number.str();
  string const postfix = postfix_build.str();

  // announce the graphics dump
  if (rtt_c4::node() == 0) {
    std::cout << ">>> ENSIGHT GRAPHICS DUMP: icycle= " << icycle << " time= " << time
              << " dt= " << dt << "\ndir= " << d_prefix << ", dump_number= " << igrdump_num << '\n';
    // write case file
    write_case();
  }

  // >>> Open the geometry file.
  if ((!d_static_geom) || (d_dump_times.size() == 1)) {
    // make output file for this timestep
    string filename = d_geo_dir + "/";
    if (d_static_geom)
      filename += "data";
    else
      filename += postfix;

    const bool geom{true};
    d_geom_out.open(filename, d_binary, geom, d_decomposed);

    // write the header
    if (rtt_c4::node() == 0) {
      d_geom_out << "Description line 1" << endl;

      ostringstream s;
      s << "probtime " << time << " cycleno " << icycle;
      d_geom_out << s.str() << endl;
      d_geom_out << "node id given" << endl;
      d_geom_out << "element id given" << endl;
    }

    d_geom_out.flush();
  }

  // >>> Open the vertex data files.
  d_vertex_out.resize(d_vdata_names.size());

  // loop over all vertex data fields and write out data for each field
  for (size_t nvd = 0; nvd < d_vdata_names.size(); nvd++) {
    // open file for this data
    std::string const filename = d_vdata_dirs[nvd] + "/" + postfix;
    const bool geom{false};
    d_vertex_out[nvd] = std::make_unique<Ensight_Stream>(filename, d_binary, geom, d_decomposed);

    if (rtt_c4::node() == 0) {
      *d_vertex_out[nvd] << d_vdata_names[nvd] << endl;
    }
    d_vertex_out[nvd]->flush();
  }

  // >>> Open the cell data files.
  d_cell_out.resize(d_cdata_names.size());

  // loop over all cell data fields
  for (size_t ncd = 0; ncd < d_cdata_names.size(); ncd++) {
    // open file for this data
    std::string const filename = d_cdata_dirs[ncd] + "/" + postfix;
    const bool geom{false};
    d_cell_out[ncd] = std::make_unique<Ensight_Stream>(filename, d_binary, geom, d_decomposed);

    if (rtt_c4::node() == 0) {
      *d_cell_out[ncd] << d_cdata_names[ncd] << endl;
    }
    d_cell_out[ncd]->flush();
  }
}
// NOLINTEND(performance-avoid-endl)

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Closes any open file streams.
 *
 * Calling this function is unnecessary if this object is destroyed.
 */
void Ensight_Translator::close() {
  if (d_geom_out.is_open())
    d_geom_out.close();

  for (auto &vo : d_vertex_out)
    if (vo->is_open())
      vo->close();

  for (auto &co : d_cell_out)
    if (co->is_open())
      co->close();
}

//------------------------------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Creates some of the file prefixes and filenames for ensight dump.
 *
 * \param[in] prefix std::string giving the name of the problem
 */
void Ensight_Translator::create_filenames(const std::string &prefix) {
  if (d_dump_dir[d_dump_dir.size() - 1] == rtt_dsxx::UnixDirSep ||
      d_dump_dir[d_dump_dir.size() - 1] == rtt_dsxx::WinDirSep) {
    // ensight directory name
    d_prefix = d_dump_dir + prefix + "_ensight";
  } else {
    // ensight directory name
    d_prefix = d_dump_dir + rtt_dsxx::dirSep + prefix + "_ensight";
  }
  // case file name
  d_case_filename = d_prefix + rtt_dsxx::dirSep + prefix + ".case";
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Common initializer for constructors.
 *
 * \param[in] graphics_continue If true, use existing ensight directory. If false, create or wipe
 *               out the existing directory.
 */
void Ensight_Translator::initialize(const bool graphics_continue) {
  using std::cerr;
  using std::endl;
  using std::strerror;

  d_num_cell_types = 17;

  // Assign values to d_cell_names. These are the official "Ensight" names that must be used in the
  // Ensight file.
  d_cell_names = {"point",  "bar2",   "bar3",    "tria3",    "tria6",     "quad4",
                  "quad8",  "tetra4", "tetra10", "pyramid5", "pyramid13", "hexa8",
                  "hexa20", "penta6", "penta15", "nsided",   "nfaced"};
  Check(d_cell_names.size() == d_num_cell_types);

  // Assign values to vrtx_count, the number of vertices in a cell.
  d_vrtx_cnt = {1, 2, 3, 3, 6, 4, 8, 4, 10, 5, 13, 8, 20, 6, 15, -1, -1};
  Check(d_vrtx_cnt.size() == d_num_cell_types);

  // Assign values to d_cell_type_index. The user will use these to identify cell types.
  d_cell_type_index = {point,
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
                       unstructured_3d};
  Check(d_cell_type_index.size() == d_num_cell_types);

  if (rtt_c4::node() == 0) {
    // Check d_dump_dir
    rtt_dsxx::draco_getstat dumpDirStat(d_dump_dir);
    if (!dumpDirStat.isdir()) {
      std::ostringstream dir_error;
      dir_error << "Error opening dump directory \"" << d_dump_dir << "\": " << strerror(errno);
      Insist(dumpDirStat.isdir(), dir_error.str());
    }

    // try to create d_prefix
    rtt_dsxx::draco_mkdir(d_prefix);
    rtt_dsxx::draco_getstat prefixDirStat(d_prefix);
    if (!prefixDirStat.isdir()) {
      std::ostringstream dir_error;
      dir_error << "Unable to create EnSight directory \"" << d_prefix << "\": " << strerror(errno);
      Insist(dumpDirStat.isdir(), dir_error.str());
    }

    // NOLINTBEGIN(misc-include-cleaner)
    // * Don't worry about including system headers needed for direct inclusion of 'stat' def.

    // See if the case file exists
    struct stat sbuf {};
    int const stat_ret = stat(d_case_filename.c_str(), &sbuf);
    // NOLINTEND(misc-include-cleaner)

    // build the ensight directory if this is not a continuation
    if (!graphics_continue) {
      // We have guaranteed that our prefix directory exists at this point.  Now, wipe out files
      // that we might have created in there...
      if (!stat_ret) {
        rtt_dsxx::draco_remove_dir(d_prefix);
        rtt_dsxx::draco_mkdir(d_prefix);
      }
    } else {
      // We were asked for a continuation.  Complain if we don't have a case file.
      if (stat_ret) {
        std::ostringstream dir_error;
        dir_error << "EnSight directory \"" << d_prefix << "\" doesn't contain a case file!";
        Insist(0, dir_error.str().c_str());
      }
    }
  }

  // Check to make sure the variable names are of acceptable length and contain no forbidden
  // characters. Ensight prohibits the characters "( ) [ ] + - @ !  # * ^ $ / space", and requires
  // the names be 19 characters or less. Moreover, since these names will be used to label output
  // and to create directories, the names should also be unique.

  size_t const nvdata = d_vdata_names.size();
  size_t const ncdata = d_cdata_names.size();

  // Create a name list for testing.
  sf_string name_tmp(nvdata + ncdata);
  for (size_t i = 0; i < nvdata; i++)
    name_tmp[i] = d_vdata_names[i];
  for (size_t i = 0; i < ncdata; i++)
    name_tmp[i + nvdata] = d_cdata_names[i];
  // Check for name lengths out of limits:
  for (std::string const &iname : name_tmp) {
    Insist(!iname.empty() && iname.size() < 20,
           std::string("Error: Size of name is not in allowable range (1-20 chars), name = \"") +
               iname + "\".");
  }
  // Check for bad characters.
  for (std::string const &iname : name_tmp) {
    Insist(name_allowed(iname), std::string("Error: Found disallowed character(s) in name \"") +
                                    iname + "\"." + "Forbidden characters: \"" +
                                    std::string(forbidden_chars) + std::string("\""));
  }
  // Check for non-unique names (generating the message is a bit messy.)
  Insist(std::unordered_set<std::string>(name_tmp.begin(), name_tmp.end()).size() ==
             name_tmp.size(),
         std::string("Error: Duplicate name found, names = {") +
             std::accumulate(
                 std::next(name_tmp.begin()), name_tmp.end(), name_tmp[0],
                 [](std::string a, std::string const &b) { return std::move(a) + "," + b; }) +
             "}.");

  // calculate and make the geometry directory if this is not a continuation
  d_geo_dir = d_prefix + "/geo";
  if (!graphics_continue && rtt_c4::node() == 0)
    rtt_dsxx::draco_mkdir(d_geo_dir);

  // make data directory names and directories
  d_vdata_dirs.resize(d_vdata_names.size());
  d_cdata_dirs.resize(d_cdata_names.size());
  for (size_t i = 0; i < d_vdata_names.size(); i++) {
    d_vdata_dirs[i] = d_prefix + rtt_dsxx::dirSep + d_vdata_names[i];

    // if this is not a continuation make the directory
    if (!graphics_continue && rtt_c4::node() == 0)
      rtt_dsxx::draco_mkdir(d_vdata_dirs[i]);
  }
  for (size_t i = 0; i < d_cdata_names.size(); i++) {
    d_cdata_dirs[i] = d_prefix + rtt_dsxx::dirSep + d_cdata_names[i];

    // if this is not a continuation make the directory (Mat_Erg, Mat_Temp, Rad_Temp, etc.)
    if (!graphics_continue && rtt_c4::node() == 0)
      rtt_dsxx::draco_mkdir(d_cdata_dirs[i]);
  }
}

//------------------------------------------------------------------------------------------------//
//! Write out case file.
void Ensight_Translator::write_case() {
  // create the case file name (directory already created)
  const char *filename = d_case_filename.c_str();

  // open the case file
  ofstream caseout(filename);

  // write the format header
  caseout << "FORMAT\ntype: ensight gold\n\n";

  // write the geometry file block
  caseout << "GEOMETRY\n";
  if (d_static_geom)
    caseout << "model: ./geo/data";
  else
    caseout << "model: 1   ./geo/data.****";

  // write the variable block header
  caseout << "\n\nVARIABLE\n";

  // write the pointer to the node variables
  for (auto &data_name : d_vdata_names)
    caseout << "scalar per node:    1  " << setw(19) << setiosflags(ios::left) << data_name
            << setw(3) << " ./" << data_name << "/data.****\n";

  // write the pointer to the cell variables
  for (auto &cdata_name : d_cdata_names)
    caseout << "scalar per element: 1  " << setw(19) << setiosflags(ios::left) << cdata_name
            << setw(3) << " ./" << cdata_name << "/data.****\n";

  // write out the time block
  caseout << "\nTIME\n"
          << "time set:              " << setw(4) << "   1\n"
          << "number of steps:       " << setw(4) << setiosflags(ios::right) << d_dump_times.size()
          << "\nfilename start number: " << setw(4) << "   1\n"
          << "filename increment:    " << setw(4) << "   1\n"
          << "time values:           \n";

  // write out times
  caseout.precision(5);
  caseout.setf(ios::scientific, ios::floatfield);
  for (double const dump_time : d_dump_times)
    caseout << setw(12) << setiosflags(ios::right) << dump_time << '\n';
}

} // namespace rtt_viz

//------------------------------------------------------------------------------------------------//
// end of Ensight_Translator.cc
//------------------------------------------------------------------------------------------------//
