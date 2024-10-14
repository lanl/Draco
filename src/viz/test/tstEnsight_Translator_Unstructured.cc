//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   viz/test/tstEnsight_Translator_Unstructured.cc
 * \author Thomas M. Evans, Ryan T. Wollaeger
 * \date   Wednesday, Oct 03, 2018, 15:27 pm
 * \brief  Ensight_Translator unstructured mesh test.
 * \note   Copyright (C) 2018-2024 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Assert.hh"
#include "dsxx/Release.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "dsxx/path.hh"
#include "c4/C4_Functions.hh"
#include "c4/ParallelUnitTest.hh"
#include "viz/Ensight_Translator.hh"
#include <climits>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using rtt_viz::Ensight_Translator;

//------------------------------------------------------------------------------------------------//
// ensight data for 2D unstructured mesh
template <typename IT>
void ensight_dump_test_unstr2d(rtt_dsxx::UnitTest &ut, std::string prefix, bool const binary,
                               bool const geom, bool const decomposed) {

  // short-cuts
  using vec_s = std::vector<std::string>;
  using vec_i = std::vector<IT>;
  using vec2_i = std::vector<vec_i>;
  using vec_d = std::vector<double>;
  using vec2_d = std::vector<vec_d>;

  if (binary)
    std::cout << "\nGenerating binary files...\n\n";
  else
    std::cout << "\nGenerating ascii files...\n\n";

  // >>> SET SCALAR ENSIGHT CTOR ARGS

  if (binary)
    prefix += "_binary";

  int const icycle = 1;
  double const time = .01;
  double const dt = .01;

  std::string const gd_wpath =
      rtt_dsxx::getFilenameComponent(ut.getTestInputPath(), rtt_dsxx::FC_NATIVE);

  // >>> INITIALIZE AND SET VECTOR DATA

  // dimensions
  size_t const ncells = 2;
  size_t const nvert = 8;
  size_t const ndim = 2;
  size_t const ndata = 2;
  std::vector<size_t> nvert_per_cell = {4, 6};

  // size the cell-vertex vector
  vec2_i ipar = {vec_i(nvert_per_cell[0]), vec_i(nvert_per_cell[1])};
  Check(ipar.size() == ncells);
  vec2_d vrtx_data(nvert, vec_d(ndata));
  vec2_d cell_data(ncells, vec_d(ndata));
  vec2_d pt_coor(nvert, vec_d(ndim));
  // set the element type to be unstructured
  vec_i const iel_type(ncells, rtt_viz::unstructured);
  // set the vertex-centered and cell-centered data names
  vec_s const vdata_names = {"Densities", "Temperatures"};
  Check(vdata_names.size() == ndata);
  vec_s const cdata_names = {"Velocity", "Pressure"};
  Check(cdata_names.size() == ndata);
  // region information
  vec_i const rgn_index = {1, 1};
  Check(rgn_index.size() == ncells);
  vec_s const rgn_name = {"RGN_A"};
  vec_i const rgn_data = {1};

  // create some arbitrary cell and vertex based data (as in structured test)
  for (size_t i = 0; i < ndata; i++) {
    // cell data
    for (size_t cell = 0; cell < ncells; cell++) {
      Check(1 + cell < INT_MAX);
      cell_data[cell][i] = static_cast<int>(1 + cell);
    }

    // vrtx data
    for (size_t v = 0; v < nvert; v++) {
      Check(1 + v < INT_MAX);
      vrtx_data[v][i] = static_cast<int>(1 + v);
    }
  }

  // Build path for the input file "cell_data_unstr2d"
  std::string const cdInputFile = ut.getTestSourcePath() + std::string("cell_data_unstr2d");
  std::ifstream input(cdInputFile.c_str());
  if (!input)
    ITFAILS;

  for (auto &vec1d : pt_coor)
    for (auto &elem : vec1d)
      input >> elem;
  for (size_t i = 0; i < ipar.size(); i++) {
    for (size_t j = 0; j < ipar[i].size(); j++)
      input >> ipar[i][j];
  }

  // build an Ensight_Translator (make sure it overwrites any existing stuff)
  Ensight_Translator translator(prefix, gd_wpath, vdata_names, cdata_names, true, geom, binary,
                                decomposed);

  translator.ensight_dump(icycle, time, dt, ipar, iel_type, rgn_index, pt_coor, vrtx_data,
                          cell_data, rgn_data, rgn_name);

  std::vector<double> dump_times = translator.get_dump_times();
  if (dump_times.size() != 1)
    ITFAILS;
  if (!rtt_dsxx::soft_equiv(dump_times[0], 0.01))
    ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("ensight_dump_test_unstr2d finished successfully.");
  else
    FAILMSG("ensight_dump_test_unstr2d did not finish successfully.");
  return;
}

//------------------------------------------------------------------------------------------------//
// ensight data for 3D unstructured mesh
template <typename IT>
void ensight_dump_test_unstr3d(rtt_dsxx::UnitTest &ut, std::string prefix, bool const binary,
                               bool const geom, bool const decomposed) {

  // short-cuts
  using vec_s = std::vector<std::string>;
  using vec_i = std::vector<IT>;
  using vec2_i = std::vector<vec_i>;
  using vec3_i = std::vector<vec2_i>;
  using vec_d = std::vector<double>;
  using vec2_d = std::vector<vec_d>;

  if (binary)
    std::cout << "\nGenerating binary files...\n\n";
  else
    std::cout << "\nGenerating ascii files...\n\n";

  // >>> SET SCALAR ENSIGHT CTOR ARGS

  if (binary)
    prefix += "_binary";

  int const icycle = 1;
  double const time = .01;
  double const dt = .01;

  std::string const gd_wpath =
      rtt_dsxx::getFilenameComponent(ut.getTestInputPath(), rtt_dsxx::FC_NATIVE);

  // >>> INITIALIZE AND SET VECTOR DATA

  // dimensions
  size_t const ncells = 3;
  size_t const nvert = 14;
  size_t const ndim = 3;
  size_t const ndata = 2;
  std::vector<size_t> nface_per_cell = {5, 5, 7};
  std::vector<std::vector<size_t>> nvert_per_face_per_cell = {
      {3, 3, 4, 4, 4}, {3, 3, 4, 4, 4}, {5, 5, 4, 4, 4, 4, 4}};

  // size the cell-vertex vector
  vec3_i ipar(3);
  for (size_t i = 0; i < ncells; ++i) {
    ipar[i].resize(nface_per_cell[i]);
    for (size_t j = 0; j < nface_per_cell[i]; ++j) {
      ipar[i][j].resize(nvert_per_face_per_cell[i][j]);
    }
  }
  Check(ipar.size() == ncells);
  vec2_d vrtx_data(nvert, vec_d(ndata));
  vec2_d cell_data(ncells, vec_d(ndata));
  vec2_d pt_coor(nvert, vec_d(ndim));
  // set the element type to be unstructured
  vec_i const iel_type(ncells, rtt_viz::unstructured_3d);
  // set the vertex-centered and cell-centered data names
  vec_s const vdata_names = {"Densities", "Temperatures"};
  Check(vdata_names.size() == ndata);
  vec_s const cdata_names = {"Velocity", "Pressure"};
  Check(cdata_names.size() == ndata);
  // region information
  vec_i const rgn_index = {1, 1, 1};
  vec_s const rgn_name = {"RGN_A"};
  vec_i const rgn_data = {1};

  // create some arbitrary cell and vertex based data (as in structured test)
  for (size_t i = 0; i < ndata; i++) {
    // cell data
    for (size_t cell = 0; cell < ncells; cell++) {
      Check(1 + cell < INT_MAX);
      cell_data[cell][i] = static_cast<int>(1 + cell);
    }

    // vrtx data
    for (size_t v = 0; v < nvert; v++) {
      Check(1 + v < INT_MAX);
      vrtx_data[v][i] = static_cast<int>(1 + v);
    }
  }

  // Build path for the input file "cell_data_unstr3d"
  // \todo: create this input file
  std::string const cdInputFile = ut.getTestSourcePath() + std::string("cell_data_unstr3d");
  std::ifstream input(cdInputFile.c_str());
  if (!input)
    ITFAILS;

  for (auto &vec1d : pt_coor)
    for (auto &elem : vec1d)
      input >> elem;
  for (size_t i = 0; i < ipar.size(); i++) {
    for (size_t j = 0; j < ipar[i].size(); j++) {
      for (size_t k = 0; k < ipar[i][j].size(); k++)
        input >> ipar[i][j][k];
    }
  }

  // build an Ensight_Translator (make sure it overwrites any existing stuff)
  Ensight_Translator translator(prefix, gd_wpath, vdata_names, cdata_names, true, geom, binary,
                                decomposed);

  translator.ensight_dump(icycle, time, dt, ipar, iel_type, rgn_index, pt_coor, vrtx_data,
                          cell_data, rgn_data, rgn_name);

  std::vector<double> dump_times = translator.get_dump_times();
  if (dump_times.size() != 1)
    ITFAILS;
  if (!rtt_dsxx::soft_equiv(dump_times[0], 0.01))
    ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("ensight_dump_test_unstr3d finished successfully.");
  else
    FAILMSG("ensight_dump_test_unstr3d did not finish successfully.");
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // ASCII dumps
    bool binary(false);
    bool const geom{false};
    if (rtt_c4::node() == 0) {

      // check serial writes
      bool const decomposed{false};

      // 2D
      std::string prefix = "unstr2d_testproblem_serial_" + std::to_string(rtt_c4::nodes());

      // Binary dumps
      binary = true;
      ensight_dump_test_unstr2d<uint32_t>(ut, prefix, binary, geom, decomposed);

      // ASCII dumps with unsigned integer data
      binary = false;
      ensight_dump_test_unstr2d<uint32_t>(ut, prefix, binary, geom, decomposed);

      // 3D
      prefix = "unstr3d_testproblem_serial_" + std::to_string(rtt_c4::nodes());

      // Binary dumps
      binary = true;
      ensight_dump_test_unstr3d<uint32_t>(ut, prefix, binary, geom, decomposed);

      // ASCII dumps with unsigned integer data
      binary = false;
      ensight_dump_test_unstr3d<uint32_t>(ut, prefix, binary, geom, decomposed);
    }

    // check decomposed writes
    bool const decomposed{true};

    // 2D
    std::string prefix = "unstr2d_testproblem_parallel_" + std::to_string(rtt_c4::nodes());

    // Binary dumps
    binary = true;
    ensight_dump_test_unstr2d<uint32_t>(ut, prefix, binary, geom, decomposed);

    // ASCII dumps with unsigned integer data
    binary = false;
    ensight_dump_test_unstr2d<uint32_t>(ut, prefix, binary, geom, decomposed);

    // 3D
    prefix = "unstr3d_testproblem_parallel_" + std::to_string(rtt_c4::nodes());

    // Binary dumps
    binary = true;
    ensight_dump_test_unstr3d<uint32_t>(ut, prefix, binary, geom, decomposed);

    // ASCII dumps with unsigned integer data
    binary = false;
    ensight_dump_test_unstr3d<uint32_t>(ut, prefix, binary, geom, decomposed);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of viz/test/tstEnsight_Translator_Unstructured.cc
//------------------------------------------------------------------------------------------------//
