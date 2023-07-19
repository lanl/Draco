//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   viz/test/tstvtk_camr_writer.cc
 * \author Ben R. Ryan
 * \date   2023/7/12
 * \brief  vtk_camr_writer_t test.
 * \note   Copyright (C) 2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "viz/vtk_camr_writer.hh"

using namespace std;
using rtt_viz::vtk_camr_writer_t;

//------------------------------------------------------------------------------------------------//
// vtk_camr_writer I/O test
void test_vtk_camr_writer(rtt_dsxx::UnitTest &ut) {

  vtk_camr_writer_t writer;

  // Simple 2x2 mesh

  const int dim = 2;
  const int num_cells = 4;
  const double dx = 0.5;
  const std::vector<double> vertices{0., 0., 0.5, 0., 0., 0.5, 0.5, 0.5};
  const std::vector<int> levels{1, 1, 1, 1};
  const std::vector<double> cell_data{1., 2., 3., 4.};
  const std::string varname = "variable";

  FAIL_IF_NOT(writer.open("vtk_out.vtk") == 0);
  FAIL_IF_NOT(writer.write_camr_unstructured_begin() == 0);
  FAIL_IF_NOT(writer.write_camr_piece_begin(dim, num_cells) == 0);
  FAIL_IF_NOT(writer.write_camr_points(dim, num_cells, dx, vertices, levels) == 0);
  FAIL_IF_NOT(writer.write_camr_cells(dim, num_cells) == 0);
  FAIL_IF_NOT(writer.write_camr_cell_data_begin() == 0);
  FAIL_IF_NOT(writer.write_camr_cell_data(varname, cell_data) == 0);
  FAIL_IF_NOT(writer.write_camr_cell_data_end() == 0);
  FAIL_IF_NOT(writer.write_camr_piece_end() == 0);
  FAIL_IF_NOT(writer.write_camr_unstructured_end() == 0);
  FAIL_IF_NOT(writer.close() == 0);

  if (ut.numFails == 0) {
    PASSMSG("test_vtk_camr_writer passes.");
  }

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // >>> UNIT TESTS
    test_vtk_camr_writer(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstvtk_camr_writer.cc
//------------------------------------------------------------------------------------------------//
