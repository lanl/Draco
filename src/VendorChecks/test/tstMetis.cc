//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   VendorChecks/test/tstMetis.cc
 * \date   Wednesday, May 11, 2016, 12:01 pm
 * \brief  Attempt to link to libmetis and run a simple problem.
 * \note   Copyright (C) 2016-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifdef _MSC_FULL_VER
// - 4005: macro redefinition - see https://github.com/microsoft/vcpkg/issues/13774
#pragma warning(push)
#pragma warning(disable : 4005)
#endif

#include "metis.h"

#ifdef _MSC_FULL_VER
#pragma warning(pop)
#endif

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include <array>
#include <vector>

// Original provided by Erik Zenker, https://gist.github.com/erikzenker/c4dc42c8d5a8c1cd3e5a

void test_metis(rtt_dsxx::UnitTest &ut) {
  idx_t nVertices = 10;
  idx_t nWeights = 1;
  idx_t nParts = 2;

  idx_t objval(0);
  std::vector<idx_t> part(nVertices, 0);

  // here's the mesh, there is only one valid cut so the expected result (or a mirror of it) should
  // alays be obtained
  //
  //  0 \       / 6
  //  1 \       / 7
  //  2 - 4 - 5 - 8
  //  3 /       \ 9

  // Indexes of starting points in adjacent array
  std::array<idx_t, 11> xadj = {0, 1, 2, 3, 4, 9, 14, 15, 16, 17, 18};

  // Adjacent vertices in consecutive index order
  // conn. for:     0, 1, 2, 3, 4            , 5,            ,6, 7, 8, 9
  // index:         0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13, 14,15,16,17
  std::array<idx_t, 18> adjncy = {4, 4, 4, 4, 0, 1, 2, 3, 5, 4, 6, 7, 8, 9, 5, 5, 5, 5};

  // Weights of vertices: if all weights are equal then can be set to NULL
  std::vector<idx_t> vwgt(nVertices * nWeights, 0);

  // Partition a graph into k parts using either multilevel recursive bisection or multilevel k-way
  // partitioning.
  int ret = METIS_PartGraphKway(&nVertices, &nWeights, xadj.data(), adjncy.data(), nullptr, nullptr,
                                nullptr, &nParts, nullptr, nullptr, nullptr, &objval, &part[0]);

  std::cout << "partition: ";
  for (int32_t i = 0; i < nVertices; ++i) {
    std::cout << part[i] << " ";
  }
  std::cout << std::endl;

  if (ret == METIS_OK)
    PASSMSG("Successfully called METIS_PartGraphKway().");
  else
    FAILMSG("Call to METIS_PartGraphKway() failed.");

  std::array<int, 10> expectedResult = {1, 1, 1, 1, 1, 0, 0, 0, 0, 0};
  std::array<int, 10> mirrorExpectedResult = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
  std::vector<idx_t> vExpectedResult(expectedResult.begin(), expectedResult.end());
  std::vector<idx_t> vMirrorExpectedResult(mirrorExpectedResult.begin(),
                                           mirrorExpectedResult.end());
  if (part == vExpectedResult || part == vMirrorExpectedResult)
    PASSMSG("Metis returned the expected result.");
  else
    FAILMSG("Metis failed to return the expected result.");

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_metis(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstMetis.cc
//------------------------------------------------------------------------------------------------//
