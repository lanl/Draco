//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   VendorChecks/test/tstParmetis.cc
 * \date   Monday, May 16, 2016, 16:30 pm
 * \brief  Attempt to link to libparmetis and run a simple problem.
 * \note   Copyright (C) 2016-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include <parmetis.h>
#include <sstream>
#include <vector>

// Example borrowed from the Parmetis manual.

void test_parmetis(rtt_c4::ParallelUnitTest &ut) {
  using std::cout;
  using std::endl;

  Insist(rtt_c4::nodes() == 3, "test_parmetis must be called with 3 MPI ranks exactly.");

  // MPI VARIABLES
  size_t const MPI_PROC_ID = rtt_c4::node();
  size_t const MPI_PROC_TOTAL_NUM = rtt_c4::nodes();

  if (MPI_PROC_ID == 0)
    cout << " Parmetis example from LiberLocus." << '\n';
  cout << "I am Proc " << MPI_PROC_ID << '\n';

  // Needed by parmetis

  // These store the (local) adjacency structure of the graph at each processor.
  std::vector<idx_t> xadj(6);
  std::vector<idx_t> adjncy;
  // These store the weights of the vertices and edges.
  idx_t *vwgt = nullptr;
  idx_t *adjwgt = nullptr;
  // This is used to indicate if the graph is weighted. (0 == no weights)
  idx_t wgtflag = 0;
  // C-style numbering that starts from 0.
  idx_t numflag = 0;
  // This is used to specify the number of weights that each vertex has. It is also the number of
  // balance constraints that must be satisfied.
  idx_t ncon = 1;
  // This is used to specify the number of sub-domains that are desired. Note that the number of
  // subdomains is independent of the number of processors that call this routine.
  idx_t nparts = 3;
  constexpr auto one = static_cast<real_t>(1);
  // An array of size ncon * that is used to specify the fraction of vertex weight that
  // should be distributed to each sub-domain for each balance constraint. If all of the sub-domains
  // are to be of the same size for every vertex weight, then each of the ncon * elements
  // should be set to a value of 1/nparts.
  std::vector<real_t> tpwgts(ncon * nparts, (one / static_cast<real_t>(nparts)));
  // An array of size ncon that is used to specify the imbalance tolerance for each vertex weight,
  // with 1 being perfect balance and nparts being perfect imbalance. A value of 1.05 for each of
  // the ncon weights is recommended.
#if REALTYPEWIDTH == 64
  real_t ubvec = 1.05;
#else
  real_t ubvec = 1.05F;
#endif
  // This is an array of integers that is used to pass additional parameters for the routine.
  std::vector<idx_t> options(4, 0);
  // Upon successful completion, the number of edges that are cut by the partitioning is written to
  // this parameter.
  idx_t edgecut(0);

  MPI_Comm_dup(MPI_COMM_WORLD, &rtt_c4::communicator);

  // This is an array of size equal to the number of locally-stored vertices. Upon successful
  // completion the partition vector of the locally-stored vertices is written to this array.
  Check(MPI_PROC_ID < INT_MAX);
  std::vector<idx_t> part(5, static_cast<idx_t>(MPI_PROC_ID));

  // This array describes how the vertices of the graph are distributed among the processors. Its
  // contents are identical for every processor.
  std::vector<idx_t> vtxdist = {0, 5, 10, 15};

  // Dependent on each processor
  if (MPI_PROC_ID == 0) {
    adjncy.resize(13);
    xadj = {0, 2, 5, 8, 11, 13};
    adjncy = {1, 5, 0, 2, 6, 1, 3, 7, 2, 4, 8, 3, 9};
  } else if (MPI_PROC_ID == 1) {
    adjncy.resize(18);
    xadj = {0, 3, 7, 11, 15, 18};
    adjncy = {0, 6, 10, 1, 5, 7, 11, 2, 6, 8, 12, 3, 7, 9, 13, 4, 8, 14};
  } else if (MPI_PROC_ID == 2) {
    adjncy.resize(13);
    xadj = {0, 2, 5, 8, 11, 13};
    adjncy = {5, 11, 6, 10, 12, 7, 11, 13, 8, 12, 14, 9, 13};
  }
  if (MPI_PROC_ID == 0)
    cout << "parmetis initialized." << '\n';

  int result = ParMETIS_V3_PartKway(&vtxdist[0], &xadj[0], &adjncy[0], vwgt, adjwgt, &wgtflag,
                                    &numflag, &ncon, &nparts, &tpwgts[0], &ubvec, &options[0],
                                    &edgecut, &part[0], &rtt_c4::communicator);

  if (result == METIS_OK) {
    std::ostringstream msg;
    msg << "[" << MPI_PROC_ID << "] ParMETIS_V3_AdaptiveRepart did not return an error.";
    PASSMSG(msg.str());
  } else {
    std::ostringstream msg;
    msg << "[" << MPI_PROC_ID << "] ParMETIS_V3_AdaptiveRepart returned an error code.";
    FAILMSG(msg.str());
  }

  if (MPI_PROC_ID == 0)
    cout << "parmetis finalized." << endl;

  for (size_t pid = 0; pid < MPI_PROC_TOTAL_NUM; ++pid) {
    rtt_c4::global_barrier();
    cout << MPI_PROC_ID << " edgecut " << edgecut << '\n';
    for (int i = 0; i < 5; i++)
      cout << "[" << MPI_PROC_ID << "] " << part[i] << endl;
  }

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_parmetis(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstParmetis.cc
//------------------------------------------------------------------------------------------------//
