//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstProcessor_Group.cc
 * \author Kent Budge
 * \date   Tue Sep 21 11:45:44 2004
 * \brief  Unit tests for the Processor_Group class and member functions.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "c4/Processor_Group.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstProcessor_Group(rtt_dsxx::UnitTest &ut) {

  unsigned const pid = rtt_c4::node();

#ifdef C4_MPI
  // Test construction
  Processor_Group comm(2);
#else
  // Test construction
  Processor_Group comm(1);
#endif // C4_MPI

  std::string pm = "Processor_Group constructed on pid " + to_string(pid) + ".";
  PASSMSG(pm);

  // Test sum
  size_t const group_pids = comm.size();

  size_t const base = pid % 2;
  vector<double> sum(1, static_cast<double>(base + 1));
  comm.sum(sum);
  if (rtt_dsxx::soft_equiv(sum[0], static_cast<double>(group_pids * (base + 1))))
    PASSMSG("Correct processor group sum");
  else
    FAILMSG("NOT correct processor group sum");

  // Test assemble_vector
  {
    vector<double> myvec;
    size_t const vlen(5);
    for (size_t i = 0; i < vlen; ++i)
      myvec.push_back(static_cast<double>(pid * 1000 + i));
    vector<double> globalvec;
    comm.assemble_vector(myvec, globalvec);

    if (globalvec.size() != group_pids * vlen)
      ITFAILS;

    // Check the the first 5 elements
    vector<double> goldglobalvec;
    for (size_t j = 0; j < group_pids; ++j)
      for (size_t i = 0; i < vlen; ++i)
        goldglobalvec.push_back(static_cast<double>((base + 2 * j) * 1000 + i));

    if (!rtt_dsxx::soft_equiv(goldglobalvec.begin(), goldglobalvec.end(), globalvec.begin(),
                              globalvec.end()))
      ITFAILS;
  }

  // Test assemble_vector - pointer signature
  {
    vector<double> myvec;
    size_t const vlen(5);
    for (size_t i = 0; i < vlen; ++i)
      myvec.push_back(static_cast<double>(pid * 1000 + i));
    vector<double> globalvec(group_pids * vlen);
    Check(myvec.size() < UINT32_MAX);
    comm.assemble_vector(&myvec[0], &globalvec[0], static_cast<unsigned>(myvec.size()));

    FAIL_IF_NOT(globalvec.size() == group_pids * vlen);

    // Check the the first 5 elements
    vector<double> goldglobalvec;
    for (size_t j = 0; j < group_pids; ++j)
      for (size_t i = 0; i < vlen; ++i)
        goldglobalvec.push_back(static_cast<double>((base + 2 * j) * 1000 + i));

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(goldglobalvec.begin(), goldglobalvec.end(), globalvec.begin(),
                                     globalvec.end()));
  }

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {

  ParallelUnitTest ut(argc, argv, release);
  try {
    tstProcessor_Group(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstProcessor_Group.cc
//------------------------------------------------------------------------------------------------//
