//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstGatherScatter.cc
 * \author Kent Budge
 * \date   Wed Apr 28 09:31:51 2010
 * \brief  Test c4::gather and c4::scatter functions
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "c4/gatherv.hh"
#include "c4/scatterv.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include <numeric>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstAllgather(UnitTest &ut) {
  int mypid = rtt_c4::node();
  vector<int> allpid(rtt_c4::nodes());
  int status = allgather(&mypid, &allpid[0], 1);
  if (status == 0) {
    PASSMSG("status is okay from allgather");
    status = 0;
    for (int i = 0; i < rtt_c4::nodes(); ++i) {
      if (allpid[i] != i) {
        status = 1;
        break;
      }
    }
    if (status == 0) {
      PASSMSG("gathered values are okay for allgather");
    } else {
      FAILMSG("gathered values are NOT okay for allgather");
    }
  } else
    FAILMSG("status is NOT okay from allgather");
}

//------------------------------------------------------------------------------------------------//
void tstDeterminateGatherScatter(UnitTest &ut) {
  unsigned pid = node();
  unsigned const number_of_processors = nodes();
  vector<unsigned> pids(number_of_processors);
  gather(&pid, &pids[0], 1);

  if (pid == 0)
    for (unsigned i = 1; i < number_of_processors; ++i)
      pids[i] += pids[i - 1];

  unsigned base;
  scatter(&pids[0], &base, 1);

  if (base == pid * (pid + 1) / 2)
    PASSMSG("correct base summation in gather/scatter");
  else
    FAILMSG("NOT correct base summation in gather/scatter");

  {
    vector<char> alphabet(number_of_processors);
    std::iota(alphabet.begin(), alphabet.end(), 'A');

    char c = alphabet[pid];
    vector<char> vc(number_of_processors);
    gather(&c, &vc[0], 1);

    if (pid == 0)
      FAIL_IF_NOT(std::equal(vc.begin(), vc.end(), alphabet.begin(), alphabet.end()));
  }

  return;
}

//------------------------------------------------------------------------------------------------//
void tstIndeterminateGatherScatterv(UnitTest &ut) {
  unsigned const pid = node();
  unsigned const number_of_processors = nodes();

  { // T=unsigned
    vector<unsigned> send(pid, pid);
    vector<vector<unsigned>> receive;

    indeterminate_gatherv(send, receive);
    PASSMSG("No exception thrown");

    if (pid == 0) {
      if (receive.size() == number_of_processors)
        PASSMSG("correct number of processors in gatherv");
      else
        FAILMSG("NOT correct number of processors in gatherv");
      for (unsigned p = 0; p < number_of_processors; ++p) {
        if (receive[p].size() != p) {
          FAILMSG("NOT correct number of elements in gatherv");
        } else {
          for (unsigned i = 0; i < p; ++i) {
            if (receive[p][i] != p)
              FAILMSG("NOT correct values in gatherv");
          }
        }
        // Prepare for next test
        receive[p].resize(0);
        receive[p].resize(2 * p, 3 * p);
      }
    }

    indeterminate_scatterv(receive, send);

    if (send.size() == 2 * pid)
      PASSMSG("correct number of processors in scatterv");
    else
      FAILMSG("NOT correct number of processors in scatterv");
    for (unsigned i = 0; i < 2 * pid; ++i) {
      if (send[i] != 3 * pid)
        FAILMSG("NOT correct values in scatterv");
    }
  }

  { // T=double
    vector<double> send(pid, pid);
    vector<vector<double>> receive;

    indeterminate_gatherv(send, receive);
    PASSMSG("No exception thrown for indeterminate_gatherv<double>");

    if (pid == 0) {
      if (receive.size() == number_of_processors)
        PASSMSG("correct number of processors in gatherv");
      else
        FAILMSG("NOT correct number of processors in gatherv");
      for (unsigned p = 0; p < number_of_processors; ++p) {
        if (receive[p].size() != p) {
          FAILMSG("NOT correct number of elements in gatherv");
        } else {
          for (unsigned i = 0; i < p; ++i) {
            if (!rtt_dsxx::soft_equiv(receive[p][i], static_cast<double>(p)))
              FAILMSG("NOT correct values in gatherv");
          }
        }
        // Prepare for next test
        receive[p].resize(0);
        receive[p].resize(2 * p, 3 * p);
      }
    }

    indeterminate_scatterv(receive, send);

    if (send.size() == 2 * pid)
      PASSMSG("correct number of processors in scatterv");
    else
      FAILMSG("NOT correct number of processors in scatterv");
    for (unsigned i = 0; i < 2 * pid; ++i) {
      if (!rtt_dsxx::soft_equiv(send[i], static_cast<double>(3 * pid)))
        FAILMSG("NOT correct values in scatterv");
    }
  }

  { // T=int
    vector<int> send(pid, pid);
    vector<vector<int>> receive;

    indeterminate_gatherv(send, receive);
    PASSMSG("No exception thrown for indeterminate_gatherv<int>");

    if (pid == 0) {
      if (receive.size() == number_of_processors)
        PASSMSG("correct number of processors in gatherv");
      else
        FAILMSG("NOT correct number of processors in gatherv");
      for (unsigned p = 0; p < number_of_processors; ++p) {
        if (receive[p].size() != p) {
          FAILMSG("NOT correct number of elements in gatherv");
        } else {
          for (unsigned i = 0; i < p; ++i) {
            if (receive[p][i] != static_cast<int>(p))
              FAILMSG("NOT correct values in gatherv");
          }
        }
        // Prepare for next test
        receive[p].resize(0);
        receive[p].resize(2 * p, 3 * p);
      }
    }

    indeterminate_scatterv(receive, send);

    if (send.size() == 2 * pid)
      PASSMSG("correct number of processors in scatterv");
    else
      FAILMSG("NOT correct number of processors in scatterv");
    for (unsigned i = 0; i < 2 * pid; ++i) {
      if (send[i] != static_cast<int>(3 * pid))
        FAILMSG("NOT correct values in scatterv");
    }
  }

  // Test with empty container
  {
    vector<unsigned> emptysend;
    vector<vector<unsigned>> emptyreceive;

    indeterminate_gatherv(emptysend, emptyreceive);
    PASSMSG("No exception thrown for indeterminate_gatherv with empty containers.");

    FAIL_IF(emptysend.size() != 0);
    FAIL_IF(emptyreceive.size() != number_of_processors);
    FAIL_IF(emptyreceive[pid].size() != 0);

    indeterminate_scatterv(emptyreceive, emptysend);

    FAIL_IF(emptysend.size() != 0);
    FAIL_IF(emptyreceive.size() != number_of_processors);
    FAIL_IF(emptyreceive[pid].size() != 0);
  }

  return;
}

//------------------------------------------------------------------------------------------------//
void tstDeterminateGatherScatterv(UnitTest &ut) {
  unsigned const pid = node();
  unsigned const number_of_processors = nodes();

  { // unsigned
    vector<unsigned> send(pid, pid);
    vector<vector<unsigned>> receive(number_of_processors);
    for (unsigned p = 0; p < number_of_processors; ++p) {
      receive[p].resize(p, p);
    }
    determinate_gatherv(send, receive);

    PASSMSG("No exception thrown <unsigned>");

    if (pid == 0) {
      if (receive.size() == number_of_processors)
        PASSMSG("correct number of processors in gatherv");
      else
        FAILMSG("NOT correct number of processors in gatherv");
      for (unsigned p = 0; p < number_of_processors; ++p) {
        if (receive[p].size() != p) {
          FAILMSG("NOT correct number of elements in gatherv");
        } else {
          for (unsigned i = 0; i < p; ++i) {
            if (receive[p][i] != p)
              FAILMSG("NOT correct values in gatherv");
          }
        }
        // Prepare for next test
        receive[p].resize(0);
        receive[p].resize(2 * p, 3 * p);
      }
    }

    send.resize(2 * pid);

    determinate_scatterv(receive, send);

    if (send.size() == 2 * pid)
      PASSMSG("correct number of processors in scatterv");
    else
      FAILMSG("NOT correct number of processors in scatterv");
    for (unsigned i = 0; i < 2 * pid; ++i) {
      if (send[i] != 3 * pid)
        FAILMSG("NOT correct values in scatterv");
    }
  }

  { // double
    vector<double> send(pid, pid);
    vector<vector<double>> receive(number_of_processors);
    for (unsigned p = 0; p < number_of_processors; ++p) {
      receive[p].resize(p, p);
    }
    determinate_gatherv(send, receive);

    PASSMSG("No exception thrown for determinate_gatherv<double>");

    if (pid == 0) {
      if (receive.size() == number_of_processors)
        PASSMSG("correct number of processors in gatherv");
      else
        FAILMSG("NOT correct number of processors in gatherv");
      for (unsigned p = 0; p < number_of_processors; ++p) {
        if (receive[p].size() != p) {
          FAILMSG("NOT correct number of elements in gatherv");
        } else {
          for (unsigned i = 0; i < p; ++i) {
            if (!rtt_dsxx::soft_equiv(receive[p][i], static_cast<double>(p)))
              FAILMSG("NOT correct values in gatherv");
          }
        }
        // Prepare for next test
        receive[p].resize(0);
        receive[p].resize(2 * p, 3 * p);
      }
    }

    send.resize(2 * pid);

    determinate_scatterv(receive, send);

    if (send.size() == 2 * pid)
      PASSMSG("correct number of processors in scatterv");
    else
      FAILMSG("NOT correct number of processors in scatterv");
    for (unsigned i = 0; i < 2 * pid; ++i) {
      if (!rtt_dsxx::soft_equiv(send[i], static_cast<double>(3 * pid)))
        FAILMSG("NOT correct values in scatterv");
    }
  }

  { // int
    vector<int> send(pid, pid);
    vector<vector<int>> receive(number_of_processors);
    for (unsigned p = 0; p < number_of_processors; ++p) {
      receive[p].resize(p, p);
    }
    determinate_gatherv(send, receive);

    PASSMSG("No exception thrown for determinate_gatherv<int>");

    if (pid == 0) {
      if (receive.size() == number_of_processors)
        PASSMSG("correct number of processors in gatherv");
      else
        FAILMSG("NOT correct number of processors in gatherv");
      for (unsigned p = 0; p < number_of_processors; ++p) {
        if (receive[p].size() != p) {
          FAILMSG("NOT correct number of elements in gatherv");
        } else {
          for (unsigned i = 0; i < p; ++i) {
            if (receive[p][i] != static_cast<int>(p))
              FAILMSG("NOT correct values in gatherv");
          }
        }
        // Prepare for next test
        receive[p].resize(0);
        receive[p].resize(2 * p, 3 * p);
      }
    }

    send.resize(2 * pid);

    determinate_scatterv(receive, send);

    if (send.size() == 2 * pid)
      PASSMSG("correct number of processors in scatterv");
    else
      FAILMSG("NOT correct number of processors in scatterv");
    for (unsigned i = 0; i < 2 * pid; ++i) {
      if (send[i] != static_cast<int>(3 * pid))
        FAILMSG("NOT correct values in scatterv");
    }
  }

  { // char
    vector<char> send(pid, 'A');
    vector<vector<char>> receive(number_of_processors);
    for (unsigned p = 0; p < number_of_processors; ++p) {
      receive[p].resize(p, 'B');
    }
    determinate_gatherv(send, receive);

    PASSMSG("No exception thrown for determinate_gatherv<char>");

    if (pid == 0) {
      if (receive.size() == number_of_processors)
        PASSMSG("correct number of processors in gatherv");
      else
        FAILMSG("NOT correct number of processors in gatherv");
      for (unsigned p = 0; p < number_of_processors; ++p) {
        if (receive[p].size() != p) {
          FAILMSG("NOT correct number of elements in gatherv");
        } else {
          for (unsigned i = 0; i < p; ++i) {
            if (receive[p][i] != 'A')
              FAILMSG("NOT correct values in gatherv");
          }
        }
        // Prepare for next test
        receive[p].resize(0);
        receive[p].resize(2 * p, 'C');
      }
    }
  }

  return;
}

//------------------------------------------------------------------------------------------------//
void topology_report(UnitTest &ut) {
  size_t const mpi_ranks = rtt_c4::nodes();
  size_t const my_mpi_rank = rtt_c4::node();

  if (my_mpi_rank == 0)
    std::cout << "\nStarting topology_report()..." << std::endl;

  // Store proc name on local proc
  std::string my_pname = rtt_c4::get_processor_name();
  size_t namelen = my_pname.size();

  // Create a container on IO proc to hold names of all nodes.
  vector<std::string> procnames(mpi_ranks);

  // Gather names into pnames on IO proc.
  rtt_c4::indeterminate_gatherv(my_pname, procnames);

  // Look at the data found on the IO proc.
  if (my_mpi_rank == 0) {

    FAIL_IF(procnames[my_mpi_rank].size() != namelen);

    // Count unique processors
    vector<string> unique_processor_names;
    for (size_t i = 0; i < mpi_ranks; ++i) {
      bool found(false);
      for (const auto &unique_procname : unique_processor_names)
        if (procnames[i] == unique_procname)
          found = true;
      if (!found)
        unique_processor_names.push_back(procnames[i]);
    }

    // Print a report
    std::cout << "\nWe are using " << mpi_ranks << " mpi rank(s) on "
              << unique_processor_names.size() << " unique nodes.";

    for (size_t i = 0; i < mpi_ranks; ++i) {
      std::cout << "\n  - MPI rank " << i << " is on " << procnames[i];
      FAIL_IF(procnames[i].size() < 1);
    }
    std::cout << std::endl;

    // Generate a map with the node name as the key and a list of MPI ranks as a
    // vector<int> of data.

    vector<vector<size_t>> map_proc_to_ranks(unique_processor_names.size());
    for (size_t i = 0; i < mpi_ranks; ++i) {
      size_t node_number(0);
      for (size_t j = 0; j < unique_processor_names.size(); ++j)
        if (procnames[i] == unique_processor_names[j]) {
          node_number = j;
          break;
        }
      map_proc_to_ranks[node_number].push_back(i);
    }

    std::cout << "\nMPI ranks per node:";
    for (size_t j = 0; j < unique_processor_names.size(); ++j) {
      std::cout << "\n  - Node " << j << " (" << unique_processor_names[j] << "): ";
      for (auto const i : map_proc_to_ranks[j])
        std::cout << map_proc_to_ranks[j][i] << ",";
    }
    std::cout << std::endl;
  }
  return;
}

//------------------------------------------------------------------------------------------------//
void tstDeterminateAllGatherv(UnitTest &ut) {
  unsigned const pid = node();
  unsigned const number_of_processors = nodes();

  { // T=unsigned
    vector<unsigned> send(pid, pid);

    // determinate_allgatherv already has the data sizes from the other MPI/C4 ranks/nodes
    vector<vector<unsigned>> receive(number_of_processors);
    for (unsigned p = 0; p < number_of_processors; ++p)
      receive[p].resize(p);

    // gather data from all nodes to all nodes at rank index in receive vector
    determinate_allgatherv(send, receive);
    PASSMSG("No exception thrown");

    // check values gathered from each rank (and check on each rank of course)
    for (unsigned p = 0; p < number_of_processors; ++p) {
      for (unsigned i = 0; i < p; ++i) {
        if (receive[p][i] != p)
          FAILMSG("NOT correct values in allgatherv");
      }
    }
  }

  { // T=double
    vector<double> send(pid, static_cast<double>(pid));

    // determinate_allgatherv already has the data sizes from the other MPI/C4 ranks/nodes
    vector<vector<double>> receive(number_of_processors);
    for (unsigned p = 0; p < number_of_processors; ++p)
      receive[p].resize(p);

    // gather data from all nodes to all nodes at rank index in receive vector
    determinate_allgatherv(send, receive);
    PASSMSG("No exception thrown");

    // check values gathered from each rank (and check on each rank of course)
    for (unsigned p = 0; p < number_of_processors; ++p) {
      const auto p_dbl = static_cast<double>(p);
      for (unsigned i = 0; i < p; ++i) {
        if (!rtt_dsxx::soft_equiv(receive[p][i], p_dbl))
          FAILMSG("NOT correct values in allgatherv");
      }
    }
  }

  // successful test output
  if (ut.numFails == 0)
    PASSMSG("tstDeterminateAllGatherv tests ok.");
  return;
}

//------------------------------------------------------------------------------------------------//
void tstIndeterminateAllGatherv(UnitTest &ut) {
  unsigned const pid = node();
  unsigned const number_of_processors = nodes();

  // T=unsigned
  vector<unsigned> send(pid, pid);
  vector<vector<unsigned>> receive;

  // gather data from all nodes to all nodes at rank index in receive vector
  indeterminate_allgatherv(send, receive);
  PASSMSG("No exception thrown");

  // check the size of the receiving vector is now the number of MPI/C4 ranks/nodes
  if (receive.size() == number_of_processors)
    PASSMSG("correct number of processors in allgatherv");
  else
    FAILMSG("NOT correct number of processors in allgatherv");

  // check values gathered from each rank (and check on each rank of course)
  for (unsigned p = 0; p < number_of_processors; ++p) {
    if (receive[p].size() != p) {
      FAILMSG("NOT correct number of elements in allgatherv");
    } else {
      for (unsigned i = 0; i < p; ++i) {
        if (receive[p][i] != p)
          FAILMSG("NOT correct values in allgatherv");
      }
    }
  }

  // successful test output
  if (ut.numFails == 0)
    PASSMSG("tstIndeterminateAllGatherv tests ok.");
  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, release);
  try {
    tstAllgather(ut);
    tstDeterminateGatherScatter(ut);
    tstIndeterminateGatherScatterv(ut);
    tstDeterminateGatherScatterv(ut);
    topology_report(ut);
    tstDeterminateAllGatherv(ut);
    tstIndeterminateAllGatherv(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstGatherScatter.cc
//------------------------------------------------------------------------------------------------//
