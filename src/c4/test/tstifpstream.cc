//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstifpstream.cc
 * \author Mathew Cleveland
 * \date   Feb. 2020
 * \brief  Tests ifpstream class for paralle read
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "c4/ifpstream.hh"
#include "c4/ofpstream.hh"
#include "ds++/Release.hh"
#include <sstream>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;

//------------------------------------------------------------------------------------------------//
// HELPER FUNCTIONS
//------------------------------------------------------------------------------------------------//

void write_stream(std::string const &filename, const std::ios_base::openmode mode) {

  unsigned const pid = rtt_c4::node();

  ofpstream out(filename, mode);

  if (pid == 0)
    out << "just rank 0" << endl;

  out << "MPI rank " << pid << " reporting ..." << endl;
  out.send();
  out.shrink_to_fit();

  out << "MPI rank " << pid << " reporting a second time ..." << endl;
  out.shrink_to_fit();
  out.send();
}

//------------------------------------------------------------------------------------------------//
void tstifpstream(UnitTest &ut, const std::ios_base::openmode mode) {

  int pid = rtt_c4::node();

  std::string filename = "tstifpstream" + std::to_string(rtt_c4::nodes()) + ".txt";
  if (mode == std::ofstream::binary)
    filename = "tstifpstream" + std::to_string(rtt_c4::nodes()) + ".bin";

  write_stream(filename, mode);
  // Read file on head rank, check for correct conversion and ordering
  ifpstream in(filename, mode);
  unsigned zeroth_buffer_size = 12;
  if (pid != 0)
    zeroth_buffer_size = 0;
  unsigned first_buffer_size = 25;
  unsigned second_buffer_size = 39;

  // read first strings only on proc 0
  std::string compare_string;
  if (pid == 0)
    compare_string = "just rank 0\n";
  std::string zeroth_string(zeroth_buffer_size, '\0');
  in.fill_buffers(zeroth_buffer_size);
  in.read(&zeroth_string[0], zeroth_buffer_size);
  FAIL_IF_NOT(zeroth_string.compare(compare_string) == 0);

  // read first strings on each rank
  compare_string = "MPI rank " + to_string(pid) + " reporting ...\n";
  std::string first_string(first_buffer_size, '\0');
  in.fill_buffers(first_buffer_size);
  in.read(&first_string[0], first_buffer_size);
  FAIL_IF_NOT(first_string.compare(compare_string) == 0);

  // read second strings on each rank
  compare_string = "MPI rank " + to_string(pid) + " reporting a second time ...\n";
  std::string second_string(second_buffer_size, '\0');
  in.fill_buffers(second_buffer_size);
  in.read(&second_string[0], second_buffer_size);
  FAIL_IF_NOT(second_string.compare(compare_string) == 0);

  if (pid == 0)
    PASSMSG("completed serialized read of " + filename + " without hanging or segfaulting");
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, release);
  try {
    tstifpstream(ut, std::ofstream::out);
    tstifpstream(ut, std::ofstream::binary);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstofpstream.cc
//------------------------------------------------------------------------------------------------//
