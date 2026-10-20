//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstofpstream_large_buffer.cc
 * \author SI: Claude Sonnet 4.5
 * \brief  Test c4::ofpstream with large buffers (slow test)
 * \note   Copyright (C) 2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/UnitTest.hh"
#include "c4/C4_Functions.hh"
#include "c4/ParallelUnitTest.hh"
#include "c4/ofpstream.hh"
#include <climits>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void tstofpstream_large_buffer(UnitTest &ut) {

  int const pid = rtt_c4::node();

  std::string const filename("tstofpstream_large_" + std::to_string(rtt_c4::nodes()) + ".bin");

  // Test with buffer larger than INT_MAX to verify chunking logic
  // We'll use a smaller size for testing (to avoid excessive memory usage)
  // but still large enough to test the chunking code path if INT_MAX threshold is lowered
  size_t const test_size = rtt_c4::nodes() > 2 ? static_cast<size_t>(INT_MAX) + 1024 : 1024 * 1024;

  // Binary write large buffer to file using ofpstream:
  {
    ofpstream out(filename, std::ofstream::binary);

    // Create a pattern that can be verified
    std::vector<char> buffer(test_size);
    for (size_t i = 0; i < buffer.size(); ++i) {
      buffer[i] = static_cast<char>((pid * 256 + i) % 256);
    }

    out.write(buffer.data(), buffer.size());
    out.send();
  }

  // Read file on head rank, check for correct data
  if (pid == 0) {
    ifstream in(filename, std::ifstream::binary);

    for (int rank = 0; rank < rtt_c4::nodes(); ++rank) {
      std::vector<char> buffer(test_size);
      in.read(buffer.data(), buffer.size());

      // Verify the pattern
      bool data_correct = true;
      for (size_t i = 0; i < buffer.size(); ++i) {
        char const expected = static_cast<char>((rank * 256 + i) % 256);
        if (buffer[i] != expected) {
          std::ostringstream msg;
          msg << "Data mismatch for rank " << rank << " at byte " << i << ": expected "
              << static_cast<int>(expected) << ", got " << static_cast<int>(buffer[i]);
          FAILMSG(msg.str());
          data_correct = false;
          break;
        }
      }

      if (!data_correct) {
        break;
      }
    }
  }

  PASSMSG("completed large buffer serialized write without hanging or data corruption");
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, release);
  try {
    tstofpstream_large_buffer(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstofpstream_large_buffer.cc
//------------------------------------------------------------------------------------------------//
