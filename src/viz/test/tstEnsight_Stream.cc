//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   viz/test/tstEnsight_Stream.cc
 * \author Rob Lowrie
 * \date   Fri Nov 12 22:52:46 2004
 * \brief  Test for Ensight_Stream.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Packing_Utils.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include "viz/Ensight_Stream.hh"
#include <array>

using namespace std;
using rtt_viz::Ensight_Stream;

//------------------------------------------------------------------------------------------------//
// Utility functions
//------------------------------------------------------------------------------------------------//

// Reads binary value from stream.
template <typename T> void binary_read(ifstream &stream, T &v) {
  char *vc = new char[sizeof(T)];
  stream.read(vc, sizeof(T));

  rtt_dsxx::Unpacker p;
  p.set_buffer(sizeof(T), vc);
  p.unpack(v);

  delete[] vc;
}

// Various overloaded read functions.

void readit(ifstream &stream, const bool binary, double &d) {
  if (binary) {
    float x;
    binary_read(stream, x);
    d = static_cast<double>(x);
  } else
    stream >> d;
}

void readit(ifstream &stream, const bool binary, int &d) {
  if (binary)
    binary_read(stream, d);
  else
    stream >> d;
}

void readit(ifstream &stream, const bool binary, string &s) {
  if (binary) {
    s.resize(80);
    for (int i = 0; i < 80; ++i)
      stream.read(&s[i], 1);
  } else
    stream >> s;
}

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void test_simple(rtt_dsxx::UnitTest &ut, bool const binary, bool const geom,
                 bool const decomposed) {
  // Dump a few values into the stream

  const int i(20323);
  const string s("dog");
  const double d(112.3);
  const string file("ensight_stream_" + std::to_string(rtt_c4::nodes()) + ".out");

  {
    Ensight_Stream f(file, binary, geom, decomposed);

    f << i << rtt_viz::endl;
    f << d << rtt_viz::endl;
    f << s << rtt_viz::endl;
    f.flush();
  }

  // Read the file back in using rank 0 and check the values.
  if (rtt_c4::node() == 0) {

    std::ios::openmode mode = std::ios::in;

    if (binary) {
      cout << "Testing binary mode." << endl;
      mode = mode | std::ios::binary;
    } else
      cout << "Testing ascii mode." << endl;

    ifstream in(file.c_str(), mode);
    // read out the "C Binary" data
    if (binary && geom) {
      std::array<char, 8> buf{};
      in.read(buf.data(), sizeof(char) * buf.size());
      FAIL_IF_NOT(std::string("C Binary") == std::string(buf.begin(), buf.end()));
    }

    int i_in;
    readit(in, binary, i_in);
    FAIL_IF_NOT(i == i_in);

    double d_in;
    readit(in, binary, d_in);
    // floats are inaccurate
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(d, d_in, 0.01));

    string s_in;
    readit(in, binary, s_in);
    for (size_t k = 0; k < s.size(); ++k)
      FAIL_IF_NOT(s[k] == s_in[k]);
  }

  if (ut.numFails == 0)
    PASSMSG("test_simple() completed successfully.");
  else
    FAILMSG("test_simple() did not complete successfully.");

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try { // >>> UNIT TESTS
    // serial/replicated use test
    if (rtt_c4::node() == 0) {
      test_simple(ut, true, false, false);  // test binary
      test_simple(ut, false, false, false); // test ascii
      test_simple(ut, true, false, false);  // test binary with geom flag
    }

    // Wait for rank 0 to finish serial testing before proceeding:
    rtt_c4::global_barrier();

    // parallel/decomposition tests in decomposition mode
    test_simple(ut, true, false, true);  // test binary
    test_simple(ut, false, false, true); // test ascii
    test_simple(ut, true, false, true);  // test binary with geom flag
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstEnsight_Stream.cc
//------------------------------------------------------------------------------------------------//
