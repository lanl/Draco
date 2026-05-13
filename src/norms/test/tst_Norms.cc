//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   norms/test/tst_Norms.cc
 * \author Rob Lowrie
 * \brief  Tests Norms.
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "c4/C4_Functions.hh"
#include "c4/ParallelUnitTest.hh"
#include "norms/Norms.hh"
#include "norms/Norms_Labeled.hh"
#include "norms/Norms_Proc.hh"
#include <cmath>
#include <cstddef>
#include <sstream>

// cSpell:ignore Linf

namespace {

using namespace std;

using rtt_dsxx::soft_equiv;
using rtt_norms::Norms;
using rtt_norms::Norms_Labeled;
using rtt_norms::Norms_Proc;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

// Test Norms
void test_Norms(rtt_dsxx::UnitTest &ut) {
  Norms norms;

  // All norms should be the same if we add just a single quantity, in this case, the processor
  // number + 1.

  size_t const x = rtt_c4::node();
  double const xp1 = static_cast<double>(x) + 1.0;

  norms.add(xp1, x);

  FAIL_IF_NOT(soft_equiv(norms.L1(), xp1));
  FAIL_IF_NOT(soft_equiv(norms.L2(), xp1));
  FAIL_IF_NOT(soft_equiv(norms.Linf(), xp1));
  FAIL_IF_NOT(norms.index_Linf() == x);

  // Check accumulating the results back to the host proc.

  norms.comm(0);

  if (x == 0) {
    size_t const num_nodes = rtt_c4::nodes();

    FAIL_IF_NOT(norms.index_Linf() == num_nodes - 1);
    FAIL_IF_NOT(soft_equiv(norms.Linf(), double(num_nodes)));

    double tL1 = 0.0;
    double tL2 = 0.0;

    for (size_t i = 0; i < num_nodes; ++i) {
      double const xt = static_cast<double>(i) + 1.0;
      tL1 += xt;
      tL2 += xt * xt;
    }

    tL1 /= static_cast<double>(num_nodes);
    tL2 = std::sqrt(tL2 / static_cast<double>(num_nodes));

    FAIL_IF_NOT(soft_equiv(norms.L1(), tL1));
    FAIL_IF_NOT(soft_equiv(norms.L2(), tL2));
  }

  { // test assignment

    Norms nc;
    nc = norms;
    FAIL_IF_NOT(nc == norms);
  }

  { // test copy ctor.

    Norms const nc(norms); // NOLINT(performance-unnecessary-copy-initialization)
    FAIL_IF_NOT(nc == norms);
  }

  // Done testing

  if (ut.numFails == 0)
    PASSMSG("test_Norms() ok.");
  return;
}

//------------------------------------------------------------------------------------------------//
// Test Norms_Labeled
void test_Norms_Labeled(rtt_dsxx::UnitTest &ut) {
  Norms_Labeled norms;

  // All norms should be the same if we add just a single quantity, in this case, the processor
  // number + 1.

  size_t const x = rtt_c4::node();
  double const xp1 = static_cast<double>(x) + 1.0;

  Norms_Labeled::Index indx(x);
  ostringstream o;
  o << "proc " << x;
  indx.label = o.str();

  norms.add(xp1, indx);

  FAIL_IF_NOT(soft_equiv(norms.L1(), xp1));
  FAIL_IF_NOT(soft_equiv(norms.L2(), xp1));
  FAIL_IF_NOT(soft_equiv(norms.Linf(), xp1));
  FAIL_IF_NOT(norms.index_Linf() == indx);

  // Check accumulating the results back to the host proc.

  norms.comm(0);

  if (x == 0) {
    size_t const num_nodes = rtt_c4::nodes();

    //FAIL_IF_NOT(norms.index_Linf() == num_nodes-1);
    FAIL_IF_NOT(soft_equiv(norms.Linf(), double(num_nodes)));

    double tL1 = 0.0;
    double tL2 = 0.0;

    for (size_t i = 0; i < num_nodes; ++i) {
      double const xt = static_cast<double>(i) + 1.0;
      tL1 += xt;
      tL2 += xt * xt;
    }

    tL1 /= static_cast<double>(num_nodes);
    tL2 = std::sqrt(tL2 / static_cast<double>(num_nodes));

    FAIL_IF_NOT(soft_equiv(norms.L1(), tL1));
    FAIL_IF_NOT(soft_equiv(norms.L2(), tL2));
  }

  { // test assignment

    Norms_Labeled nc;
    nc = norms;
    FAIL_IF_NOT(nc == norms);
  }

  { // test copy ctor.

    Norms_Labeled const nc(norms); // NOLINT(performance-unnecessary-copy-initialization)
    FAIL_IF_NOT(nc == norms);
  }

  // Done testing

  if (ut.numFails == 0)
    PASSMSG("test_Norms_Labeled() ok.");
  return;
}

//------------------------------------------------------------------------------------------------//
// Test Norms_Proc
void test_Norms_Proc(rtt_dsxx::UnitTest &ut) {
  Norms_Proc norms;

  // All norms should be the same if we add just a single quantity, in this case, the processor
  // number + 1.

  size_t const x = rtt_c4::node();
  double const xp1 = static_cast<double>(x) + 1.0;

  Norms_Proc::Index const indx(x);

  norms.add(xp1, indx);

  FAIL_IF_NOT(soft_equiv(norms.L1(), xp1));
  FAIL_IF_NOT(soft_equiv(norms.L2(), xp1));
  FAIL_IF_NOT(soft_equiv(norms.Linf(), xp1));
  FAIL_IF_NOT(norms.index_Linf() == indx);

  // Check accumulating the results back to the host proc.

  norms.comm(0);

  if (x == 0) {
    size_t const num_nodes = rtt_c4::nodes();

    //FAIL_IF_NOT(norms.index_Linf() == num_nodes-1);
    FAIL_IF_NOT(soft_equiv(norms.Linf(), double(num_nodes)));

    double tL1 = 0.0;
    double tL2 = 0.0;

    for (size_t i = 0; i < num_nodes; ++i) {
      double const xt = static_cast<double>(i) + 1.0;
      tL1 += xt;
      tL2 += xt * xt;
    }

    tL1 /= static_cast<double>(num_nodes);
    tL2 = std::sqrt(tL2 / static_cast<double>(num_nodes));

    FAIL_IF_NOT(soft_equiv(norms.L1(), tL1));
    FAIL_IF_NOT(soft_equiv(norms.L2(), tL2));
  }

  { // test assignment

    Norms_Proc nc;
    nc = norms;
    FAIL_IF_NOT(nc == norms);
  }

  { // test copy ctor.

    Norms_Proc const nc(norms); // NOLINT(performance-unnecessary-copy-initialization)
    FAIL_IF_NOT(nc == norms);
  }

  // Done testing

  if (ut.numFails == 0)
    PASSMSG("test_Norms_Proc() ok.");
  return;
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_Norms(ut);
    test_Norms_Labeled(ut);
    test_Norms_Proc(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tst_Norms.cc
//------------------------------------------------------------------------------------------------//
