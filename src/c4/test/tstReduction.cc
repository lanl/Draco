//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/test/tstReduction.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 15:41:00 2002
 * \brief  C4 Reduction test.
 * \note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "dsxx/UnitTest.hh"
#include "c4/C4_Functions.hh"
#include "c4/ParallelUnitTest.hh"
#include <algorithm> // std::equal
#include <cstddef>   // size_t
#include <cstdint>   // uint32_t
#include <iostream>
#include <limits>
#include <vector>

namespace {

using namespace std;

using rtt_c4::C4_Req;
using rtt_c4::global_and;
using rtt_c4::global_isum;
using rtt_c4::global_max;
using rtt_c4::global_min;
using rtt_c4::global_prod;
using rtt_c4::global_sum;
using rtt_c4::prefix_sum;
using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void elemental_reduction(rtt_dsxx::UnitTest &ut) {

  // test bools with with blocking mpi_land
  bool are_all_proc_even = (rtt_c4::node() + 1) % 2 == 0;
  global_and(are_all_proc_even);

  // at least one processor index must be odd (adding 1)
  FAIL_IF(are_all_proc_even);

  // test ints with blocking and non-blocking sums
  int xint = rtt_c4::node() + 1;
  global_sum(xint);

  int xint_send = rtt_c4::node() + 1;
  int xint_recv = 0;
  C4_Req int_request;
  global_isum(xint_send, xint_recv, int_request);
  int_request.wait();

  int int_answer = 0;
  for (int i = 0; i < rtt_c4::nodes(); i++)
    int_answer += i + 1;

  FAIL_IF_NOT(xint == int_answer);
  FAIL_IF_NOT(xint_recv == int_answer);
  if (!rtt_c4::node())
    cout << "int: Global non-blocking sum: " << xint_recv << " answer: " << int_answer << '\n';

  // Test with deprecated form of global_sum
  xint = rtt_c4::node() + 1;
  global_sum(xint);
  FAIL_IF(xint != int_answer);

  // test longs for blocking and non-blocking sums
  long const max_long(numeric_limits<long>::max());
  int64_t const ten_billion(10000000000L); // 1e10 > MAX_INT
  int32_t const one_billion(1000000000L);  // 1e9 < MAX_INT

  long x_long = rtt_c4::node() + (max_long > ten_billion ? ten_billion : one_billion);
  global_sum(x_long);

  long x_long_send = rtt_c4::node() + (max_long > ten_billion ? ten_billion : one_billion);
  long x_long_recv = 0;
  C4_Req long_request;
  global_isum(x_long_send, x_long_recv, long_request);
  long_request.wait();

  long long_answer = 0;
  for (int i = 0; i < rtt_c4::nodes(); i++)
    long_answer += i + (max_long > ten_billion ? ten_billion : one_billion);

  FAIL_IF(x_long != long_answer);

  FAIL_IF(x_long_recv != long_answer);

  // test doubles for blocking and non-blocking sums
  double xdbl = static_cast<double>(rtt_c4::node()) + 0.1;
  global_sum(xdbl);

  double xdouble_send = static_cast<double>(rtt_c4::node()) + 0.1;
  double xdouble_recv = 0;
  C4_Req double_request;
  global_isum(xdouble_send, xdouble_recv, double_request);
  double_request.wait();

  double dbl_answer = 0.0;
  for (int i = 0; i < rtt_c4::nodes(); i++)
    dbl_answer += static_cast<double>(i) + 0.1;

  FAIL_IF(!soft_equiv(xdbl, dbl_answer));
  FAIL_IF(!soft_equiv(xdouble_recv, dbl_answer));

  // test product
  x_long = rtt_c4::node() + 1;
  global_prod(x_long);

  long_answer = 1;
  for (int i = 0; i < rtt_c4::nodes(); i++)
    long_answer *= (i + 1);
  FAIL_IF_NOT(x_long == long_answer);

  // Test with deprecated form of global_prod
  x_long = rtt_c4::node() + 1;
  global_prod(x_long);
  FAIL_IF_NOT(x_long == long_answer);

  // test min
  xdbl = 0.5 + rtt_c4::node();
  global_min(xdbl);
  FAIL_IF_NOT(soft_equiv(xdbl, 0.5));

  // Test with deprecated form of global_min
  xdbl = rtt_c4::node() + 0.5;
  global_min(xdbl);
  FAIL_IF_NOT(soft_equiv(xdbl, 0.5));

  // test max
  xdbl = 0.7 + rtt_c4::node();
  global_max(xdbl);
  FAIL_IF_NOT(soft_equiv(xdbl, rtt_c4::nodes() - 0.3));

  // Test with deprecated form of global_max
  xdbl = 0.7 + rtt_c4::node();
  global_max(xdbl);
  FAIL_IF_NOT(soft_equiv(xdbl, rtt_c4::nodes() - 0.3));

  { // T = float

    float xflt = static_cast<float>(rtt_c4::node()) + 0.1F;
    global_sum(xflt);

    float x_float_send = static_cast<float>(rtt_c4::node()) + 0.1F;
    float x_float_recv = 0;
    C4_Req float_request;
    global_isum(x_float_send, x_float_recv, float_request);
    float_request.wait();

    float flt_answer = 0.0;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      flt_answer += static_cast<float>(i) + 0.1F;

    FAIL_IF_NOT(soft_equiv(xflt, flt_answer));
    FAIL_IF_NOT(soft_equiv(x_float_recv, flt_answer));

    // test product
    xflt = static_cast<float>(rtt_c4::node()) + 0.1F;
    global_prod(xflt);

    flt_answer = 1.0F;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      flt_answer *= static_cast<float>(i) + 0.1F;
    FAIL_IF_NOT(soft_equiv(xflt, flt_answer));

    // test min
    xflt = static_cast<float>(rtt_c4::node()) + 0.5F;
    global_min(xflt);
    FAIL_IF_NOT(soft_equiv(xflt, 0.5F));

    // test max
    xflt = 0.7F + static_cast<float>(rtt_c4::node());
    global_max(xflt);

    FAIL_IF_NOT(soft_equiv(xflt, static_cast<float>(rtt_c4::nodes()) - 0.3F));
  }
  { // T = double

    double xflt = static_cast<double>(rtt_c4::node()) + 0.1;
    global_sum(xflt);

    double x_float_send = static_cast<double>(rtt_c4::node()) + 0.1;
    double x_float_recv = 0;
    C4_Req float_request;
    global_isum(x_float_send, x_float_recv, float_request);
    float_request.wait();

    double flt_answer = 0.0;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      flt_answer += static_cast<double>(i) + 0.1;

    FAIL_IF_NOT(soft_equiv(xflt, flt_answer));
    FAIL_IF_NOT(soft_equiv(x_float_recv, flt_answer));

    // test product
    xflt = static_cast<double>(rtt_c4::node()) + 0.1;
    global_prod(xflt);

    flt_answer = 1.0;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      flt_answer *= static_cast<double>(i) + 0.1;
    FAIL_IF_NOT(soft_equiv(xflt, flt_answer));

    // Test with deprecated form of global_prod
    xflt = static_cast<double>(rtt_c4::node()) + 0.1;
    global_prod(xflt);
    FAIL_IF_NOT(soft_equiv(xflt, flt_answer));

    // test min
    xflt = static_cast<double>(rtt_c4::node()) + 0.5;
    global_min(xflt);

    FAIL_IF_NOT(soft_equiv(xflt, 0.5));

    // Test with deprecated form of global_min
    xflt = static_cast<double>(rtt_c4::node()) + 0.5;
    global_min(xflt);
    FAIL_IF_NOT(soft_equiv(xflt, 0.5));

    // test max
    xflt = 0.7 + static_cast<double>(rtt_c4::node());
    global_max(xflt);

    FAIL_IF_NOT(soft_equiv(xflt, static_cast<double>(rtt_c4::nodes()) - 0.3));

    // Test with deprecated form of global_max
    xflt = 0.7 + static_cast<double>(rtt_c4::node());
    global_max(xflt);
    FAIL_IF_NOT(soft_equiv(xflt, static_cast<double>(rtt_c4::nodes()) - 0.3));
  }
#ifndef CRAYPE_CCE
  { // T = long double

    long double xld = static_cast<long double>(rtt_c4::node()) + 0.1L;
    global_sum(xld);

    long double xld_send = static_cast<long double>(rtt_c4::node()) + 0.1L;
    long double xld_recv = 0;
    C4_Req ld_request;
    global_isum(xld_send, xld_recv, ld_request);
    ld_request.wait();

    auto ld_answer = static_cast<long double>(0.0);
    for (int i = 0; i < rtt_c4::nodes(); i++)
      ld_answer += static_cast<long double>(i) + 0.1L;

    FAIL_IF_NOT(soft_equiv(xld, ld_answer));
    FAIL_IF_NOT(soft_equiv(xld_recv, ld_answer));

    // test product
    xld = static_cast<long double>(rtt_c4::node()) + 0.1L;
    global_prod(xld);

    ld_answer = 1.0L;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      ld_answer *= static_cast<long double>(i) + 0.1L;
    FAIL_IF_NOT(soft_equiv(xld, ld_answer));

    // test min
    xld = static_cast<long double>(rtt_c4::node()) + 0.5L;
    global_min(xld);
    FAIL_IF_NOT(soft_equiv(xld, 0.5L));

    // test max
    xld = 0.7L + static_cast<long double>(rtt_c4::node());
    global_max(xld);

    FAIL_IF_NOT(soft_equiv(xld, static_cast<long double>(rtt_c4::nodes()) - 0.3L));
  }
#endif
  { // T = int

    xint = rtt_c4::node() + 1;
    global_sum(xint);

    xint_send = rtt_c4::node() + 1;
    xint_recv = 0;
    // C4_Req int_request;
    global_isum(xint_send, xint_recv, int_request);
    int_request.wait();

    int_answer = 0;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      int_answer += i + 1;
    FAIL_IF_NOT(xint == int_answer);
    FAIL_IF_NOT(xint_recv == int_answer);

    // test product
    xint = rtt_c4::node() + 1;
    global_prod(xint);

    int_answer = 1;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      int_answer *= (i + 1);
    FAIL_IF_NOT(xint == int_answer);

    // Test with deprecated form of global_prod
    xint = rtt_c4::node() + 1;
    global_prod(xint);
    FAIL_IF_NOT(xint == int_answer);

    // test min
    xint = rtt_c4::node() + 1;
    global_min(xint);
    FAIL_IF_NOT(xint == 1);

    // test max
    xint = rtt_c4::node() + 1;
    global_max(xint);
    FAIL_IF_NOT(xint == rtt_c4::nodes());
  }

  { // T = long

    x_long = rtt_c4::node() + 1;
    global_sum(x_long);

    x_long_send = rtt_c4::node() + 1;
    x_long_recv = 0;
    // C4_Req long_request;
    global_isum(x_long_send, x_long_recv, long_request);
    long_request.wait();

    long_answer = 0;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      long_answer += i + 1;
    FAIL_IF_NOT(x_long == long_answer);
    FAIL_IF_NOT(x_long_recv == long_answer);

    // test product
    x_long = rtt_c4::node() + 1;
    global_prod(x_long);

    long_answer = 1;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      long_answer *= (i + 1);
    FAIL_IF_NOT(x_long == long_answer);

    // Test with deprecated form of global_prod
    x_long = rtt_c4::node() + 1;
    global_prod(x_long);
    FAIL_IF_NOT(x_long == long_answer);

    // test min
    x_long = rtt_c4::node() + 1;
    global_min(x_long);
    FAIL_IF_NOT(x_long == 1);

    // test max
    x_long = rtt_c4::node() + 1;
    global_max(x_long);
    FAIL_IF_NOT(x_long == rtt_c4::nodes());
  }

  { // T = short
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#endif

    auto x_short = static_cast<short>(rtt_c4::node() + 1);
    global_sum(x_short);

    auto x_short_send = static_cast<short>(rtt_c4::node() + 1);
    short x_short_recv = 0;
    C4_Req short_request;
    global_isum(x_short_send, x_short_recv, short_request);
    short_request.wait();

    short short_answer = 0;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      short_answer += static_cast<short int>(i + 1);
    FAIL_IF_NOT(x_short == short_answer);
    FAIL_IF_NOT(x_short_recv == short_answer);

    // test product
    x_short = static_cast<short>(rtt_c4::node() + 1);
    global_prod(x_short);

    short_answer = 1;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      short_answer *= static_cast<short>(i + 1);
    FAIL_IF_NOT(x_short == short_answer);

    // Test with deprecated form of global_prod
    x_short = static_cast<short>(rtt_c4::node() + 1);
    global_prod(x_short);
    FAIL_IF_NOT(x_short == short_answer);

    // test min
    x_short = static_cast<short>(rtt_c4::node() + 1);
    global_min(x_short);
    FAIL_IF_NOT(x_short == 1);

    // test max
    x_short = static_cast<short>(rtt_c4::node() + 1);
    global_max(x_short);
    FAIL_IF_NOT(x_short == rtt_c4::nodes());

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
  }

  { // T = unsigned

    auto x_unsigned = static_cast<unsigned>(rtt_c4::node() + 1);
    global_sum(x_unsigned);

    auto x_unsigned_send = static_cast<unsigned>(rtt_c4::node() + 1);
    unsigned x_unsigned_recv = 0;
    C4_Req unsigned_request;
    global_isum(x_unsigned_send, x_unsigned_recv, unsigned_request);
    unsigned_request.wait();

    unsigned unsigned_answer = 0;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      unsigned_answer += static_cast<unsigned int>(i + 1);
    FAIL_IF_NOT(x_unsigned == unsigned_answer);
    FAIL_IF_NOT(x_unsigned_recv == unsigned_answer);

    // test product
    x_unsigned = static_cast<unsigned>(rtt_c4::node() + 1);
    global_prod(x_unsigned);

    unsigned_answer = 1;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      unsigned_answer *= static_cast<unsigned>(i + 1);
    FAIL_IF_NOT(x_unsigned == unsigned_answer);

    // Test with deprecated form of global_prod
    x_unsigned = static_cast<unsigned>(rtt_c4::node() + 1);
    global_prod(x_unsigned);
    FAIL_IF_NOT(x_unsigned == unsigned_answer);

    // test min
    x_unsigned = static_cast<unsigned>(rtt_c4::node() + 1);
    global_min(x_unsigned);
    FAIL_IF_NOT(x_unsigned == 1);

    // test max
    x_unsigned = static_cast<unsigned>(rtt_c4::node() + 1);
    global_max(x_unsigned);
    FAIL_IF_NOT(x_unsigned == rtt_c4::nranks());
  }

  { // T = unsigned long

    auto x_unsigned = static_cast<unsigned long>(rtt_c4::node() + 1);
    global_sum(x_unsigned);

    auto x_unsigned_send = static_cast<unsigned long>(rtt_c4::node() + 1);
    unsigned long x_unsigned_recv = 0;
    C4_Req unsigned_request;
    global_isum(x_unsigned_send, x_unsigned_recv, unsigned_request);
    unsigned_request.wait();

    unsigned long unsigned_answer = 0;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      unsigned_answer += static_cast<unsigned long>(i + 1);
    FAIL_IF_NOT(x_unsigned == unsigned_answer);
    FAIL_IF_NOT(x_unsigned_recv == unsigned_answer);

    // test product
    x_unsigned = static_cast<unsigned long>(rtt_c4::node() + 1);
    global_prod(x_unsigned);

    unsigned_answer = 1;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      unsigned_answer *= static_cast<unsigned long>(i + 1);
    FAIL_IF_NOT(x_unsigned == unsigned_answer);

    // Test with deprecated form of global_prod
    x_unsigned = static_cast<unsigned long>(rtt_c4::node() + 1);
    global_prod(x_unsigned);
    FAIL_IF_NOT(x_unsigned == unsigned_answer);

    // test min
    x_unsigned = static_cast<unsigned long>(rtt_c4::node() + 1);
    global_min(x_unsigned);
    FAIL_IF_NOT(x_unsigned == 1);

    // test max
    x_unsigned = static_cast<unsigned long>(rtt_c4::node() + 1);
    global_max(x_unsigned);
    FAIL_IF_NOT(x_unsigned == rtt_c4::nranks());
  }

  { // T = long long

    auto x_unsigned = static_cast<long long>(rtt_c4::node() + 1);
    global_sum(x_unsigned);

    auto x_unsigned_send = static_cast<long long>(rtt_c4::node() + 1);
    long long x_unsigned_recv = 0;
    C4_Req unsigned_request;
    global_isum(x_unsigned_send, x_unsigned_recv, unsigned_request);
    unsigned_request.wait();

    long long unsigned_answer = 0;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      unsigned_answer += static_cast<long long>(i + 1);
    FAIL_IF_NOT(x_unsigned == unsigned_answer);
    FAIL_IF_NOT(x_unsigned_recv == unsigned_answer);

    // test product
    x_unsigned = static_cast<long long>(rtt_c4::node() + 1);
    global_prod(x_unsigned);

    unsigned_answer = 1;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      unsigned_answer *= static_cast<long long>(i + 1);
    FAIL_IF_NOT(x_unsigned == unsigned_answer);

    // Test with deprecated form of global_prod
    x_unsigned = static_cast<long long>(rtt_c4::node() + 1);
    global_prod(x_unsigned);
    FAIL_IF_NOT(x_unsigned == unsigned_answer);

    // test min
    x_unsigned = static_cast<long long>(rtt_c4::node() + 1);
    global_min(x_unsigned);
    FAIL_IF_NOT(x_unsigned == 1);

    // test max
    x_unsigned = static_cast<long long>(rtt_c4::node() + 1);
    global_max(x_unsigned);
    FAIL_IF_NOT(x_unsigned == rtt_c4::nodes());
  }

  { // T = unsigned long long

    auto x_unsigned = static_cast<unsigned long long>(rtt_c4::node() + 1);
    global_sum(x_unsigned);

    auto x_unsigned_send = static_cast<unsigned long long>(rtt_c4::node() + 1);
    unsigned long long x_unsigned_recv = 0;
    C4_Req unsigned_request;
    global_isum(x_unsigned_send, x_unsigned_recv, unsigned_request);
    unsigned_request.wait();

    unsigned long long unsigned_answer = 0;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      unsigned_answer += static_cast<unsigned long long>(i + 1);
    FAIL_IF_NOT(x_unsigned == unsigned_answer);
    FAIL_IF_NOT(x_unsigned_recv == unsigned_answer);

    // test product
    x_unsigned = static_cast<unsigned long long>(rtt_c4::node() + 1);
    global_prod(x_unsigned);

    unsigned_answer = 1;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      unsigned_answer *= static_cast<unsigned long long>(i + 1);
    FAIL_IF_NOT(x_unsigned == unsigned_answer);

    // Test with deprecated form of global_prod
    x_unsigned = static_cast<unsigned long long>(rtt_c4::node() + 1);
    global_prod(x_unsigned);
    FAIL_IF_NOT(x_unsigned == unsigned_answer);

    // test min
    x_unsigned = static_cast<unsigned long long>(rtt_c4::node() + 1);
    global_min(x_unsigned);
    FAIL_IF_NOT(x_unsigned == 1);

    // test max
    x_unsigned = static_cast<unsigned long long>(rtt_c4::node() + 1);
    global_max(x_unsigned);
    FAIL_IF_NOT(x_unsigned == rtt_c4::nranks());
  }

  { // T = unsigned short
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#endif

    auto x_unsigned = static_cast<unsigned short>(rtt_c4::node() + 1);
    global_sum(x_unsigned);

    auto x_unsigned_send = static_cast<unsigned short>(rtt_c4::node() + 1);
    unsigned short x_unsigned_recv = 0;
    C4_Req unsigned_request;
    global_isum(x_unsigned_send, x_unsigned_recv, unsigned_request);
    unsigned_request.wait();

    unsigned short unsigned_answer = 0;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      unsigned_answer += static_cast<unsigned short>(i + 1);
    FAIL_IF_NOT(x_unsigned == unsigned_answer);
    FAIL_IF_NOT(x_unsigned_recv == unsigned_answer);

    // test product
    x_unsigned = static_cast<unsigned short>(rtt_c4::node() + 1);
    global_prod(x_unsigned);

    unsigned_answer = 1;
    for (int i = 0; i < rtt_c4::nodes(); i++)
      unsigned_answer *= static_cast<unsigned short>(i + 1);
    FAIL_IF_NOT(x_unsigned == unsigned_answer);

    // Test with deprecated form of global_prod
    x_unsigned = static_cast<unsigned short>(rtt_c4::node() + 1);
    global_prod(x_unsigned);
    FAIL_IF_NOT(x_unsigned == unsigned_answer);

    // test min
    x_unsigned = static_cast<unsigned short>(rtt_c4::node() + 1);
    global_min(x_unsigned);
    FAIL_IF_NOT(x_unsigned == 1);

    // test max
    x_unsigned = static_cast<unsigned short>(rtt_c4::node() + 1);
    global_max(x_unsigned);
    FAIL_IF_NOT(x_unsigned == rtt_c4::nranks());

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
  }

  if (ut.numFails == 0)
    PASSMSG("Elemental reductions ok.");
  return;
}

//------------------------------------------------------------------------------------------------//
void array_reduction(rtt_dsxx::UnitTest &ut) {
  {
    // make a vector of doubles
    vector<double> x(100);
    vector<double> prod(100, 1.0);
    vector<double> sum(100, 0.0);
    vector<double> lmin(100, 0.0);
    vector<double> lmax(100, 0.0);

    // fill it
    for (int i = 0; i < 100; i++) {
      x[i] = rtt_c4::node() + 0.11;
      for (int j = 0; j < rtt_c4::nodes(); j++) {
        sum[i] += (j + 0.11);
        prod[i] *= (j + 0.11);
      }
      lmin[i] = 0.11;
      lmax[i] = rtt_c4::nodes() + 0.11 - 1.0;
    }

    vector<double> c;

    {
      c = x;
      global_sum(&c[0], 100);
      FAIL_IF_NOT(soft_equiv(c.begin(), c.end(), sum.begin(), sum.end()));

      c = x;
      global_sum(&c[0], static_cast<size_t>(100));
      FAIL_IF_NOT(soft_equiv(c.begin(), c.end(), sum.begin(), sum.end()));

      c = x;
      global_prod(&c[0], 100);
      FAIL_IF_NOT(soft_equiv(c.begin(), c.end(), prod.begin(), prod.end()));

      c = x;
      global_min(&c[0], 100);
      FAIL_IF_NOT(soft_equiv(c.begin(), c.end(), lmin.begin(), lmin.end()));

      c = x;
      global_max(&c[0], 100);
      FAIL_IF_NOT(soft_equiv(c.begin(), c.end(), lmax.begin(), lmax.end()));
    }
  }
  { // T = float

    // make a vector of floats
    vector<float> x(100);
    vector<float> prod(100, 1.0);
    vector<float> sum(100, 0.0);
    vector<float> lmin(100, 0.0);
    vector<float> lmax(100, 0.0);

    // fill it
    for (int i = 0; i < 100; i++) {
      x[i] = static_cast<float>(rtt_c4::node()) + 0.11F;
      for (int j = 0; j < rtt_c4::nodes(); j++) {
        sum[i] += (static_cast<float>(j) + 0.11F);
        prod[i] *= (static_cast<float>(j) + 0.11F);
      }
      lmin[i] = 0.11F;
      lmax[i] = static_cast<float>(rtt_c4::nodes()) + 0.11F - 1.0F;
    }

    vector<float> c;
    float const eps = 1.0e-6F;

    {
      c = x;
      global_sum(&c[0], 100);
      FAIL_IF_NOT(soft_equiv(c.begin(), c.end(), sum.begin(), sum.end(), eps));

      c = x;
      global_sum(&c[0], static_cast<size_t>(100));
      FAIL_IF_NOT(soft_equiv(c.begin(), c.end(), sum.begin(), sum.end(), eps));

      c = x;
      global_prod(&c[0], 100);
      FAIL_IF_NOT(soft_equiv(c.begin(), c.end(), prod.begin(), prod.end(), eps));

      c = x;
      global_min(&c[0], 100);
      FAIL_IF_NOT(soft_equiv(c.begin(), c.end(), lmin.begin(), lmin.end(), eps));

      c = x;
      global_max(&c[0], 100);
      FAIL_IF_NOT(soft_equiv(c.begin(), c.end(), lmax.begin(), lmax.end(), eps));
    }
  }
#ifndef CRAYPE_CCE
  { // T = long double

    // make a vector of long doubles
    vector<long double> x(100);
    vector<long double> prod(100, static_cast<long double>(1.0));
    vector<long double> sum(100, static_cast<long double>(0.0));
    vector<long double> lmin(100, static_cast<long double>(0.0));
    vector<long double> lmax(100, static_cast<long double>(0.0));

    // fill it
    for (int i = 0; i < 100; i++) {
      x[i] = static_cast<long double>(rtt_c4::node()) + static_cast<long double>(0.11F);
      for (int j = 0; j < rtt_c4::nodes(); j++) {
        sum[i] += (static_cast<long double>(j) + static_cast<long double>(0.11F));
        prod[i] *= (static_cast<long double>(j) + static_cast<long double>(0.11F));
      }
      lmin[i] = static_cast<long double>(0.11F);
      lmax[i] = static_cast<long double>(rtt_c4::nodes()) + static_cast<long double>(0.11F - 1.0F);
    }

    vector<long double> c;
    auto const eps = static_cast<long double>(1.0e-6F);

    {
      c = x;
      global_sum(&c[0], 100);
      FAIL_IF_NOT(soft_equiv(c.begin(), c.end(), sum.begin(), sum.end(), eps));

      c = x;
      global_sum(&c[0], static_cast<size_t>(100));
      FAIL_IF_NOT(soft_equiv(c.begin(), c.end(), sum.begin(), sum.end(), eps));

      c = x;
      global_prod(&c[0], 100);
      FAIL_IF_NOT(soft_equiv(c.begin(), c.end(), prod.begin(), prod.end(), eps));

      c = x;
      global_min(&c[0], 100);
      FAIL_IF_NOT(soft_equiv(c.begin(), c.end(), lmin.begin(), lmin.end(), eps));

      c = x;
      global_max(&c[0], 100);
      FAIL_IF_NOT(soft_equiv(c.begin(), c.end(), lmax.begin(), lmax.end(), eps));
    }
  }
#endif
  { // T = int

    // make a vector of ints
    vector<int> x(100);
    vector<int> prod(100, 1);
    vector<int> sum(100, 0);
    vector<int> lmin(100, 0);
    vector<int> lmax(100, 0);

    // fill it
    for (int i = 0; i < 100; i++) {
      x[i] = rtt_c4::node() + 1;
      for (int j = 0; j < rtt_c4::nodes(); j++) {
        sum[i] += (j + 1);
        prod[i] *= (j + 1);
      }
      lmin[i] = 1;
      lmax[i] = rtt_c4::nodes();
    }

    vector<int> c;

    {
      c = x;
      global_sum(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, sum));

      c = x;
      global_sum(&c[0], static_cast<size_t>(100));
      FAIL_IF_NOT(std::ranges::equal(c, sum));

      c = x;
      global_prod(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, prod));

      c = x;
      global_min(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, lmin));

      c = x;
      global_max(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, lmax));
    }
  }
  { // T = unsigned int

    // make a vector of unsigned ints
    vector<unsigned int> x(100);
    vector<unsigned int> prod(100, 1);
    vector<unsigned int> sum(100, 0);
    vector<unsigned int> lmin(100, 0);
    vector<unsigned int> lmax(100, 0);

    // fill it
    for (int i = 0; i < 100; i++) {
      x[i] = static_cast<unsigned int>(rtt_c4::node()) + 1;
      for (int j = 0; j < rtt_c4::nodes(); j++) {
        sum[i] += (static_cast<unsigned int>(j) + 1);
        prod[i] *= (static_cast<unsigned int>(j) + 1);
      }
      lmin[i] = 1;
      lmax[i] = static_cast<unsigned int>(rtt_c4::nodes());
    }

    vector<unsigned int> c;

    {
      c = x;
      global_sum(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, sum));

      c = x;
      global_sum(&c[0], static_cast<size_t>(100));
      FAIL_IF_NOT(std::ranges::equal(c, sum));

      c = x;
      global_prod(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, prod));

      c = x;
      global_min(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, lmin));

      c = x;
      global_max(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, lmax));
    }
  }
  { // T = unsigned long

    // make a vector of unsigned longs
    vector<unsigned long> x(100);
    vector<unsigned long> prod(100, 1);
    vector<unsigned long> sum(100, 0);
    vector<unsigned long> lmin(100, 0);
    vector<unsigned long> lmax(100, 0);

    // fill it
    for (int i = 0; i < 100; i++) {
      x[i] = static_cast<unsigned long>(rtt_c4::node()) + 1;
      for (int j = 0; j < rtt_c4::nodes(); j++) {
        sum[i] += (static_cast<unsigned long>(j) + 1);
        prod[i] *= (static_cast<unsigned long>(j) + 1);
      }
      lmin[i] = 1;
      lmax[i] = static_cast<unsigned long>(rtt_c4::nodes());
    }

    vector<unsigned long> c;

    {
      c = x;
      global_sum(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, sum));

      c = x;
      global_sum(&c[0], static_cast<size_t>(100));
      FAIL_IF_NOT(std::ranges::equal(c, sum));

      c = x;
      global_prod(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, prod));

      c = x;
      global_min(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, lmin));

      c = x;
      global_max(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, lmax));
    }
  }
  { // T = unsigned short

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#endif

    // make a vector of unsigned shorts
    vector<unsigned short> x(100);
    vector<unsigned short> prod(100, 1);
    vector<unsigned short> sum(100, 0);
    vector<unsigned short> lmin(100, 0);
    vector<unsigned short> lmax(100, 0);

    // fill it
    for (int i = 0; i < 100; i++) {
      x[i] = static_cast<unsigned short>(rtt_c4::node() + 1);
      for (int j = 0; j < rtt_c4::nodes(); j++) {
        sum[i] += (static_cast<unsigned short>(j + 1));
        prod[i] *= (static_cast<unsigned short>(j + 1));
      }
      lmin[i] = 1;
      lmax[i] = static_cast<unsigned short>(rtt_c4::nodes());
    }

    vector<unsigned short> c;

    {
      c = x;
      global_sum(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, sum));

      c = x;
      global_sum(&c[0], static_cast<size_t>(100));
      FAIL_IF_NOT(std::ranges::equal(c, sum));

      c = x;
      global_prod(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, prod));

      c = x;
      global_min(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, lmin));

      c = x;
      global_max(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, lmax));
    }

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
  }
  { // T = unsigned long long

    // make a vector of unsigned long longs
    vector<unsigned long long> x(100);
    vector<unsigned long long> prod(100, 1);
    vector<unsigned long long> sum(100, 0);
    vector<unsigned long long> lmin(100, 0);
    vector<unsigned long long> lmax(100, 0);

    // fill it
    for (int i = 0; i < 100; i++) {
      x[i] = static_cast<unsigned long long>(rtt_c4::node()) + 1;
      for (int j = 0; j < rtt_c4::nodes(); j++) {
        sum[i] += (static_cast<unsigned long long>(j) + 1);
        prod[i] *= (static_cast<unsigned long long>(j) + 1);
      }
      lmin[i] = 1;
      lmax[i] = static_cast<unsigned long long>(rtt_c4::nodes());
    }

    vector<unsigned long long> c;

    {
      c = x;
      global_sum(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, sum));

      c = x;
      global_sum(&c[0], static_cast<size_t>(100));
      FAIL_IF_NOT(std::ranges::equal(c, sum));

      c = x;
      global_prod(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, prod));

      c = x;
      global_min(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, lmin));

      c = x;
      global_max(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, lmax));
    }
  }
  { // T = long long

    // make a vector of long longs
    vector<long long> x(100);
    vector<long long> prod(100, 1);
    vector<long long> sum(100, 0);
    vector<long long> lmin(100, 0);
    vector<long long> lmax(100, 0);

    // fill it
    for (int i = 0; i < 100; i++) {
      x[i] = static_cast<long long>(rtt_c4::node()) + 1;
      for (int j = 0; j < rtt_c4::nodes(); j++) {
        sum[i] += (static_cast<long long>(j) + 1);
        prod[i] *= (static_cast<long long>(j) + 1);
      }
      lmin[i] = 1;
      lmax[i] = static_cast<long long>(rtt_c4::nodes());
    }

    vector<long long> c;

    {
      c = x;
      global_sum(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, sum));

      c = x;
      global_sum(&c[0], static_cast<size_t>(100));
      FAIL_IF_NOT(std::ranges::equal(c, sum));

      c = x;
      global_prod(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, prod));

      c = x;
      global_min(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, lmin));

      c = x;
      global_max(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, lmax));
    }
  }
  { // T = long

    // make a vector of longs
    vector<long> x(100);
    vector<long> prod(100, 1);
    vector<long> sum(100, 0);
    vector<long> lmin(100, 0);
    vector<long> lmax(100, 0);

    // fill it
    for (int i = 0; i < 100; i++) {
      x[i] = static_cast<long>(rtt_c4::node()) + 1;
      for (int j = 0; j < rtt_c4::nodes(); j++) {
        sum[i] += (static_cast<long>(j) + 1);
        prod[i] *= (static_cast<long>(j) + 1);
      }
      lmin[i] = 1;
      lmax[i] = static_cast<long>(rtt_c4::nodes());
    }

    vector<long> c;

    {
      c = x;
      global_sum(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, sum));

      c = x;
      global_sum(&c[0], static_cast<size_t>(100));
      FAIL_IF_NOT(std::ranges::equal(c, sum));

      c = x;
      global_prod(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, prod));

      c = x;
      global_min(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, lmin));

      c = x;
      global_max(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, lmax));
    }
  }
  { // T = short

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#endif

    // make a vector of shorts
    vector<short> x(100);
    vector<short> prod(100, 1);
    vector<short> sum(100, 0);
    vector<short> lmin(100, 0);
    vector<short> lmax(100, 0);

    // fill it
    for (int i = 0; i < 100; i++) {
      x[i] = static_cast<short>(rtt_c4::node()) + 1;
      for (int j = 0; j < rtt_c4::nodes(); j++) {
        sum[i] += (static_cast<short>(j) + 1);
        prod[i] *= (static_cast<short>(j) + 1);
      }
      lmin[i] = 1;
      lmax[i] = static_cast<short>(rtt_c4::nodes());
    }

    vector<short> c;

    {
      c = x;
      global_sum(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, sum));

      c = x;
      global_sum(&c[0], static_cast<size_t>(100));
      FAIL_IF_NOT(std::ranges::equal(c, sum));

      c = x;
      global_prod(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, prod));

      c = x;
      global_min(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, lmin));

      c = x;
      global_max(&c[0], 100);
      FAIL_IF_NOT(std::ranges::equal(c, lmax));
    }

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
  }

  if (ut.numFails == 0)
    PASSMSG("Array reductions ok.");
  return;
}

//------------------------------------------------------------------------------------------------//
void test_prefix_sum(rtt_dsxx::UnitTest &ut) {

  // Calculate prefix sums on rank ID with MPI call and by hand and compare the output. The prefix
  // sum on a node includes all previous node's value and the value of the current node

  // test ints
  int const xint = rtt_c4::node();
  int const xint_prefix_sum = prefix_sum(xint);

  int int_answer = 0;
  for (int i = 0; i < rtt_c4::nodes(); i++) {
    if (i < rtt_c4::node())
      int_answer += i + 1;
  }

  std::cout << "int: Prefix sum on this node: " << xint_prefix_sum;
  std::cout << " Answer: " << int_answer << '\n';

  FAIL_IF_NOT(xint_prefix_sum == int_answer);

  // test unsigned ints (start at max of signed int)
  uint32_t xuint = rtt_c4::node();
  if (rtt_c4::node() == 0)
    xuint = std::numeric_limits<int>::max();
  uint32_t const xuint_prefix_sum = prefix_sum(xuint);

  uint32_t uint_answer = std::numeric_limits<int>::max();
  for (int i = 0; i < rtt_c4::nodes(); i++) {
    if (i < rtt_c4::node())
      uint_answer += i + 1;
  }

  std::cout << "uint32_t: Prefix sum on this node: " << xuint_prefix_sum
            << " Answer: " << uint_answer << '\n';

  FAIL_IF_NOT(xuint_prefix_sum == uint_answer);

  // test longs
  long const x_long = rtt_c4::node() + 1000;
  long const x_long_prefix_sum = prefix_sum(x_long);

  long long_answer = 0;
  for (int i = 0; i < rtt_c4::nodes(); i++) {
    if (i <= rtt_c4::node() || i == 0)
      long_answer += i + 1000;
  }

  std::cout << "long: Prefix sum on this node: " << x_long_prefix_sum << " Answer: " << long_answer
            << '\n';

  FAIL_IF_NOT(x_long_prefix_sum == long_answer);

  // test long long
  long long const xl_long = rtt_c4::node() + 1000L;
  long long const xl_long_prefix_sum = prefix_sum(xl_long);

  long long llong_answer = 0;
  for (auto i = 0; i < rtt_c4::nodes(); i++) {
    if (i <= rtt_c4::node() || i == 0)
      llong_answer += i + 1000L;
  }

  std::cout << "long long: Prefix sum on this node: " << xl_long_prefix_sum
            << " Answer: " << llong_answer << '\n';

  FAIL_IF_NOT(xl_long_prefix_sum == llong_answer);

  // test unsigned longs (start at max of unsigned int)
  uint64_t xu_long = rtt_c4::node();
  if (rtt_c4::node() == 0)
    xu_long = std::numeric_limits<uint32_t>::max();
  uint64_t const xu_long_prefix_sum = prefix_sum(xu_long);

  uint64_t ulong_answer = std::numeric_limits<uint32_t>::max();
  for (int i = 0; i < rtt_c4::nodes(); i++) {
    if (i < rtt_c4::node())
      ulong_answer += i + 1;
  }

  std::cout << "uint64_t: Prefix sum on this node: " << xu_long_prefix_sum
            << " Answer: " << ulong_answer << '\n';

  FAIL_IF_NOT(xu_long_prefix_sum == ulong_answer);

  // test floats
  auto x_float = static_cast<float>(rtt_c4::node() + 0.01);
  float const x_float_prefix_sum = prefix_sum(x_float);

  float float_answer = 0.0;
  for (int i = 0; i < rtt_c4::nodes(); i++) {
    if (i <= rtt_c4::node() || i == 0)
      float_answer += static_cast<float>(i + 0.01);
  }

  std::cout << "float: Prefix sum on this node: " << x_float_prefix_sum
            << " Answer: " << float_answer << '\n';

  FAIL_IF_NOT(soft_equiv(x_float_prefix_sum, float_answer));

  // test doubles
  double const xdbl = static_cast<double>(rtt_c4::node()) + 1.0e-9;
  double const xdbl_prefix_sum = prefix_sum(xdbl);

  double dbl_answer = 0.0;
  for (int i = 0; i < rtt_c4::nodes(); i++) {
    if (i <= rtt_c4::node() || i == 0)
      dbl_answer += static_cast<double>(i) + 1.0e-9;
  }

  std::cout.precision(16);
  std::cout << "double: Prefix sum on this node: " << xdbl_prefix_sum << " Answer: " << dbl_answer
            << '\n';

  FAIL_IF_NOT(soft_equiv(xdbl_prefix_sum, dbl_answer));

  if (ut.numFails == 0)
    PASSMSG("Prefix sum ok.");
  return;
}

//------------------------------------------------------------------------------------------------//
void test_array_prefix_sum(rtt_dsxx::UnitTest &ut) {

  // Calculate prefix sums on rank ID with MPI call and by hand and compare the output. The prefix
  // sum on a node includes all previous node's value and the value of the current node

  const int array_size = 12;

  // test ints
  vector<int> xint(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i)
    xint[i] = rtt_c4::node() * 10 + i;

  prefix_sum(&xint[0], array_size);

  vector<int> int_answer(array_size, 0);
  for (int i = 0; i < array_size; ++i) {
    for (int r = 0; r < rtt_c4::nodes(); ++r) {
      if (r <= rtt_c4::node())
        int_answer[i] += r * 10 + i;
    }
  }

  for (uint32_t i = 0; i < xint.size(); ++i) {
    std::cout << "int: Prefix sum on this node: " << xint[i] << " Answer: " << int_answer[i]
              << '\n';
    FAIL_IF_NOT(xint[i] == int_answer[i]);
  }

  // test unsigned ints (use the maximum int value to make sure all types are handled correctly in
  // the calls)
  vector<uint32_t> xuint(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i)
    xuint[i] = static_cast<uint32_t>(std::numeric_limits<int>::max()) + rtt_c4::node() * 10 + i;

  prefix_sum(&xuint[0], array_size);

  vector<uint32_t> uint_answer(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i) {
    for (int32_t r = 0; r < rtt_c4::nodes(); ++r) {
      if (r <= rtt_c4::node())
        uint_answer[i] += static_cast<uint32_t>(std::numeric_limits<int>::max()) + r * 10 + i;
    }
  }

  for (uint32_t i = 0; i < xuint.size(); ++i) {
    std::cout << "uint32_t: Prefix sum on this node: " << xuint[i] << " Answer: " << uint_answer[i]
              << '\n';
    FAIL_IF_NOT(xuint[i] == uint_answer[i]);
  }

  // test long ints (use the maximum uint32_t value to make sure all types are handled correctly in
  // the calls)
  vector<int64_t> x_long(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i)
    x_long[i] = std::numeric_limits<uint32_t>::max() + rtt_c4::node() * 10 + i;

  prefix_sum(&x_long[0], array_size);

  vector<int64_t> long_answer(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i) {
    for (int32_t r = 0; r < rtt_c4::nodes(); ++r) {
      if (r <= rtt_c4::node())
        long_answer[i] += std::numeric_limits<uint32_t>::max() + r * 10 + i;
    }
  }

  for (uint32_t i = 0; i < x_long.size(); ++i) {
    std::cout << "int64_t: Prefix sum on this node: " << x_long[i] << " Answer: " << long_answer[i]
              << '\n';
    FAIL_IF_NOT(x_long[i] == long_answer[i]);
  }

  // test unsigned long ints (use the maximum int64_t value to make sure all types are handled
  // correctly in the calls)
  vector<uint64_t> xu_long(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i)
    xu_long[i] =
        static_cast<uint64_t>(std::numeric_limits<int64_t>::max()) + rtt_c4::node() * 10 + i;

  prefix_sum(&xu_long[0], array_size);

  vector<uint64_t> ulong_answer(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i) {
    for (int32_t r = 0; r < rtt_c4::nodes(); ++r) {
      if (r <= rtt_c4::node())
        ulong_answer[i] += static_cast<uint64_t>(std::numeric_limits<int64_t>::max()) + r * 10 + i;
    }
  }

  for (uint32_t i = 0; i < xu_long.size(); ++i) {
    std::cout << "uint64_t: Prefix sum on this node: " << xu_long[i]
              << " Answer: " << ulong_answer[i] << '\n';
    FAIL_IF_NOT(xu_long[i] == ulong_answer[i]);
  }

  // test floats
  vector<float> x_float(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i)
    x_float[i] = static_cast<float>(rtt_c4::node() * 9.99 + i);

  prefix_sum(&x_float[0], array_size);

  vector<float> float_answer(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i) {
    for (int32_t r = 0; r < rtt_c4::nodes(); ++r) {
      if (r <= rtt_c4::node())
        float_answer[i] += static_cast<float>(r * 9.99 + i);
    }
  }

  // comparison between floats after operations needs soft_equiv with a loose tolerance
  for (uint32_t i = 0; i < x_float.size(); ++i) {
    std::cout << "float: Prefix sum on this node: " << x_float[i];
    std::cout << " Answer: " << float_answer[i] << '\n';
    FAIL_IF_NOT(soft_equiv(x_float[i], float_answer[i], float(1.0e-6)));
  }

  // test doubles, try to express precision beyond float to test type handling in function calls
  vector<double> xdouble(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i)
    xdouble[i] = rtt_c4::node() * 9.000000000002 + i;

  prefix_sum(&xdouble[0], array_size);

  vector<double> double_answer(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i) {
    for (int32_t r = 0; r < rtt_c4::nodes(); ++r) {
      if (r <= rtt_c4::node())
        double_answer[i] += r * 9.000000000002 + i;
    }
  }

  for (uint32_t i = 0; i < xdouble.size(); ++i) {
    std::cout << "double: Prefix sum on this node: " << xdouble[i]
              << " Answer: " << double_answer[i] << '\n';
    FAIL_IF_NOT(soft_equiv(xdouble[i], double_answer[i]));
  }

  if (ut.numFails == 0)
    PASSMSG("Array prefix sum ok.");
  return;
}

} // namespace

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    elemental_reduction(ut);
    array_reduction(ut);
    test_prefix_sum(ut);
    test_array_prefix_sum(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstReduction.cc
//------------------------------------------------------------------------------------------------//
