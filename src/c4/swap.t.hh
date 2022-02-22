//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/swap.t.hh
 * \author Thomas M. Evans
 * \date   Thu Mar 21 16:56:17 2002
 * \brief  C4 MPI template implementation.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef c4_swap_t_hh
#define c4_swap_t_hh

#include "C4_Functions.hh"
#include "swap.hh"
#include "c4/config.h"
#include "ds++/Assert.hh"

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
// EXCHANGE
//------------------------------------------------------------------------------------------------//

#ifdef C4_MPI

template <typename T>
void determinate_swap(std::vector<unsigned> const &outgoing_pid,
                      std::vector<std::vector<T>> const &outgoing_data,
                      std::vector<unsigned> const &incoming_pid,
                      std::vector<std::vector<T>> &incoming_data, int tag) {
  Require(outgoing_pid.size() == outgoing_data.size());
  Require(incoming_pid.size() == incoming_data.size());
  Require(&outgoing_data != &incoming_data);

  Check(incoming_pid.size() < UINT_MAX);
  auto incoming_processor_count = static_cast<unsigned>(incoming_pid.size());
  Check(outgoing_pid.size() < UINT_MAX);
  auto outgoing_processor_count = static_cast<unsigned>(outgoing_pid.size());

  // This block is a no-op for with-c4=scalar. Ditto when the vectors are of zero-length.
  {

    // Post the asynchronous sends.
    std::vector<C4_Req> outgoing_C4_Req(outgoing_processor_count);
    for (unsigned p = 0; p < outgoing_processor_count; ++p) {
      Check(outgoing_data[p].size() < INT_MAX);
      Check(outgoing_pid[p] < INT_MAX);
      outgoing_C4_Req[p] = rtt_c4::send_async(
          (outgoing_data[p].size() > 0 ? &outgoing_data[p][0] : nullptr),
          static_cast<int>(outgoing_data[p].size()), static_cast<int>(outgoing_pid[p]), tag);
    }

    // Post the asynchronous receives
    std::vector<C4_Req> incoming_C4_Req(incoming_processor_count);
    for (unsigned p = 0; p < incoming_processor_count; ++p) {
      Check(incoming_data[p].size() < INT_MAX);
      Check(incoming_pid[p] < INT_MAX);
      incoming_C4_Req[p] = receive_async(
          (incoming_data[p].size() > 0 ? &incoming_data[p][0] : nullptr),
          static_cast<int>(incoming_data[p].size()), static_cast<int>(incoming_pid[p]), tag);
    }

    // Wait for all the receives to complete.
    wait_all(incoming_processor_count,
             incoming_processor_count > 0 ? &incoming_C4_Req[0] : nullptr);

    // Wait until all the posted sends have completed.
    wait_all(outgoing_processor_count,
             outgoing_processor_count > 0 ? &outgoing_C4_Req[0] : nullptr);
  }
  return;
}

//------------------------------------------------------------------------------------------------//
template <typename T>
void determinate_swap(std::vector<std::vector<T>> const &outgoing_data,
                      std::vector<std::vector<T>> &incoming_data, int tag) {
  Require(static_cast<int>(outgoing_data.size()) == rtt_c4::nodes());
  Require(static_cast<int>(incoming_data.size()) == rtt_c4::nodes());
  Require(&outgoing_data != &incoming_data);

  { // This block is a no-op for with-c4=scalar

    unsigned const N = rtt_c4::nodes();

    // Post the asynchronous sends.
    std::vector<C4_Req> outgoing_C4_Req(N);
    for (unsigned p = 0; p < N; ++p) {
      if (outgoing_data[p].size() > 0) {
        Check(outgoing_data[p].size() < INT_MAX);
        Check(p < INT_MAX);
        outgoing_C4_Req[p] =
            rtt_c4::send_async((outgoing_data[p].size() > 0 ? &outgoing_data[p][0] : nullptr),
                               static_cast<int>(outgoing_data[p].size()), static_cast<int>(p), tag);
      }
    }

    // Post the asynchronous receives
    std::vector<C4_Req> incoming_C4_Req(N);
    for (unsigned p = 0; p < N; ++p) {
      if (incoming_data[p].size() > 0) {
        Check(incoming_data[p].size() < INT_MAX);
        Check(p < INT_MAX);
        incoming_C4_Req[p] =
            receive_async((incoming_data[p].size() > 0 ? &incoming_data[p][0] : nullptr),
                          static_cast<int>(incoming_data[p].size()), static_cast<int>(p), tag);
      }
    }

    // Wait for all the receives to complete.
    wait_all(N, &incoming_C4_Req[0]);

    // Wait until all the posted sends have completed.
    wait_all(N, &outgoing_C4_Req[0]);
  }

  return;
}
//------------------------------------------------------------------------------------------------//
template <typename T>
void semideterminate_swap(std::vector<unsigned> const &outgoing_pid,
                          std::vector<std::vector<T>> const &outgoing_data,
                          std::vector<unsigned> const &incoming_pid,
                          std::vector<std::vector<T>> &incoming_data, int tag) {
  Require(outgoing_pid.size() == outgoing_data.size());
  Require(&outgoing_data != &incoming_data);

  Check(incoming_pid.size() < UINT_MAX);
  auto incoming_processor_count = static_cast<unsigned>(incoming_pid.size());
  Check(outgoing_pid.size() < UINT_MAX);
  auto outgoing_processor_count = static_cast<unsigned>(outgoing_pid.size());

  { // This block is a no-op for with-c4=scalar

    // Send the sizing information using determinate_swap
    std::vector<std::vector<unsigned>> outgoing_size(outgoing_pid.size(), std::vector<unsigned>(1));
    std::vector<std::vector<unsigned>> incoming_size(incoming_pid.size(), std::vector<unsigned>(1));

    for (unsigned p = 0; p < outgoing_processor_count; ++p) {
      Check(outgoing_data.size() < UINT_MAX);
      outgoing_size[p][0] = static_cast<unsigned>(outgoing_data[p].size());
    }

    determinate_swap(outgoing_pid, outgoing_size, incoming_pid, incoming_size, tag);

    // Post the asynchronous sends.
    std::vector<C4_Req> outgoing_C4_Req(outgoing_processor_count);
    for (unsigned p = 0; p < outgoing_processor_count; ++p) {
      Check(outgoing_data[p].size() < INT_MAX);
      Check(outgoing_pid[p] < INT_MAX);
      outgoing_C4_Req[p] = rtt_c4::send_async(
          (outgoing_data[p].size() > 0 ? &outgoing_data[p][0] : nullptr),
          static_cast<int>(outgoing_data[p].size()), static_cast<int>(outgoing_pid[p]), tag);
    }

    // Post the asynchronous receives
    incoming_data.resize(incoming_pid.size());
    std::vector<C4_Req> incoming_C4_Req(incoming_processor_count);
    for (unsigned p = 0; p < incoming_processor_count; ++p) {
      incoming_data[p].resize(incoming_size[p][0]);
      Check(incoming_data[p].size() < INT_MAX);
      Check(incoming_pid[p] < INT_MAX);
      incoming_C4_Req[p] = receive_async(
          (incoming_data[p].size() > 0 ? &incoming_data[p][0] : nullptr),
          static_cast<int>(incoming_data[p].size()), static_cast<int>(incoming_pid[p]), tag);
    }

    // Wait for all the receives to complete.
    if (incoming_C4_Req.size() > 0)
      wait_all(incoming_processor_count, &incoming_C4_Req[0]);

    // Wait until all the posted sends have completed.
    if (outgoing_C4_Req.size() > 0)
      wait_all(outgoing_processor_count, &outgoing_C4_Req[0]);
  }

  return;
}

//------------------------------------------------------------------------------------------------//
// These functions do nothing if there is no communicator (C4_SCALAR=1)
//------------------------------------------------------------------------------------------------//
#else

template <typename T>
void determinate_swap(std::vector<unsigned> const & /*outgoing_pid*/,
                      std::vector<std::vector<T>> const & /*outgoing_data*/,
                      std::vector<unsigned> const & /*incoming_pid*/,
                      std::vector<std::vector<T>> & /*incoming_data*/, int /*tag*/) {
  return;
}
template <typename T>
void determinate_swap(std::vector<std::vector<T>> const & /*outgoing_data*/,
                      std::vector<std::vector<T>> & /*incoming_data*/, int /*tag*/) {
  return;
}
template <typename T>
void semideterminate_swap(std::vector<unsigned> const & /*outgoing_pid*/,
                          std::vector<std::vector<T>> const & /*outgoing_data*/,
                          std::vector<unsigned> const & /*incoming_pid*/,
                          std::vector<std::vector<T>> & /*incoming_data*/, int /*tag*/) {
  return;
}
#endif // C4_MPI

} // end namespace rtt_c4

#endif // c4_swap_t_hh

//------------------------------------------------------------------------------------------------//
// end of c4/swap.t.hh
//------------------------------------------------------------------------------------------------//
