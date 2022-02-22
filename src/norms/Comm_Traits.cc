//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   norms/Comm_Traits.cc
 * \author Rob Lowrie
 * \date   Fri Jan 14 12:45:49 2005
 * \brief  Implementations for Comm_Traits.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Comm_Traits.hh"

namespace rtt_norms {

// Index_Labeled send/receive

void Comm_Traits<Index_Labeled>::send(const Index_Labeled &x, const size_t n) {
  Check(n < INT_MAX);
  rtt_c4::send(&(x.index), 1, static_cast<int>(n));
  rtt_c4::send(&(x.processor), 1, static_cast<int>(n));

  size_t const s = x.label.size();
  Check(s < INT_MAX);
  rtt_c4::send(&s, 1, static_cast<int>(n));
  rtt_c4::send(x.label.c_str(), static_cast<int>(s), static_cast<int>(n));
}

void Comm_Traits<Index_Labeled>::receive(Index_Labeled &x, const size_t n) {
  Check(n < INT_MAX);
  rtt_c4::receive(&(x.index), 1, static_cast<int>(n));
  rtt_c4::receive(&(x.processor), 1, static_cast<int>(n));

  size_t s(0);
  rtt_c4::receive(&s, 1, static_cast<int>(n));
  Check(s < INT_MAX);
  x.label.resize(s);
  rtt_c4::receive(&(x.label[0]), static_cast<int>(s), static_cast<int>(n));
}

// Index_Proc send/receive

void Comm_Traits<Index_Proc>::send(const Index_Proc &x, const size_t n) {
  Check(n < INT_MAX);
  rtt_c4::send(&(x.index), 1, static_cast<int>(n));
  rtt_c4::send(&(x.processor), 1, static_cast<int>(n));
}

void Comm_Traits<Index_Proc>::receive(Index_Proc &x, const size_t n) {
  Check(n < INT_MAX);
  rtt_c4::receive(&(x.index), 1, static_cast<int>(n));
  rtt_c4::receive(&(x.processor), 1, static_cast<int>(n));
}

} // end namespace rtt_norms

//------------------------------------------------------------------------------------------------//
// end of norms/Comm_Traits.hh
//------------------------------------------------------------------------------------------------//
