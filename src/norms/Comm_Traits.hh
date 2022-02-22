//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   norms/Comm_Traits.hh
 * \author Rob Lowrie
 * \date   Fri Jan 14 12:45:49 2005
 * \brief  Header for Comm_Traits.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_norms_Comm_Traits_hh
#define rtt_norms_Comm_Traits_hh

#include "Index_Labeled.hh"
#include "Index_Proc.hh"
#include "c4/C4_Functions.hh"

namespace rtt_norms {

//================================================================================================//
/*!
  \class Comm_Traits
  \brief Communicates the Norms index type.
  
  This class contains the necessary functions used to communicate the index
  type of class Norms.  The following static functions must be provided:

  // Sends \a x to processor \a n.
  void send(const Index_t x, const size_t n)

  // Receives \a x from processor \a n.
  void receive(Index_t &x, const size_t n)
 */
//================================================================================================//

template <typename Index_t> class Comm_Traits {
  // The default implementation of Comm_Traits assumes that Index_t is
  // supported as an argument type to c4's send/receive.

public:
  static void send(const Index_t x, const size_t n) {
    Check(n < INT_MAX);
    rtt_c4::send(&x, 1, static_cast<int>(n));
  }

  static void receive(Index_t &x, const size_t n) {
    Check(n < INT_MAX);
    rtt_c4::receive(&x, 1, static_cast<int>(n));
  }
};

//------------------------------------------------------------------------------------------------//
// SPECIALIZATIONS
//------------------------------------------------------------------------------------------------//

// For Index_Labeled
template <> class Comm_Traits<Index_Labeled> {
public:
  static void send(const Index_Labeled &x, const size_t n);
  static void receive(Index_Labeled &x, const size_t n);
};

// For Index_Proc
template <> class Comm_Traits<Index_Proc> {
public:
  static void send(const Index_Proc &x, const size_t n);
  static void receive(Index_Proc &x, const size_t n);
};

} // end namespace rtt_norms

#endif // rtt_norms_Comm_Traits_hh

//------------------------------------------------------------------------------------------------//
// end of norms/Comm_Traits.hh
//------------------------------------------------------------------------------------------------//
