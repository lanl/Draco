//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   Norms_Index.t.hh
 * \author Rob Lowrie
 * \date   Fri Jan 14 13:00:47 2005
 * \brief  Implemention for Norms_Index class.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Comm_Traits.hh"
#include "Norms_Index.hh"

namespace rtt_norms {

//------------------------------------------------------------------------------------------------//
//! Default constructor.
template <typename Index_t>
Norms_Index<Index_t>::Norms_Index() : Norms_Base(), d_index_Linf(Index_t()) {
  /* empty */
}

//------------------------------------------------------------------------------------------------//
/*!
  \brief Adds to the norm values.

  \param v      The value to be added to the norms.
  \param index  The location of \a v.
  \param weight The weight factor for \a v.
*/
//------------------------------------------------------------------------------------------------//
template <typename Index_t>
void Norms_Index<Index_t>::add(const double v, const Index_t &index, const double weight) {
  double vabs = std::fabs(v);

  d_sum_weights += weight;
  d_sum_L1 += weight * vabs;
  d_sum_L2 += weight * v * v;

  if (vabs >= d_Linf) {
    d_Linf = vabs;
    d_index_Linf = index;
  }
}

//------------------------------------------------------------------------------------------------//
//! Re-initializes the norm values.
// template <typename Index_t> void Norms_Index<Index_t>::reset() {
//   Norms_Base::reset();
//   d_index_Linf = Index_t();
// }

//------------------------------------------------------------------------------------------------//
/*!
  \brief Accumulates results to proc \a n.

  After calling this function, processor \a n contains the norms over all
  processors.  All processors other than \a n still contain their same
  norm values.

  \param n Processor on which norms are summed.
*/
//------------------------------------------------------------------------------------------------//
template <typename Index_t> void Norms_Index<Index_t>::comm(const size_t n) {
  const size_t num_nodes = rtt_c4::nodes();

  if (num_nodes == 1)
    return;

  const size_t node = rtt_c4::node();

  if (node == n) {
    // Accumulate the results onto this node.

    Index_t pe_index(0); // temporary for index

    for (size_t i = 0; i < num_nodes; ++i) {
      if (i != n) {
        double x(0);
        Check(i < INT_MAX);
        rtt_c4::receive(&x, 1, static_cast<int>(i));
        d_sum_L1 += x;
        rtt_c4::receive(&x, 1, static_cast<int>(i));
        d_sum_L2 += x;
        rtt_c4::receive(&x, 1, static_cast<int>(i));
        d_sum_weights += x;
        rtt_c4::receive(&x, 1, static_cast<int>(i));
        Comm_Traits<Index_t>::receive(pe_index, i);
        if (x > d_Linf) {
          d_Linf = x;
          d_index_Linf = pe_index;
        }
      }
    }
  } else {
    // Send this proc's result back to node n.

    Check(n < INT_MAX);
    rtt_c4::send(&d_sum_L1, 1, static_cast<int>(n));
    rtt_c4::send(&d_sum_L2, 1, static_cast<int>(n));
    rtt_c4::send(&d_sum_weights, 1, static_cast<int>(n));
    rtt_c4::send(&d_Linf, 1, static_cast<int>(n));
    Comm_Traits<Index_t>::send(d_index_Linf, n);
  }
}

//------------------------------------------------------------------------------------------------//
//! Equality operator.
template <typename Index_t> bool Norms_Index<Index_t>::operator==(const Norms_Index &n) const {
  bool b = Norms_Base::operator==(n);
  return b && (d_index_Linf == n.d_index_Linf);
}

} // namespace rtt_norms

//------------------------------------------------------------------------------------------------//
// end of norms/Norms_Index.t.hh
//------------------------------------------------------------------------------------------------//
