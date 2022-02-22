//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   Norms_Index.hh
 * \author Rob Lowrie
 * \date   Fri Jan 14 13:00:32 2005
 * \brief  Header file for Norms_Index.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_norms_Norms_Index_hh
#define rtt_norms_Norms_Index_hh

#include "Norms_Base.hh"

namespace rtt_norms {

//================================================================================================//
/*!
  \class Norms_Index

  \brief Computes \f$ L_{\infty} \f$, and \f$ L_1 \f$ \f$ L_2 \f$ norms.

  The norms are computed as
  \f[
      L_1(v) = \frac{\sum_i w_i |v_i|}{\sum_i w_i}, \quad
      L_2(v) = \sqrt{\frac{\sum_i w_i v_i^2}{\sum_i w_i}}, \quad
      L_{\infty}(v) = \max_i(|v_i|).
  \f]
  where \f$ i \f$ is the index, \f$ w_i \f$ the weight factor, and
  \f$ v_i \f$ the corresponding value.

  The norms are computed by making consecutive calls to Norms_Index<>::add.
  Each call to Norms_Index<>::add adds a term to the above summations.
  Norms_Index<>::reset() re-initializes the sums to zero.  The member functions
  Norms_Base::L1, Norms_Base::L2, and Norms_Base::Linf compute their values
  on-processor.  In order to compute these norms across processors, the results
  must be accumulated to a single processor, using a call to
  Norms_Index<>::comm.

  Some member functions are documented in the base class Norms_Base.

  \param Index_t The index type for labeling the location of the Linf norm.
*/
//================================================================================================//
template <typename Index_t> class Norms_Index : public Norms_Base {
public:
  //! Expose the template parameter
  using Index = Index_t;

private:
  // DATA

  // index location of max norm
  Index_t d_index_Linf;

public:
  // CREATORS

  Norms_Index();

  // Use default copy ctor, dtor, and assignment.

  // MANIPULATORS

  // Adds v onto the norms with weight factor.
  void add(const double v, const Index_t &index, const double weight = 1.0);

  // Accumulates results to processor \a n.
  void comm(const size_t n = 0);

  // Equality operator.
  bool operator==(const Norms_Index &n) const;

  // ACCESSORS

  /// Returns the index location of the \f$ L_{\infty} \f$ norm.
  Index_t index_Linf() const { return d_index_Linf; }
};

} // namespace rtt_norms

#endif // rtt_norms_Norms_Index_hh

//------------------------------------------------------------------------------------------------//
// end of norms/Norms_Index.hh
//------------------------------------------------------------------------------------------------//
