//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   compton_tools/Sparse_Compton_Matrix.cc
 * \author Andrew Till
 * \date   11 May 2020
 * \brief  Implementation file for compton interface
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "compton_tools/Sparse_Compton_Matrix.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include <iostream>

namespace rtt_compton_tools {

#ifdef 0
/*!
 * \brief Data container class for a sparse Compton matrix representation. Not yet implemented.
 */
Sparse_Compton_Matrix::Sparse_Compton_Matrix(size_t num_points, size_t num_groups,
                                             size_t num_leg_moments, size_t data_len)
    : num_points_(num_points), num_groups_(num_groups), num_leg_moments_(num_leg_moments),
      first_groups_(num_points * num_groups, 0U), indexes_(num_points * num_groups + 1U, 0U),
      data_(data_len, 0.0) {
  Require(num_points >= 1U);
  Require(num_groups >= 1U);
  Require(num_leg_moments >= 1U);
  Require(data_len >= (num_points * num_groups * num_leg_moments));

  std::cout << "SCM constructor not yet fully implemented\n";
  Ensure(check_class_invariants());
}

void Sparse_Compton_Matrix::resize_data(size_t data_len) {
  Require(data_len >= (num_points_ * num_groups_ * num_leg_moments_));

  data_.resize(data_len, 0.0);

  Ensure(check_class_invariants());
}

void Sparse_Compton_Matrix::matvec(std::vector<double> &x, bool zeroth_moment_only) const {
  // TODO: Implement!
  std::cout << "SCM matvec not yet implemented\n";
}
#endif

} // namespace rtt_compton_tools

//------------------------------------------------------------------------------------------------//
// End compton_tools/Sparse_Compton_Matrix.cc
//------------------------------------------------------------------------------------------------//
