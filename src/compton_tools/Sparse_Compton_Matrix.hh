//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   compton_tools/Sparse_Compton_Matrix.hh
 * \author Andrew Till
 * \date   22 May 2020
 * \brief  Header file for compton interface
 * \note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_compton_tools_Sparse_Compton_Matrix_hh
#define rtt_compton_tools_Sparse_Compton_Matrix_hh

#include <vector>

namespace rtt_compton_tools {

//================================================================================================//
/*!
 * \class Sparse_Compton_Matrix
 *
 * \brief Stores a Compton scattering matrix sparsely and allows matrix multiplication; can store an
 *        arbitrary number of points, where each point represents an M x G x G matrix (moments,
 *        groups)
 *
 * THIS CLASS IS A STUB (not yet fully implemented or used)
*/
//================================================================================================//
class Sparse_Compton_Matrix {

#ifdef 0
public:
  // Construct with zeros and known sizes
  Sparse_Compton_Matrix(size_t num_points, size_t num_groups, size_t num_leg_moments,
                        size_t data_len);

  // Resize data
  void resize_data(size_t data_len);

  // Access internal data for manual filling

  // reference to first groups
  std::vector<size_t> &ref_first_groups() { return first_groups_; }

  // reference to indexes
  std::vector<size_t> &ref_indexes() { return indexes_; }

  // reference to data
  std::vector<double> &ref_data() { return data_; }

  // Multiply against a vector in-place
  // x := matrix * x
  void matvec(std::vector<double> &x, bool zeroth_moment_only = false) const;

  // Accessor functions
  size_t get_num_points() const { return num_points_; }
  size_t get_num_groups() const { return num_groups_; }
  size_t get_num_leg_moments() const { return num_leg_moments_; }
  size_t get_highest_leg_moment() const { return num_leg_moments_ - 1U; }

  // Size checks for valid state
  bool check_class_invariants() const {
    bool all_good = (num_points_ > 0U) && (num_groups_ > 0U) && (num_leg_moments_ > 0U) &&
                    (first_groups_.size() == (num_points_ * num_groups_)) &&
                    (indexes_.size() == (num_points_ * num_groups_ + 1U)) &&
                    (data_.size() >= num_points_ * num_groups_ * num_leg_moments_) && true;
    return all_good;
  }

private:
  size_t num_points_;
  size_t num_groups_;
  size_t num_leg_moments_;

  // sparse data storage

  // first group-to with nonzero value
  // 1D array of [points, group-from]
  std::vector<size_t> first_groups_;

  // cumulative sum of row offsets into data_ and derivs_
  // 1D array of [points, group-from]
  std::vector<size_t> indexes_;

  // csk data
  // 1D array of [moment, point, group-from, group-to]
  std::vector<double> data_;
#endif
};

} // namespace rtt_compton_tools

#endif // rtt_compton_tools_Sparse_Compton_Matrix_hh

//------------------------------------------------------------------------------------------------//
// End compton_tools/Sparse_Compton_Matrix.hh
//------------------------------------------------------------------------------------------------//
