//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   viz/Viz_Traits.hh
 * \author Thomas M. Evans
 * \date   Fri Jan 21 17:10:54 2000
 * \brief  Viz_Traits header file.
 * \note   Copyright (C) 2014-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_viz_Viz_Traits_hh
#define rtt_viz_Viz_Traits_hh

#include "ds++/Assert.hh"
#include <vector>

namespace rtt_viz {

//================================================================================================//
/*!
 * \class Viz_Traits
 *
 * \brief Traits that are used by the rtt_viz package.
 *
 * These traits provide a common way to access 2D-styles arrays/fields in the viz package.
 * Essentially, they allow vector<vector<T> > types to access data using (i,j) operator overloading.
 * There is a general field templated type class; specializations exist for vector<vector>.  Other
 * specializations can be added as needed.
 *
 * The generalized class requires the Field Type (FT) template argument to have the following
 * services:
 *
 * \arg operator()(int i, int j) where the range is [0:N-1, 0:N-1];
 * \arg nrows() returns the number of rows (i index);
 * \arg ncols(int row) returns the number of columns in row (j index);
 * \arg FT::value_type defined to the type returned by the field (int, double, etc).
 */
//================================================================================================//

template <typename FT> class Viz_Traits {
private:
  //! Reference to the field.
  const FT &field;

public:
  //! Constructor.
  explicit Viz_Traits(const FT &field_in) : field(field_in) {}

  //! Overloaded operator().
  typename FT::value_type operator()(size_t i, size_t j) const { return field(i, j); }

  // Overloaded operator().
  typename FT::value_type operator()(size_t i, size_t j, size_t /*k*/) const { return field(i, j); }

  //! Row size accessor.
  size_t nrows() const { return field.nrows(); }

  //! Column size accessor.
  size_t ncols(size_t row) const { return field.ncols(row); }

  // Depth size accessor.
  size_t ndpth(size_t /*row*/, size_t /*col*/) const { return 0; }
};

//------------------------------------------------------------------------------------------------//
// Specialization for std::vector<std::vector>

template <typename T> class Viz_Traits<std::vector<std::vector<T>>> {
public:
  // Type traits
  using elementType = T;

private:
  // Reference to vector<vector> field.
  const std::vector<std::vector<T>> &field;

public:
  // Constructor.
  explicit Viz_Traits(const std::vector<std::vector<T>> &fin) : field(fin) {
    // Nothing to do here
  }

  // Overloaded operator().
  T operator()(size_t i, size_t j) const {
    Require(i < field.size());
    Require(j < field[i].size());
    return field[i][j];
  }

  // Overloaded operator().
  T operator()(size_t i, size_t j, size_t /*k*/) const {
    Require(i < field.size());
    Require(j < field[i].size());
    return field[i][j];
  }

  // Row size accessor.
  size_t nrows() const { return field.size(); }

  // Column size accessor.
  size_t ncols(size_t row) const {
    Require(row < field.size());
    return field[row].size();
  }

  // Depth size accessor.
  size_t ndpth(size_t /*row*/, size_t /*col*/) const { return 0; }
};

//------------------------------------------------------------------------------------------------//
// Specialization for std::vector<std::vector<std::vector>>

template <typename T> class Viz_Traits<std::vector<std::vector<std::vector<T>>>> {
public:
  // Type traits
  using elementType = T;

private:
  // Reference to vector<vector> field.
  const std::vector<std::vector<std::vector<T>>> &field;

public:
  // Constructor.
  explicit Viz_Traits(const std::vector<std::vector<std::vector<T>>> &fin) : field(fin) {
    // Nothing to do here
  }

  // Overloaded operator().
  std::vector<T> operator()(size_t i, size_t j) const {
    Require(i < field.size());
    Require(j < field[i].size());
    return field[i][j];
  }

  // Overloaded operator().
  T operator()(size_t i, size_t j, size_t k) const {
    Require(i < field.size());
    Require(j < field[i].size());
    Require(k < field[i][j].size());
    return field[i][j][k];
  }

  // Row size accessor.
  size_t nrows() const { return field.size(); }

  // Column size accessor.
  size_t ncols(size_t row) const {
    Require(row < field.size());
    return field[row].size();
  }

  // Depth size accessor.
  size_t ndpth(size_t row, size_t col) const {
    Require(row < field.size());
    Require(col < field[row].size());
    return field[row][col].size();
  }
};
} // namespace rtt_viz

#endif // rtt_viz_Viz_Traits_hh

//------------------------------------------------------------------------------------------------//
// end of viz/Viz_Traits.hh
//------------------------------------------------------------------------------------------------//
