//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   viz/Viz_Traits.hh
 * \author Thomas M. Evans
 * \date   Fri Jan 21 17:10:54 2000
 * \brief  Viz_Traits header file.
 * \note   Copyright (C) 2014-2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_viz_Viz_Traits_hh
#define rtt_viz_Viz_Traits_hh

#include "ds++/Assert.hh"
#include <set>
#include <vector>

namespace rtt_viz {

// helper for SFINAE
template <bool C, typename T> using enable_if_t = typename std::enable_if<C, T>::type;

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
 *
 * https://stackoverflow.com/questions/61562051/c-sfinae-for-array-vector
 *
 * \tparam FT field type must be a 2D or 3D STL container.  Eg. vector<vector<int>> or
 *          vector<vector<set<float>>>.
 * \tparam elementTypeOuter SFINAE failure if FT is not an STL container.
 * \tparam elementType SFINAE failure if elementTypeOuter is not an STL container.
 */
//================================================================================================//

template <typename FT, typename elementTypeOuter = typename FT::value_type,
          typename elementType = typename elementTypeOuter::value_type>
class Viz_Traits {

  //! Reference to the field.
  const FT &field;

public:
  //! Constructor.
  explicit Viz_Traits(const FT &field_in) : field(field_in) {}

  //! Overloaded operator() with 2 args.
  auto operator()(size_t i, size_t j) const {
    Require(i < field.size());
    Require(j < field[i].size());
    return field[i][j];
  }

  //! Overloaded operator() with 3 args used on 2D data
  template <typename U = elementType, enable_if_t<std::is_pod<U>::value, bool> = true>
  auto operator()(size_t i, size_t j, size_t /*k*/) const {
    return this->operator()(i, j);
  }
  /*
  * \brief Overloaded operator() with 3 args used on 3D data.
  * \bug If the innermost container is \c std::set, then the [k] operator is invalid.  You cannot
  *      access the elements of a \c std::set with an index.
  */
  template <typename U = elementType, enable_if_t<!std::is_pod<U>::value, bool> = true>
  auto operator()(size_t i, size_t j, size_t k) const {
    Require(i < field.size());
    Require(j < field[i].size());
    Require(k < field[i][j].size());
    return field[i][j][k];
  }

  //! Row size accessor.
  size_t nrows() const { return field.size(); }

  //! Column size accessor.
  size_t ncols(size_t row) const {
    Require(row < field.size());
    return field[row].size();
  }

  //! Depth size accessor for 2D fields
  template <typename U = elementType, enable_if_t<std::is_pod<U>::value, bool> = true>
  size_t ndpth(size_t /*row*/, size_t /*col*/) const {
    return 0;
  }

  //! Depth size accessor for 3D fields
  template <typename U = elementType, enable_if_t<!std::is_pod<U>::value, bool> = true>
  size_t ndpth(size_t row, size_t col) const {
    Require(row < field.size());
    Require(col < field[row].size());
    return field[row][col].size();
  }

}; // class Viz_Trais

} // namespace rtt_viz

#endif // rtt_viz_Viz_Traits_hh

//------------------------------------------------------------------------------------------------//
// end of viz/Viz_Traits.hh
//------------------------------------------------------------------------------------------------//
