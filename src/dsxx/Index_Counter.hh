//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   dsxx/Index_Counter.hh
 * \author Mike Buksas
 * \note   Copyright (C) 2024-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef dsxx_Index_Counter_hh
#define dsxx_Index_Counter_hh

#include "dsxx/Index_Set.hh"
#include <array>
#include <vector>

namespace rtt_dsxx {

// forward declaration
template <unsigned D, int OFFSET> class Index_Converter;

//================================================================================================//
/*!
 * \class Index_Counter
 * \brief Utility for iterating over multi-dimensional indices
 *
 * \details The Index_Counter class provides a mechanism for iterating through a multi-dimensional
 *          grid of indices while maintaining both the individual coordinates (multi-dimensional
 *          indices) and a single linear index that represents the position in a flattened array.
 *
 * \tparam D The number of dimensions for the counter
 * \tparam OFFSET The starting index (typically 0 or 1)
 *
 * \section index_counter_features Key Features
 *
 * \subsection index_counter_feature1 Multi-dimensional Iteration
 * Handles the complexities of incrementing and decrementing across multiple dimensions,  similar to
 * how a car's odometer works - when one digit reaches its limit, it rolls over and increments the
 * next digit.
 *
 * \subsection index_counter_feature2 Coordinate Mapping
 * Maintains both the N-dimensional coordinates (in the indices array) and a single linear index
 * (in the index variable), allowing access to data in either representation.
 *
 * \subsection feature3 Range Checking
 * The in_range flag keeps track of whether the iterator is still within the valid range of indices.
 *
 * \section index_counter_usage Usage
 *
 * \code
 * // Create an Index_Set representing a 3D grid of size 10x10x10
 * Index_Set<3, 0> index_set;
 * index_set.set_size(10); // All dimensions of size 10
 *
 * // Create an iterator over this set
 * Index_Counter<3, 0> counter(index_set);
 *
 * // Iterate through all points in the 3D grid
 * while (counter.is_in_range()) {
 *     // Get the flattened index
 *     int linear_idx = counter.get_index();
 *
 *     // Get the 3D coordinates
 *     int x = counter.get_index(0);
 *     int y = counter.get_index(1);
 *     int z = counter.get_index(2);
 *
 *     // Process the point at (x,y,z)
 *     // ...
 *
 *     // Move to the next point
 *     ++counter;
 * }
 * \endcode
 *
 * \note Index_Counter works in tandem with Index_Converter which handles conversion between linear
 *       indices and multi-dimensional coordinates. While Index_Counter focuses on iteration,
 *       Index_Converter specializes in the mathematical conversion between these two
 *       representations.
 *
 * \sa Index_Converter, Index_Set
 *
 * \example dsxx/test/tstIndex_Counter.cc
 */
//================================================================================================//
template <unsigned D, int OFFSET> class Index_Counter {
public:
  friend class Index_Converter<D, OFFSET>;

  // CREATORS

  //! Default constructors.
  explicit Index_Counter(const Index_Set<D, OFFSET> &index_set);

  //! Copy constructor
  Index_Counter(Index_Counter const &rhs) = default;
  //! Move constructor
  Index_Counter(Index_Counter &&rhs) noexcept = default;

  //! Destructor.
  ~Index_Counter() = default;

  // MANIPULATORS

  //! Assignment operator for Index_Counter.
  Index_Counter &operator=(const Index_Counter &rhs);
  Index_Counter &operator=(Index_Counter &&rhs) noexcept = delete;

  // ACCESSORS

  Index_Counter &operator++();
  Index_Counter &operator--();

  // Accessors for the 1-index
  [[nodiscard]] int get_index() const { return index; }

  // Accessors for the N-indices
  [[nodiscard]] int get_index(unsigned d) const {
    Require(d < D);
    return indices[d];
  }

  [[nodiscard]] std::vector<int> get_indices() const {
    return std::vector<int>(indices.begin(), indices.begin() + D);
  }

  template <typename IT> void get_indices(IT out) const {
    std::copy(indices.begin(), indices.begin() + D, out);
  }

  [[nodiscard]] bool is_in_range() const { return in_range; }

private:
  // DATA

  /*! Reference to the Index_Set that defines the range of valid indices for each dimension. Used to
   * determine boundaries when incrementing or decrementing the counter. */
  const Index_Set<D, OFFSET> &index_set;

  /*! Stores the current position as N-dimensional coordinates. Each element indices[d] represents
   * the current position in dimension d. */
  std::array<int, D> indices;

  /*! The linear index corresponding to the current N-dimensional position. This value provides a
   * one-dimensional representation of the current position that can be used to access elements in a
   * flattened array. */
  int index;

  /*! Indicates whether the counter is currently within the valid range of indices. Set to false
   * when the counter has been incremented past the end of the range or decremented past the
   * beginning. */
  bool in_range{true};
};

//------------------------------------------------------------------------------------------------//
//! Construct from an Index_Set object
template <unsigned D, int OFFSET>
Index_Counter<D, OFFSET>::Index_Counter(const Index_Set<D, OFFSET> &index_set_arg)
    : index_set(index_set_arg), index(OFFSET) {
  indices.fill(OFFSET);
}

//------------------------------------------------------------------------------------------------//
// IMPLEMENTATION
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
//! Decrement the iterator
template <unsigned D, int OFFSET> Index_Counter<D, OFFSET> &Index_Counter<D, OFFSET>::operator--() {
  Require(in_range);

  --indices[0];
  --index;

  // Use constexpr to potentially optimize compile-time evaluation
  constexpr size_t last_dim = D - 1;

  // Early check for single dimension case
  if constexpr (D == 1) {
    if (indices[0] < index_set.min_of_index(0)) {
      indices[0] = index_set.max_of_index(0);
      in_range = false;
    }

  } else {

    // Multi-dimensional case
    for (size_t d = 0; d < last_dim; ++d) {
      if (indices[d] < index_set.min_of_index(d)) {
        indices[d] = index_set.max_of_index(d);
        --indices[d + 1];
      } else
        break;
    }

    if (indices[last_dim] < index_set.min_of_index(last_dim)) {
      indices[last_dim] = index_set.max_of_index(last_dim);
      in_range = false;
    }
  }
  return *this;
}

//------------------------------------------------------------------------------------------------//
//! Increment the iterator
template <unsigned D, int OFFSET> Index_Counter<D, OFFSET> &Index_Counter<D, OFFSET>::operator++() {
  Require(in_range);

  ++indices[0];
  ++index;

  // Use constexpr to potentially optimize compile-time evaluation
  constexpr size_t last_dim = D - 1;

  // Early check for single dimension case
  if constexpr (D == 1) {
    if (indices[0] > index_set.max_of_index(0)) {
      indices[0] = index_set.min_of_index(0);
      in_range = false;
    }

  } else {

    // Multi-dimensional case
    for (size_t d = 0; d < last_dim; ++d) {
      if (indices[d] <= index_set.max_of_index(d)) {
        // No carry needed, done with the increment
        break;
      }

      // Handle carry to next dimension
      ++indices[d + 1];
      indices[d] = index_set.min_of_index(d);
    }

    // Check if we've gone beyond the last dimension's range
    if (indices[last_dim] > index_set.max_of_index(last_dim)) {
      indices[last_dim] = index_set.min_of_index(last_dim);
      in_range = false;
    }
  }
  return *this;
}

} // end namespace rtt_dsxx
#endif // dsxx_Index_Counter_hh

//------------------------------------------------------------------------------------------------//
// end of dsxx/Index_Counter.hh
//------------------------------------------------------------------------------------------------//
