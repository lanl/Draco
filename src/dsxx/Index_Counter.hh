//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   dsxx/Index_Counter.hh
 * \author Mike Buksas
 * \date   Tue Jan 31 16:45:39 2006
 * \note   Copyright (C) 2024-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef dsxx_Index_Counter_hh
#define dsxx_Index_Counter_hh

#include "Index_Set.hh"
#include <array>
#include <vector>

namespace rtt_dsxx {

// forward declaration
template <unsigned D, int OFFSET> class Index_Converter;

//================================================================================================//
/*!
 * \class Index_Counter
 * \brief Facilitates iterating over a multi-dimensional range of indices.
 * \sa Index_Counter.cc for detailed descriptions.
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
  Index_Counter(Index_Counter const &rhs)
      : index_set(rhs.index_set), indices(rhs.indices), index(rhs.index), in_range(rhs.in_range) {}

  //! Move constructor
  Index_Counter(Index_Counter const &&rhs) noexcept
      : index_set(rhs.index_set), indices(rhs.indices), index(rhs.index), in_range(rhs.in_range) {}

  //! Destructor.
  ~Index_Counter() = default;

  // MANIPULATORS

  //! Assignment operator for Index_Counter.
  Index_Counter &operator=(const Index_Counter &rhs);
  Index_Counter &operator=(Index_Counter &&rhs) noexcept = delete;

  // ACCESSORS

  Index_Counter &operator++() {
    increment();
    return *this;
  }
  Index_Counter &operator--() {
    decrement();
    return *this;
  }

  // Accessors for the 1-index
  [[nodiscard]] int get_index() const { return index; }

  // Accessors for the N-indices
  [[nodiscard]] int get_index(unsigned d) const {
// gcc RWDI builds erroneously flag these lines as accessing array subscript 3.
// We simply suppress the array-bounds warning here.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#endif
    Check(dimension_okay(d));
    return indices[d];
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
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

  const Index_Set<D, OFFSET> &index_set;

  std::array<int, D> indices;
  int index;
  bool in_range{true};

  // IMPLEMENTATION

  void increment();
  void decrement();

  [[nodiscard]] bool dimension_okay(size_t d) const { return d < D; }
};

//------------------------------------------------------------------------------------------------//
//! Construct from an Index_Set object
template <unsigned D, int OFFSET>
Index_Counter<D, OFFSET>::Index_Counter(const Index_Set<D, OFFSET> &converter)
    : index_set(converter), index(OFFSET) {
  for (size_t d = 0; d < D; ++d)
    indices[d] = OFFSET;
}

//------------------------------------------------------------------------------------------------//
// IMPLEMENTATION
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
//! Increment the iterator
template <unsigned D, int OFFSET> void Index_Counter<D, OFFSET>::increment() {

  Require(in_range);

  ++indices[0];
  ++index;

  for (size_t d = 0; d < static_cast<size_t>(D) - 1; ++d) {
    if (indices[d] > index_set.max_of_index(d)) {
      ++indices[d + 1];
      indices[d] = index_set.min_of_index(d);
    } else {
      break;
    }
  }

  if (indices[D - 1] > index_set.max_of_index(D - 1)) {
    indices[D - 1] = index_set.min_of_index(D - 1);
    in_range = false;
  }
}

//------------------------------------------------------------------------------------------------//
//! Decrement the iterator
template <unsigned D, int OFFSET> void Index_Counter<D, OFFSET>::decrement() {

  Require(in_range);

  --indices[0];
  --index;

  for (size_t d = 0; d < static_cast<size_t>(D) - 1; ++d) {
    if (indices[d] < index_set.min_of_index(d)) {
      indices[d] = index_set.max_of_index(d);
      --indices[d + 1];
    } else
      break;
  }

  if (indices[D - 1] < index_set.min_of_index(D - 1)) {
    indices[D - 1] = index_set.max_of_index(D - 1);
    in_range = false;
  }
}

} // end namespace rtt_dsxx

#endif // dsxx_Index_Counter_hh

//------------------------------------------------------------------------------------------------//
// end of dsxx/Index_Counter.hh
//------------------------------------------------------------------------------------------------//
