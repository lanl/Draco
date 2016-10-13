//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Index_Counter.hh
 * \author Mike Buksas
 * \date   Tue Jan 31 16:45:39 2006
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef dsxx_Index_Counter_hh
#define dsxx_Index_Counter_hh

#include "Index_Set.hh"
#include <vector>

namespace rtt_dsxx {

// forward declaration
template <unsigned D, int OFFSET> class Index_Converter;

//===========================================================================//
/*!
 * \class Index_Counter
 * \brief Facilitates iterating over a multi-dimensional range of indices.
 * \sa Index_Counter.cc for detailed descriptions.
 */
/*!
 * \example ds++/test/tstIndex_Counter.cc
 */
//===========================================================================//
template <unsigned D, int OFFSET> class Index_Counter {
public:
  friend class Index_Converter<D, OFFSET>;

  // NESTED CLASSES AND TYPEDEFS

  // CREATORS

  //! Default constructors.
  Index_Counter(const Index_Set<D, OFFSET> &index_set);

  //! Destructor.
  ~Index_Counter() { /* ... */
  }

  // MANIPULATORS

  //! Assignment operator for Index_Counter.
  Index_Counter &operator=(const Index_Counter &rhs);

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
  int get_index() const { return index; }
  //    operator int()  const { return index; }

  // Accessors for the N-indices
  int get_index(unsigned d) const {
    Check(dimension_okay(d));
    return indices[d];
  }

  std::vector<int> get_indices() const {
    return std::vector<int>(indices, indices + D);
  }

  template <typename IT> void get_indices(IT out) const {
    std::copy(indices, indices + D, out);
  }

  bool is_in_range() const { return in_range; }

private:
  // DATA

  const Index_Set<D, OFFSET> &index_set;

  int indices[D];
  int index;
  bool in_range;

  // IMPLEMENTATION

  void increment();
  void decrement();

  bool dimension_okay(size_t d) const { return d < D; }

  // Private constructor for use by Index_Converter
  // kgbudge (091201): Appears to be dead code
  //     Index_Counter(const Index_Set<D,OFFSET>& base,
  //                   const int index, int const * const indices);
};

//---------------------------------------------------------------------------//
/**
 * \brief Construct from an Index_Set object
 *
 */
template <unsigned D, int OFFSET>
Index_Counter<D, OFFSET>::Index_Counter(const Index_Set<D, OFFSET> &converter)
    : index_set(converter), index(OFFSET), in_range(true) {
  for (size_t d = 0; d < D; ++d)
    indices[d] = OFFSET;
}

//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/**
 * \brief Increment the iterator
 *
 */
template <unsigned D, int OFFSET> void Index_Counter<D, OFFSET>::increment() {

  Require(in_range);

  ++indices[0];
  ++index;

  for (size_t d = 0; d < D - 1; ++d) {
    if (indices[d] > index_set.max_of_index(d)) {
      ++indices[d + 1];
      indices[d] = index_set.min_of_index(d);
    } else
      break;
  }

  if (indices[D - 1] > index_set.max_of_index(D - 1)) {
    indices[D - 1] = index_set.min_of_index(D - 1);
    in_range = false;
  }
}

//---------------------------------------------------------------------------//
/**
 * \brief Decrement the iterator
 *
 */
template <unsigned D, int OFFSET> void Index_Counter<D, OFFSET>::decrement() {

  Require(in_range);

  --indices[0];
  --index;

  for (size_t d = 0; d < D - 1; ++d) {
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

//---------------------------------------------------------------------------//
//              end of ds++/Index_Counter.hh
//---------------------------------------------------------------------------//
