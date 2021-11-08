//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/Index_Converter.hh
 * \author Mike Buksas
 * \date   Fri Jan 20 14:51:51 2006
 * \brief  Decleration and Definition of Index_Converter
 * \note   Copyright 2016-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef dsxx_Index_Converter_hh
#define dsxx_Index_Converter_hh

#include "Index_Counter.hh"
#include <array>

namespace rtt_dsxx {

//================================================================================================//
/*!
 * \class Index_Converter
 * \brief Utiltity class for converting one dimension indicies to and from N-dimensional ones.
 *
 * \sa Index_Converter.cc for detailed descriptions.
 *
 * \example ds++/test/tstIndex_Converter.cc
 * Example use of Index_Converter
 */
//================================================================================================//
template <unsigned D, int OFFSET> class Index_Converter : public Index_Set<D, OFFSET> {
public:
  using Base = Index_Set<D, OFFSET>;
  using Counter = Index_Counter<D, OFFSET>;

  //! default constructor
  Index_Converter() = default;

  //! Construct with just a pointer to the sizes
  Index_Converter(const unsigned *dimensions) { set_size(dimensions); }

  //! Construct a with all dimensions equal
  Index_Converter(const unsigned dimension) { set_size(dimension); }

  //! Copy constructor
  Index_Converter(Index_Converter const &rhs)
      : Index_Set<D, OFFSET>(rhs), sub_sizes(rhs.sub_sizes) {}

  //! Destructor
  ~Index_Converter() override = default;

  //! Re-assignment operator
  void set_size(const unsigned *dimensions);

  //! Uniform size re-assignment operator
  void set_size(unsigned size);

  //! Re-implement base function.
  int limit_of_index(unsigned d, bool pos) const { return Base::limit_of_index(d, pos); }

  // ACCESSORS

  //! Convert N-index to 1-index
  template <typename IT> int get_index(IT indices) const;

  //! Convert 1-index to N-index
  std::vector<int> get_indices(int index) const;

  //! Convert 1-index to N-index and store in provided iterator.
  template <typename IT> void get_indices(int index, IT begin) const;

  //! Extract a single N-index from the 1-index
  int get_single_index(int index, unsigned dimension) const;

  //! Get the next index from a 1-index and direction
  int get_next_index(int index, int direction) const;

  //! Get the next index from an Index_Counter and direction
  int get_next_index(const typename Index_Converter<D, OFFSET>::Counter &counter,
                     int direction) const;

  //! Create an iterator over the index set
  Counter counter() const { return Counter(*this); }

private:
  // DATA

  //! Sizes of sub-grids of increasing dimension.
  std::array<unsigned, D> sub_sizes;

  // IMPLEMENTATION

  void compute_sub_sizes();
};

//------------------------------------------------------------------------------------------------//
// Function Definitions
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
//! Resize the index converter object with new dimensions.
template <unsigned D, int OFFSET>
inline void Index_Converter<D, OFFSET>::set_size(const unsigned *dimensions) {
  Base::set_size(dimensions);
  compute_sub_sizes();
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Resize the index converter with a uniform size
 *
 * \arg dimension The new size
 */
template <unsigned D, int OFFSET>
inline void Index_Converter<D, OFFSET>::set_size(unsigned dimension) {
  Base::set_size(dimension);
  compute_sub_sizes();
}

//------------------------------------------------------------------------------------------------//
//! Convert an N-index to a 1-index
template <unsigned D, int OFFSET>
template <typename IT>
int Index_Converter<D, OFFSET>::get_index(IT indices) const {

  Check(Base::indices_in_range(indices));

  int one_index_value = 0;
  int dimension = 0;
  for (unsigned index = 0; index < D; ++index, ++dimension) {
    one_index_value += (indices[index] - OFFSET) * this->sub_sizes[dimension];
  }

  one_index_value += OFFSET;

  Ensure(Base::index_in_range(one_index_value));

  return one_index_value;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Convert a 1-index to an N-index. Store in provided pointer
 *
 * \arg index The index
 * \arg iterator The iterator pointing to the place to store the results.
 */
template <unsigned D, int OFFSET>
template <typename IT>
void Index_Converter<D, OFFSET>::get_indices(int index, IT iter) const {
  Check(Base::index_in_range(index));
  index -= OFFSET;

  for (unsigned d = 0; d <= D - 1; ++d) {
    const unsigned dim_size = Base::get_size(d);
    *(iter++) = index % dim_size + OFFSET;
    // Ensure that conversion of dim_size from unsigned to int is safe.
    Check(dim_size < static_cast<unsigned>(std::numeric_limits<int>::max()));
    index /= static_cast<int>(dim_size);
  }

  Ensure(index == 0);
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Convert a 1-index to an N-index
 *
 * \arg index The 1-index value
 *
 * This function dispatches to the write-in-place version of the function and stores the result in a
 * local int[] array. It then constructs the return vector in the return statement in order to allow
 * the compiler to perform return value optimization (RVO). This can potentially eliminate the
 * creation of a temporary return object.
 */
template <unsigned D, int OFFSET>
std::vector<int> Index_Converter<D, OFFSET>::get_indices(int index) const {
  Check(Base::index_in_range(index));

  static std::array<int, D> indices;

  get_indices(index, indices.data());

  // Construct in return statement for RVO.
  return std::vector<int>(indices.begin(), indices.begin() + D);
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Extract a single N-index from a 1-index
 *
 * \arg index The 1-index
 * \arg dimension The desired index dimension
 */
template <unsigned D, int OFFSET>
int Index_Converter<D, OFFSET>::get_single_index(int index, unsigned dimension) const {

  Check(Base::index_in_range(index));
  Check(Base::dimension_okay(dimension));

  index -= OFFSET;
  index /= this->sub_sizes[dimension];

  return index % Base::get_size(dimension) + OFFSET;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return the next index in a given direction. Return -1 if this direction is outside the
 *        range of indices
 *
 * \arg index     The index in question
 * \arg direction The direction, 1-based numbered (negative,positive) by dimension.
 */
template <unsigned D, int OFFSET>
int Index_Converter<D, OFFSET>::get_next_index(int index, int direction) const {
  Check(Base::index_in_range(index));
  Check(Base::direction_okay(direction));

  --direction;

  unsigned dimension = direction / 2;
  int sign = 2 * (direction % 2) - 1;

  const int sub_index = get_single_index(index, dimension) + sign;

  if (!Base::index_in_range(sub_index, dimension))
    return -1;

  return index + sign * sub_sizes[dimension];
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Return the next index in a given direction. Return -1 if this direction if outside the
 *        range of indices.
 *
 * \arg counter   An Index_Counter pointing to the desired space in the indices.
 * \arg direction The direction, 1-based numbered (negative, positive) by dimension.
 */
template <unsigned D, int OFFSET>
int Index_Converter<D, OFFSET>::get_next_index(
    const typename Index_Converter<D, OFFSET>::Counter &counter, int direction) const {
  --direction;

  unsigned dimension = direction / 2;
  int sign = 2 * (direction % 2) - 1;

  const int index = counter.get_index();
  const int sub_index = counter.get_index(dimension) + sign;

  if (!Base::index_in_range(sub_index, dimension))
    return -1;

  return index + sign * sub_sizes[dimension];
}

//------------------------------------------------------------------------------------------------//
// IMPLEMENTATION ROUTINES
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/**
 * \brief Assign the internal data members.
 *
 * Used once the dimensions array has been filled.
 */
template <unsigned D, int OFFSET> void Index_Converter<D, OFFSET>::compute_sub_sizes() {

  Require(Base::sizes_okay());

  sub_sizes[0] = 1;

  std::array<unsigned, D> const dimensions = Base::get_dimensions();
  Remember(unsigned *end =) std::partial_sum(dimensions.begin(), dimensions.begin() + D - 1,
                                             sub_sizes.data() + 1, std::multiplies<unsigned>());

  Ensure(end == sub_sizes.data() + D);
}

} // end namespace rtt_dsxx

#endif // dsxx_Index_Converter_hh

//------------------------------------------------------------------------------------------------//
// end of ds++/Index_Converter.hh
//------------------------------------------------------------------------------------------------//
