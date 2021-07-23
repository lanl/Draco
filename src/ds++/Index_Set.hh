//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/Index_Set.hh
 * \author Mike Buksas
 * \date   Thu Feb  2 10:01:46 2006
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef dsxx_Index_Set_hh
#define dsxx_Index_Set_hh

#include "Assert.hh"
#include <algorithm>
#include <array>

namespace rtt_dsxx {

//================================================================================================//
/*!
 * \class Index_Set
 * \brief Represents a D-dimensional set if indices.
 * \sa Index_Set.cc for detailed descriptions.
 */
/*!
 * \example ds++/test/tstIndex_Set.cc
 */
//================================================================================================//
template <unsigned D, int OFFSET> class Index_Set {
public:
  // CREATORS

  //! Default constructors.
  Index_Set() = default;

  //! Construct with pointer to sizes
  explicit Index_Set(unsigned const *const dimensions) { set_size(dimensions); }

  //! Construct with all dimensions equal
  explicit Index_Set(const unsigned dimension) { set_size(dimension); }

  //! Destructor
  virtual ~Index_Set() = default;

  //! Comparison operator
  bool operator==(const Index_Set &rhs) const;

  //! Negative comparison operator
  bool operator!=(const Index_Set &rhs) const { return !(*this == rhs); }

  //! Re-assignment operator
  void set_size(unsigned const *const dimensions);

  //! Uniform size re-assignment operator
  void set_size(const unsigned size);

  bool index_in_range(int index) const {
    return (index >= OFFSET) && (index < static_cast<int>(m_array_size) + OFFSET);
  }
  bool index_in_range(int index, unsigned dimension) const;

  template <typename IT> bool indices_in_range(IT indices) const;

  unsigned get_size() const { return m_array_size; }
  int min_of_index() const { return OFFSET; }
  int max_of_index() const { return OFFSET + m_array_size - 1; }
  int limit_of_index(const bool positive) const {
    return positive ? max_of_index() : min_of_index();
  }

  unsigned get_size(const unsigned d) const {
    Check(dimension_okay(d));
    return m_dimensions[d];
  }
  int min_of_index(const unsigned Remember(d)) const {
    Check(dimension_okay(d));
    return OFFSET;
  }
  int max_of_index(const unsigned d) const {
    Check(dimension_okay(d));
    return OFFSET + m_dimensions[d] - 1;
  }
  int limit_of_index(const unsigned d, const bool positive) const {
    return positive ? max_of_index(d) : min_of_index(d);
  }

  static bool direction_okay(const size_t d) { return (d > 0) && (d <= 2 * D); }
  static bool dimension_okay(const size_t d) { return d < D; }

private:
  void compute_size();

  unsigned m_array_size{0};             //!< Sizes of the whole index range
  std::array<unsigned, D> m_dimensions; //!< Sizes of each dimension

protected:
  // Make sure the index sizes are all positive when creating or resizing:
  bool sizes_okay() const {
    return (std::find(m_dimensions.begin(), m_dimensions.begin() + D, 0u) ==
            m_dimensions.begin() + D);
  }

  // Allow derived classes const access to the dimensions.
  std::array<unsigned, D> const get_dimensions() const { return m_dimensions; }
};

} // end namespace rtt_dsxx

#include "Index_Set.i.hh"

#endif // dsxx_Index_Set_hh

//------------------------------------------------------------------------------------------------//
// end of ds++/Index_Set.hh
//------------------------------------------------------------------------------------------------//
