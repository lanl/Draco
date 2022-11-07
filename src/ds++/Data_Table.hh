//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file    Data_Table.hh
 * \author  Paul Henning
 * \brief   Declaration of class Data_Table
 * \note    Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//
#ifndef dsxx_Data_Table_hh
#define dsxx_Data_Table_hh

#include "Assert.hh"
#include "ds++/config.h"
#include <vector>

namespace rtt_dsxx {

//================================================================================================//
/*!
 * \class Data_Table
 *
 * Provide const array-style access to an actual array or a scalar.
 *
 * Data_Table provides read-only, DBC-checked, container-like access to a sequential range of memory
 * locations, or a scalar.  This is useful in situations where the amount of data changes depending
 * on compile-time factors, but you always want to access it as an array. Because this class is read
 * only all of it's functions and constructor can be marked constexpr.
 */
//================================================================================================//
template <typename T> class Data_Table {
public:
  using const_iterator = const T *;

public:
  //! copy constructor
  constexpr inline Data_Table(Data_Table const &rhs);
  constexpr inline Data_Table(Data_Table const &&rhs) = delete;
  //! assignment operator
  constexpr Data_Table &operator=(Data_Table const &rhs);
  constexpr Data_Table &operator=(Data_Table const &&rhs) noexcept = delete;
  //! Constructor
  constexpr inline explicit Data_Table(std::vector<T> const &v);
  //! Constructor
  constexpr inline Data_Table(const_iterator const begin, const_iterator const end);
  //! Constructor
  constexpr inline explicit Data_Table(T const &value);
  //! Default constructor
  constexpr inline Data_Table();
  //! Destructor
  inline ~Data_Table() = default;
  //! Access operator
  constexpr inline T const &operator[](const unsigned i) const;
  //! begin iterator
  constexpr inline const_iterator begin() const { return d_begin; }
  //! end iterator
  constexpr inline const_iterator end() const { return d_end - 1; }
  constexpr inline uint64_t size() const { return d_end - d_begin; }
  constexpr inline T const &front() const;
  constexpr inline T const &back() const;
  constexpr inline T *access();

private:
  const_iterator const d_begin;
  const_iterator const d_end; //!< one past the end of the container

  /*! We hold a copy of the scalar to prevent the problems that would arise if you took a pointer to
   *  a function-return temporary. */
  T const d_value;
};

//------------------------------------------------------------------------------------------------//
/*!
 * Copy constructor, but update the pointers to point to the local d_value if they pointed to the
 * d_value in the RHS.
*/
template <typename T>
constexpr Data_Table<T>::Data_Table(Data_Table<T> const &rhs)
    : d_begin(rhs.d_begin), d_end(rhs.d_end), d_value(rhs.d_value) {
  if (rhs.d_begin == &(rhs.d_value)) {
    const_cast<const_iterator &>(d_begin) = &d_value;
    const_cast<const_iterator &>(d_end) = d_begin + 1;
  }
}

//------------------------------------------------------------------------------------------------//
template <typename T> constexpr Data_Table<T> &Data_Table<T>::operator=(Data_Table<T> const &rhs) {
  if (&rhs != this) {
    if (rhs.d_begin == &(rhs.d_value)) {
      const_cast<T &>(d_value) = rhs.d_value;
      const_cast<const_iterator &>(d_begin) = &d_value;
      const_cast<const_iterator &>(d_end) = d_begin + 1;
    } else {
      const_cast<const_iterator &>(d_begin) = rhs.d_begin;
      const_cast<const_iterator &>(d_end) = rhs.d_end;
    }
  }
  return *this;
}

//------------------------------------------------------------------------------------------------//
template <typename T>
constexpr inline Data_Table<T>::Data_Table(const_iterator const begin, const_iterator const end)
    : d_begin(begin), d_end(end + 1), d_value() {
  Require(!(begin > end));
}

//------------------------------------------------------------------------------------------------//
//! Copy the scalar into a local variable, and set the pointers to that copy.
template <typename T>
constexpr inline Data_Table<T>::Data_Table(T const &value)
    : d_begin(&d_value), d_end(&d_value + 1), d_value(value) {}

//------------------------------------------------------------------------------------------------//
template <typename T>
constexpr inline Data_Table<T>::Data_Table() : d_begin(nullptr), d_end(nullptr), d_value() {}

//------------------------------------------------------------------------------------------------//
template <typename T> constexpr inline T const &Data_Table<T>::operator[](unsigned const i) const {
  Require(i < size());
  return d_begin[i];
}

//------------------------------------------------------------------------------------------------//
template <typename T> constexpr inline T const &Data_Table<T>::front() const {
  Require((d_end - d_begin) > 0);
  return *d_begin;
}

//------------------------------------------------------------------------------------------------//
template <typename T> constexpr inline T const &Data_Table<T>::back() const {
  Require((d_end - d_begin) > 0);
  return *(d_end - 1);
}

//------------------------------------------------------------------------------------------------//
template <typename T> constexpr inline T *Data_Table<T>::access() {
  Require((d_end - d_begin) > 0);
  return const_cast<T *>(d_begin);
}

} // namespace rtt_dsxx

#endif // dsxx_Data_Table_hh

//------------------------------------------------------------------------------------------------//
// end of Data_Table.hh
//------------------------------------------------------------------------------------------------//
