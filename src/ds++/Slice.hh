//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/Slice.hh
 * \author Kent Budge
 * \date   Thu Jul  8 08:06:53 2004
 * \brief  Definition of Slice template
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef container_Slice_hh
#define container_Slice_hh

#include "Assert.hh"
#include <iterator>

namespace rtt_dsxx {

//================================================================================================//
/*!
 * \class Slice
 * \brief A reference to a Slice from a random access sequence.
 *
 * This class has reference semantics to the underlying random access sequence, and imposes a subset
 * on the sequence based on a stride.  It can also be used to give standard container semantics to a
 * sequence, which is particularly useful with built-in types like double* that have iterator
 * semantics.
 *
 * \arg \a Ran A random access iterator type, such as double* or std::vector<double>::iterator.
 */
//================================================================================================//

template <typename Ran> class Slice {
private:
  using traits = typename std::iterator_traits<Ran>;

public:
  // NESTED CLASSES AND TYPEDEFS

  using value_type = typename traits::value_type;
  // no allocator type
  using size_type = std::size_t;
  using difference_type = typename traits::difference_type;

  class iterator {
  public:
    using iterator_category = typename traits::iterator_category;
    using value_type = typename traits::value_type;
    using difference_type = typename traits::difference_type;
    using pointer = typename traits::pointer;
    using reference = typename traits::reference;

    Ran first() const { return first_; }
    difference_type offset() const { return offset_; }
    size_type stride() const { return stride_; }

    iterator &operator++() {
      offset_ += stride_;
      return *this;
    }

    iterator operator++(int) {
      offset_ += stride_;
      return iterator(first_, offset_ - stride_, stride_);
    }

    reference operator*() const { return first_[offset_]; }

    reference operator[](difference_type i) const { return first_[offset_ + i]; }

    iterator operator+(difference_type i) const {
      return iterator(first_, offset_ + i * stride_, stride_);
    }

    difference_type operator-(iterator i) const {
      Require((first_ - i.first_) % stride_ == 0);
      Require(stride_ == i.stride_);

      return ((first_ - i.first_) + offset_ - i.offset_) / stride_;
    }

  private:
    friend class Slice;

    iterator(Ran const first, difference_type const offset, size_type const stride)
        : first_(first), offset_(offset), stride_(stride) {
      Require(stride > 0);
    }

    friend bool operator<(iterator a, iterator b) {
      return a.first_ - b.first_ < b.offset_ - a.offset_;
    }

    friend bool operator!=(iterator a, iterator b) {
      return a.first_ - b.first_ != b.offset_ - a.offset_;
    }

    Ran first_;
    difference_type offset_;
    size_type stride_;
  };

  class const_iterator {
  public:
    using iterator_category = typename traits::iterator_category;
    using value_type = typename traits::value_type;
    using difference_type = typename traits::difference_type;
    using pointer = typename traits::pointer;
    using reference = typename traits::reference;

    Ran first() const { return first_; }
    difference_type offset() const { return offset_; }
    size_type stride() const { return stride_; }

    const_iterator &operator++() {
      offset_ += stride_;
      return *this;
    }

    const_iterator operator++(int) {
      offset_ += stride_;
      return const_iterator(first_, offset_ - stride_, stride_);
    }

    value_type const &operator*() const { return first_[offset_]; }

    value_type const &operator[](difference_type i) const { return first_[offset_ + i]; }

    const_iterator operator+(difference_type i) const {
      return const_iterator(first_, offset_ + i * stride_, stride_);
    }

    difference_type operator-(const_iterator i) const {
      Require((first_ - i.first_) % stride_ == 0);
      Require(stride_ == i.stride_);

      return ((first_ - i.first_) + offset_ - i.offset_) / stride_;
    }

    const_iterator(iterator const &i) // NOLINT [hicpp-explicit-conversions]
        : first_(i.first()), offset_(i.offset()), stride_(i.stride()) {}

  private:
    friend class Slice;

    const_iterator(Ran const first, difference_type const offset, size_type const stride)
        : first_(first), offset_(offset), stride_(stride) {
      Require(stride > 0);
    }

    friend bool operator<(const_iterator a, const_iterator b) {
      return a.first_ - b.first_ < b.offset_ - a.offset_;
    }

    friend bool operator!=(const_iterator a, const_iterator b) {
      return a.first_ - b.first_ != b.offset_ - a.offset_;
    }

    Ran first_;
    difference_type offset_;
    size_type stride_;
  };

  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  using pointer = typename traits::pointer;
  using reference = typename traits::reference;
  using const_pointer = const value_type *;
  using const_reference = const value_type &;

  // CREATORS

  //----------------------------------------------------------------------------------------------//
  /*!
   * \brief Construct from a sequence.
   *
   * The constructed Slice has reference semantics to the sequence.  That is, modifications to
   * elements of the Slice are actually modifications to the elements of the underlying sequence,
   * and the Slice becomes invalid if the underlying sequence becomes invalid.  For example, a
   * sequence based on a vector becomes invalid if the vector is resized, and so will a Slice based
   * on that sequence.
   *
   * A stride is applied to all indexing, so that an index \c i applied to the Slice references the
   * element whose index is \c stride*i in the underlying sequence.  No offset is necessary, since
   * the same result is easily obtained by modifying the starting iterator of the underlying
   * sequence. Iterators behave in a comparable manner, so that \c begin() of the Slice points to \c
   * begin()+offset of the underlying container, and the increment operator actually increments by
   * \c stride.
   *
   * \param first_in Iterator to the beginning of a sequence.
   * \param length_in Length of the constructed Slice.
   * \param stride_in Stride to apply to the sequence.
   */
  Slice(Ran const first_in, size_type const length_in, size_type const stride_in = 1)
      : first(first_in), length(length_in), stride(stride_in) {
    Require(stride_in > 0);
  }

  // MANIPULATORS

  // ACCESSORS

  iterator begin() { return iterator(first, 0, stride); }
  const_iterator begin() const { return const_iterator(first, 0, stride); }
  iterator end() { return iterator(first, length * stride, stride); }
  const_iterator end() const { return const_iterator(first, length * stride, stride); }

  reverse_iterator rbegin();
  const_reverse_iterator rbegin() const;
  reverse_iterator rend();
  const_reverse_iterator rend() const;

  reference operator[](size_type n) {
    Require(n < size());
    return first[stride * n];
  }
  const_reference operator[](size_type n) const {
    Require(n < size());
    return first[stride * n];
  }

  reference at(size_type n);
  const_reference at(size_type n) const;

  reference front();
  const_reference front() const { return first[0]; }
  reference back();
  const_reference back() const { return first[stride * (size() - 1)]; }

  size_type size() const { return length; }
  bool empty() const { return size() == 0; }
  size_type max_size() const;

private:

  // DATA

  Ran first;
  size_type length;
  size_type stride;
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return a Slice.
 *
 * This is a convenient interface that avoids the necessity for writing a lot of type expressions in
 * user code.
 *
 * \tparam Ran A random access iterator.
 *
 * \param first Start of a sequence
 * \param length Length of the Slice.
 * \param stride Stride into the container
 *
 * \return The desired Slice.
 */
template <typename Ran>
inline Slice<Ran> slice(Ran const first, std::size_t const length, std::size_t const stride = 1) {
  Require(stride > 0);

  return Slice<Ran>(first, length, stride);
}

} // end namespace rtt_dsxx

#endif // container_Slice_hh

//------------------------------------------------------------------------------------------------//
// end of container/Slice.hh
//------------------------------------------------------------------------------------------------//
