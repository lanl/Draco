//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/Packing_Utils.hh
 * \author Thomas M. Evans, Tim Kelley <tkelley@lanl.gov>
 * \date   Thu Jul 19 11:27:46 2001
 * \brief  Packing Utilities, classes for packing stuff.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 *
 * This file contains classes and utilities that are used to "pack" data into byte-streams. The
 * byte-streams are represented by the char* type.  The following classes are:
 *
 * \arg \b Packer packing class
 * \arg \b Unpacker unpacking class
 */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_ds_Packing_Utils_hh
#define rtt_ds_Packing_Utils_hh

#include "Assert.hh"
#include "Endian.hh"
#include <map>
#include <vector>

namespace rtt_dsxx {

//================================================================================================//
/*!
 * \class Packer
 *
 * \brief Pack data into a byte stream.
 *
 * This class allows clients to \em "register" a \c char* stream and then load it with data of any
 * type.  This assumes that the \c sizeof(T) operator works and has meaning for the type.  Under the
 * hood it uses \c std::memcpy to perform the loading.  This class is easily understood by checking
 * the examples.
 *
 * No memory allocation is performed by the Packer.  However, the memory requirements may be
 * computed by putting the Packer into \c compute_buffer_size_mode().
 *
 * The benefit of using the Packer class is that byte copies are isolated into this one section of
 * code, thus obviating the need for reinterpret_cast statements in client code.  In fact, this
 * functionality conforms exactly to the ANSI C++ standard for copying byte-streams of data
 * (sec. 3.9). Additionally, bounds checking is performed on all stream packing operations.  This
 * bounds checking is always on.
 *
 * This class returns real char * pointers through its query functions.  We do not use the STL
 * iterator notation, even though that is how the pointers are used, so as not to confuse the fact
 * that these char * streams are \e continuous \e data byte-streams.  The pointers that are used to
 * "iterate" through the streams are real pointers, not an abstract iterator class.  So one could
 * think of these as iterators (they act like iterators) but they are real pointers into a
 * contiguous memory \c char* stream.
 *
 * Data can be unpacked using the Unpacker class.
 *
 * \example ds++/test/tstPacking_Utils.cc
 * Test the Packer and Unpacker classes.
 */
//================================================================================================//

class Packer {
public:
  // Typedefs.
  using pointer = char *;
  using const_pointer = const char *;

private:
  //! Size of packed stream.
  uint64_t stream_size{0};

  //! Pointer (mutable) into data stream.
  pointer ptr{nullptr};

  //! Pointers to begin and end of buffers.
  pointer begin_ptr{nullptr};
  pointer end_ptr{nullptr};

  //! If true, compute the stream_size required and do no packing.
  bool size_mode{false};

public:
  //! Constructor.
  Packer() = default;

  // Sets the buffer and puts the packer into pack mode.
  inline void set_buffer(uint64_t size_in, pointer buffer);

  //! Put the packer into compute buffer size mode.
  void compute_buffer_size_mode() {
    stream_size = 0;
    size_mode = true;
  }

  /*! In pack mode, pack values into the buffer.  In size mode, adds the size of the type into the
   *  total buffer size required. */
  template <typename T> inline void pack(const T &value);

  //! Accept data from another character stream.
  template <typename IT> void accept(uint64_t bytes, IT data);

  //! Advance the pointer without adding data. Useful for byte-aligning.
  inline void pad(uint64_t bytes);

  // >>> ACCESSORS

  //! Get a pointer to the current position of the data stream.
  const_pointer get_ptr() const {
    Require(!size_mode);
    return ptr;
  }

  //! Get a pointer to the beginning position of the data stream.
  const_pointer begin() const {
    Require(!size_mode);
    return begin_ptr;
  }

  //! Get a pointer to the ending position of the data stream.
  const_pointer end() const {
    Require(!size_mode);
    return end_ptr;
  }

  //! Get the size of the data stream.
  uint64_t size() const { return stream_size; }
};

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Set an allocated buffer to write data into.
 *
 * If \c compute_buffer_size_mode() is on, this function turns it off.
 *
 * This function accepts an allocated \c char* buffer.  It assigns begin and end pointers and a
 * mutable position pointer that acts like an iterator.  The Packer will write POD (Plain Old Data)
 * data into this buffer starting at the beginning address of the buffer.  This function must be
 * called before any \c Packer::pack calls can be made.
 *
 * Once \c Packer::set_buffer is called, all subsequent calls to \c Packer::pack will write data
 * incrementally into the buffer set by set_buffer.  To write data into a different buffer, call \c
 * Packer::set_buffer again; at this point the Packer no longer has any knowledge about the old
 * buffer.
 *
 * Note, the buffer must be allocated large enough to hold all the data that the client intends to
 * load into it.  There is no memory allocation performed by the Packer class; thus, the buffer
 * cannot be increased in size if a value is written past the end of the buffer.  Optionally, the
 * required buffer size may also be computed using the \c compute_buffer_size_mode().  See the \c
 * Packer::pack function for more details.
 *
 * \param size_in size of the buffer
 * \param buffer pointer to the char * buffer
 */
void Packer::set_buffer(uint64_t size_in, pointer buffer) {
  Require(buffer);

  size_mode = false;

  // set the size, begin and end pointers, and iterator
  stream_size = size_in;
  ptr = &buffer[0];
  begin_ptr = &buffer[0];
  end_ptr = begin_ptr + stream_size;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Depending on mode, pack data into a buffer, or compute increment to buffer size.
 *
 * This function's behavior depends on whether in compute_buffer_size_mode(), or not.
 *
 * In compute_buffer_size_mode(), the sizeof(T) operator is used to add the size of the data to the
 * total stream size.  Once this function is called for all of the data to be packed, the size()
 * member function may be used to retrieve the buffer size required.
 *
 * Note that using compute_buffer_size_mode() is optional.  See examples below.
 *
 * Regardless, once the user allocates the buffer, set_buffer() may then be called, which turns off
 * compute_buffer_size_mode (if on).  A call to pack() then actually packs its argument into the
 * buffer.  It also advances the pointer (iterator) location to the next location automatically.  It
 * uses the sizeof(T) operator to get the size of the data; thus, only data where sizeof() has
 * meaning will be properly written to the buffer.
 *
 * Packer::pack() does bounds checking to ensure that the buffer and buffer size defined by
 * Packer::set_buffer are consistent.  This bounds-checking is always on as the Packer is not
 * normally used in compute-intensive calculations.
 *
 * \param value data of type T to pack into the buffer; the data size must be accessible using the
 *              sizeof() operator.
 *
 * Example using compute_buffer_size_mode():
 \code
 double d1 = 5.0, d2 = 10.3;         // data to be packed
 Packer p;
 p.compute_buffer_size_mode();
 p << d1 << d2;                      // computes required size
 vector<char> buffer(p.size());      // allocate buffer
 p.set_buffer(p.size(), &buffer[0]);
 p << d1 << d2;                      // packs d1 and d2 into buffer
 \endcode
 *
 * Example not using compute_buffer_size_mode():
 \code
 double d1 = 5.0, d2 = 10.3;
 Packer p;
 uint64_t bsize = 2 * sizeof(double);  // compute buffer size
 vector<char> buffer(bsize);
 p.set_buffer(bsize, &buffer[0]);
 p << d1 << d2;                            // packs d1 and d2 into buffer
 \endcode
*/
template <typename T> void Packer::pack(T const &value) {
  if (size_mode)
    stream_size += sizeof(T);
  else {
    Require(begin_ptr);
    Ensure(ptr >= begin_ptr);
    Ensure(ptr + sizeof(T) <= end_ptr);

    // copy value into the buffer
    std::memcpy(ptr, &value, sizeof(T));

    // advance the iterator pointer to the next location
    ptr += sizeof(T);
  }
  return;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Add data from another character stream of a given size.
 *
 * \param bytes Number of bytes of data to copy.
 * \param data The data.
 */
template <typename IT> void Packer::accept(uint64_t bytes, IT data) {

  if (size_mode)
    stream_size += bytes;
  else {
    Require(begin_ptr);
    Require(ptr >= begin_ptr);
    Require(ptr + bytes <= end_ptr);

    while (bytes-- > 0)
      *(ptr++) = *(data++);
  }
  return;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Add the given number of blank bytes to the stream.
 *
 * \param bytes Number of bytes of blank data to add.
 */
void Packer::pad(uint64_t bytes) {
  while (bytes-- > 0)
    pack(char(0));
  return;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Stream out (<<) operator for packing data.
 *
 * The overloaded stream out operator can be used to pack data into streams (Packer p;
 * p.set_buffer(i,b); p << data;).  It simply calls the Packer::pack function.  It returns a
 * reference to the Packer object so that stream out operations can be strung together.
 *
 * This function also works when compute_buffer_size_mode() is on, in which case the total required
 * stream size is incremented.
 */
template <typename T> inline Packer &operator<<(Packer &p, const T &value) {
  // pack the value
  p.pack(value);

  // return the packer object
  return p;
}

//================================================================================================//
/*!
 * \class Unpacker
 *
 * \brief Unpack data from a byte stream.
 *
 * This class allows clients to "register" a char* stream and then unload data from it.  This
 * assumes that the sizeof(T) operator works and has meaning for the type.  Under the hood it uses
 * std::memcpy to perform the unloading.  This class is easily understood by checking the examples.
 *
 * No memory allocation is performed by the Unpacker.
 *
 * The benefit of using the Unpacker class is that byte copies are isolated into this one section of
 * code, thus obviating the need for reinterpret_cast statements in client code.  In fact, this
 * functionality conforms exactly to the ANSI C++ standard for copying byte-streams of data
 * (sec. 3.9). Additionally, bounds checking is performed on all stream packing operations.  This
 * bounds checking is always on.
 *
 * This class returns real char * pointers through its query functions.  We do not use the STL
 * iterator notation, even though that is how the pointers are used, so as not to confuse the fact
 * that these char * streams are \e continuous \e data byte-streams.  The pointers that are used to
 * "iterate" through the streams are real pointers, not an abstract iterator class.  So one could
 * think of these as iterators (they act like iterators) but they are real pointers into a
 * contiguous memory char * stream.
 *
 * This class is the complement to the Packer class.
 */
//================================================================================================//

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-private-field"
#endif

class Unpacker {
public:
  // Typedefs.
  using pointer = char *;
  using const_pointer = const char *;

private:
  //! Size of packed stream.
  uint64_t stream_size{0};

  //! Pointer (mutable) into data stream.
  const_pointer ptr{nullptr};

  //! Pointers to begin and end of buffers.
  const_pointer begin_ptr{nullptr};
  const_pointer end_ptr{nullptr};

  //! Should we convert the endian nature of the data?
  bool do_byte_swap{false};

public:
  //! Constructor.
  explicit Unpacker(bool byte_swap = false) : do_byte_swap(byte_swap) { /*...*/
  }

  //! Set the buffer.
  inline void set_buffer(uint64_t size_in, const_pointer buffer);

  //! Unpack value from buffer.
  template <typename T> inline void unpack(T &value);

  // >>> ACCESSORS

  //! Skip a specific number of bytes
  inline void skip(uint64_t bytes);

  //! Copy data to a provided iterator
  template <typename IT> void extract(uint64_t bytes, IT destination);

  //! Get a pointer to the current position of the data stream.
  const_pointer get_ptr() const { return ptr; }

  //! Get a pointer to the beginning position of the data stream.
  const_pointer begin() const { return begin_ptr; }

  //! Get a pointer to the ending position of the data stream.
  const_pointer end() const { return end_ptr; }

  //! Get the size of the data stream.
  uint64_t size() const { return stream_size; }
};

#if defined(__clang__)
#pragma clang diagnostic pop
#endif

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Set an allocated buffer to read data from.
 *
 * This function accepts an allocated char* buffer.  It assigns begin and end pointers and a mutable
 * position pointer that acts like an iterator.  The Unpacker will read POD data from this buffer
 * starting at the beginning address of the buffer.  This function must be called before any
 * Unpacker::unpack calls can be made.
 *
 * Once Unpacker::set_buffer is called, all subsequent calls to Unpacker::unpack will read data
 * incrementally from the buffer set by set_buffer.  To read data from a different buffer, call
 * Unpacker::set_buffer again; at this point the Unpacker no longer has any knowledge about the old
 * buffer.
 *
 * Note, there is no memory allocation performed by the Unpacker class.  Also, the client must know
 * how much data to read from the stream (of course checks can be made telling where the end of the
 * stream is located using the Unpacker::get_ptr, Unpacker::begin, and Unpacker::end functions).
 *
 * \param size_in size of the buffer
 * \param buffer const_pointer to the char * buffer
 */
void Unpacker::set_buffer(uint64_t size_in, const_pointer buffer) {
  Require(buffer);

  // set the size, begin and end pointers, and iterator
  stream_size = size_in;
  ptr = &buffer[0];
  begin_ptr = &buffer[0];
  end_ptr = begin_ptr + stream_size;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Unpack data from the buffer.
 *
 * This function unpacks a piece of data (single datum) from the buffer set by Unpacker::set_buffer.
 * It advances the pointer (iterator) location to the next location automatically.  It uses the
 * sizeof(T) operator to get the size of the data; thus, only data where sizeof() has meaning will
 * be properly read from the buffer.POLYNOMIAL_Specific_Heat_ANALYTIC_EoS_MODEL
 *
 * Unpacker::unpack() does bounds checking to ensure that the buffer and buffer size defined by
 * Unpacker::set_buffer are consistent.  This bounds-checking is always on as this should not be
 * used in computation intensive parts of the code.
 *
 * \param value data of type T to unpack from the buffer; the data size must be accessible using the
 *           sizeof() operator
 */
template <typename T> void Unpacker::unpack(T &value) {
  Require(begin_ptr);
  Ensure(ptr >= begin_ptr);
  Ensure(ptr + sizeof(T) <= end_ptr);

  std::memcpy(&value, ptr, sizeof(T));

  if (do_byte_swap)
    byte_swap(value);

  ptr += sizeof(T);
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Skip a specified number of bytes forward in the data stream
 *
 * \param bytes The number of bytes to skip.
 *
 * This is useful for data streams which have space inserted for alignment purposes.
 */
void Unpacker::skip(uint64_t bytes) {
  Require(begin_ptr);
  Check(ptr >= begin_ptr);
  Check(ptr + bytes <= end_ptr);

  ptr += bytes;
}

//------------------------------------------------------------------------------------------------//
/**
 * \brief Copy a piece of the data to memory referenced by the provided iterator.
 *
 * \param bytes The number of bytes to copy.
 * \param it The destination iterator. Must model ForwardIterator
 */
template <typename T> void Unpacker::extract(uint64_t bytes, T it) {
  Require(begin_ptr);
  Check(ptr >= begin_ptr);
  Check(ptr + bytes <= end_ptr);

  while (bytes-- > 0)
    *(it++) = *(ptr++);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Stream in (>>) operator for unpacking data.
 *
 * The overloaded stream in operator can be used to unpack data from streams (Unpacker u;
 * u.set_buffer(i,b); u >> data;).  It simply calls the Unpacker::unpack function.  It returns a
 * reference to the Unpacker object so that stream in operations can be strung together.
 */
template <typename T> inline Unpacker &operator>>(Unpacker &u, T &value) {
  u.unpack(value);
  return u;
}

//================================================================================================//
// PACKING/UNPACKING SHORTCUT FUNCTIONS
//================================================================================================//
/*!
 * \brief Packing function.
 *
 * This function uses the rtt_dsxx::Packer to pack a given field into a \c vector<char>.  The field
 * type is represented by the template argument \b FT.  The field type must have the following
 * members defined:
 *
 * \arg FT::value_type type stored in the field
 * \arg FT::const_iterator const iterator for the field
 * \arg FT::size() returns the number of elements in the field
 * \arg FT::begin() returns an iterator to the beginning of the field
 * \arg FT::end() returns an iterator to the end of the field
 * \arg FT::empty() determines if a container is empty
 *
 * Given these constraints, the function cannot be used to pack up a pointer array; however, this is
 * accomplished easily enough with the Packer class alone.
 *
 * The data in the field is packed into a \c vector<char>. The \c vector<char> passed to the
 * function must be empty; an assertion is thrown if it is not.  We do this for usage protection; we
 * want the user to be aware that data in the \c vector<char> would be destroyed.
 *
 * In summary, this is a simple function that is a shortcut for using the Packer class for fields
 * (\c vector<double>, \c list<int>, \c string, etc).  The complement of this function is
 * unpack_data, which takes a packed \c vector<char> and writes data into the field.
 *
 * The resulting size of the data stored in the \c vector<char> argument is
 * \code
 *   sizeof(FT::value_type) * field.size() + sizeof(int).
 * \endcode
 * The extra integer is for storing the size of the data field.
 *
 * \sa rtt_dsxx::Packer, tstPacking_Utils.cc, and rtt_dsxx::unpack_data
 *
 * \param field container or string
 * \param packed vector<char> that is empty; data will be packed into it
 */
template <typename FT> void pack_data(FT const &field, std::vector<char> &packed) {
  Require(packed.empty());

  // determine the size of the field
  Check(field.size() < INT_MAX);
  auto const field_size = static_cast<int>(field.size());

  // determine the number of bytes in the field
  Check(field_size * sizeof(typename FT::value_type) + sizeof(int) < INT32_MAX);
  auto const size = static_cast<int>(field_size * sizeof(typename FT::value_type) + sizeof(int));

  // make a vector<char> large enough to hold the packed field
  packed.resize(size);

  // make an unpacker and set it
  Packer packer;
  packer.set_buffer(size, &packed[0]);

  // pack up the number of elements in the field
  packer << field_size;

  // iterate and pack
  for (auto itr = field.begin(); itr != field.end(); itr++) {
    packer << *itr;
  }

  Ensure(packer.get_ptr() == &packed[0] + size);
  return;
}

//------------------------------------------------------------------------------------------------//
template <typename keyT, typename dataT>
void pack_data(std::map<keyT, dataT> const &map, std::vector<char> &packed) {
  Require(packed.empty());

  // determine the size of the field
  size_t const numkeys = map.size();

  // determine the number of bytes in the field
  size_t const key_size = numkeys * sizeof(keyT);
  size_t const data_size = numkeys * sizeof(dataT);

  // make a vector<char> large enough to hold the packed field
  size_t const size(sizeof(size_t) + key_size + data_size);
  packed.resize(size);

  // make an unpacker and set it
  Packer packer;
  packer.set_buffer(size, &packed[0]);

  // pack up the number of elements in the field
  packer << numkeys;

  // iterate and pack
  for (auto itr = map.begin(); itr != map.end(); itr++) {
    packer << (*itr).first;
    packer << (*itr).second;
  }

  Ensure(packer.get_ptr() == &packed[0] + size);
  return;
}

//------------------------------------------------------------------------------------------------//
template <typename keyT, typename dataT>
void pack_data(std::map<keyT, std::vector<dataT>> const &map, std::vector<char> &packed) {
  Require(packed.empty());

  // determine the size of the field
  size_t const numkeys = map.size();

  // determine the number of bytes in the field
  size_t const key_size = numkeys * sizeof(keyT) + sizeof(size_t);
  size_t data_size(0);
  for (auto itr = map.begin(); itr != map.end(); ++itr) {
    // Size of data plus a size_t that indicates the vector length.
    data_size += (*itr).second.size() * sizeof(dataT) + sizeof(size_t);
  }

  // make a vector<char> large enough to hold the packed field
  size_t const size(key_size + data_size);
  packed.resize(size);

  // make an unpacker and set it
  Packer packer;
  packer.set_buffer(size, &packed[0]);

  // pack up the number of keys in the map.
  packer << numkeys;

  // iterate and pack:
  // 1. the keys
  for (auto itr = map.begin(); itr != map.end(); itr++) {
    packer << (*itr).first;
  }
  // 2. The vector size and vector data for each key.
  for (auto itr = map.begin(); itr != map.end(); itr++) {
    // The size of the fector associated with one key.
    packer << (*itr).second.size();
    // pack the data found in the vector for one key.
    for (auto it = (*itr).second.begin(); it != (*itr).second.end(); it++) {
      packer << *it;
    }
  }

  Ensure(packer.get_ptr() == &packed[0] + size);
  return;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Unpacking function.
 *
 * This function uses the rtt_dsxx::Unpacker to unpack a given field from a vector<char>.  The field
 * type is represented by the template argument FT.  The field type must have the following members
 * defined:
 *
 * \arg FT::iterator const iterator for the field
 * \arg FT::resize() returns the number of elements in the field
 * \arg FT::begin() returns an iterator to the beginning of the field
 * \arg FT::end() returns an iterator to the end of the field
 * \arg FT::empty() determines if a container is empty
 *
 * Given these constraints, the function cannot be used to unpack a pointer array; however, this is
 * accomplished easily enough with the Unpacker class alone.
 *
 * The data in the field is unpacked from a vector<char>. The data in the vector<char> must be
 * packed in a manner consistent with pack_data.  The function checks for this.  Additionally, the
 * field given to the function must be empty.
 *
 * In summary, this is a simple function that is a shortcut for using the Unpacker class for fields
 * (vector<double>, list<int>, string, etc).  The complement of this function is pack_data, which
 * packs fields into a vector<char>
 *
 * The correct size of the vector<char> containing the data is:
 * \code
 *   sizeof(FT::value_type)*field_size + sizeof(int)
 * \endcode
  * where field_size is the size of the resulting field. So you'd better know this in advance.
 *
 * \sa rtt_dsxx::Unpacker, tstPacking_Utils.cc, and rtt_dsxx::pack_data
 *
 * \param field container or string that is empty; data will be unpacked into it
 * \param packed vector<char> created by pack_data function (or in a manner
 *           analogous)
 */
template <typename FT> void unpack_data(FT &field, std::vector<char> const &packed) {
  Require(field.empty());
  Require(packed.size() >= sizeof(int));

  // make an unpacker and set it
  Unpacker unpacker;
  unpacker.set_buffer(packed.size(), &packed[0]);

  // unpack the number of elements in the field
  int field_size = 0;
  unpacker >> field_size;

  // make a field big enough to hold all the elements
  field.resize(field_size);

  // unpack the data
  for (auto itr = field.begin(); itr != field.end(); itr++)
    unpacker >> *itr;

  Ensure(unpacker.get_ptr() == &packed[0] + packed.size());
  return;
}

//------------------------------------------------------------------------------------------------//
template <typename keyT, typename dataT>
void unpack_data(std::map<keyT, dataT> &unpacked_map, std::vector<char> const &packed) {
  Require(unpacked_map.empty());
  Require(packed.size() >= sizeof(int));

  // make an unpacker and set it
  Unpacker unpacker;
  unpacker.set_buffer(packed.size(), &packed[0]);

  // unpack the number of elements in the field
  size_t numkeys(0);
  unpacker >> numkeys;

  // unpack the keys
  keyT key;
  dataT data;
  for (size_t i = 0; i < numkeys; ++i) {
    unpacker >> key;
    unpacker >> data;
    unpacked_map[key] = data;
  }

  Ensure(unpacker.get_ptr() == &packed[0] + packed.size());
  return;
}

//------------------------------------------------------------------------------------------------//
template <typename keyT, typename dataT>
void unpack_data(std::map<keyT, std::vector<dataT>> &unpacked_map,
                 std::vector<char> const &packed) {
  Require(unpacked_map.empty());
  Require(packed.size() >= sizeof(size_t));

  // make an unpacker and set it
  Unpacker unpacker;
  unpacker.set_buffer(packed.size(), &packed[0]);

  // unpack the number of elements in the field
  size_t numkeys(0);
  unpacker >> numkeys;

  // unpack the keys
  keyT key(0);
  for (size_t i = 0; i < numkeys; ++i) {
    unpacker >> key;
    unpacked_map[key] = std::vector<dataT>();
  }

  // unpack the data
  for (auto it = unpacked_map.begin(); it != unpacked_map.end(); ++it) {
    size_t numdata(0);
    unpacker >> numdata;
    unpacked_map[(*it).first].resize(numdata);
    for (auto itr = (*it).second.begin(); itr != (*it).second.end(); itr++)
      unpacker >> *itr;
  }

  Ensure(unpacker.get_ptr() == &packed[0] + packed.size());
  return;
}

//------------------------------------------------------------------------------------------------//
// GLOBAL scope functions
//------------------------------------------------------------------------------------------------//

/*!
 * \brief Pack an array into a char buffer while honoring endianess.
 * \param[in]  start iterator to the start of the buffer that will be packed
 * \param[out] dest iterator to the start of the buffer holding the packed data.
 * \param[in]  num_elements The number of elements to be packed starting with start
 * \param[in]  byte_swap (default: false)
 */
inline void pack_vec_double(double const *start, char *dest, uint32_t num_elements,
                            bool byte_swap = false) {
  rtt_dsxx::Packer packer;
  packer.set_buffer(num_elements * sizeof(double), dest);

  if (!byte_swap) {
    for (size_t i = 0; i < num_elements; ++i) {
      packer << *start++;
    }
  } else {
    for (size_t i = 0; i < num_elements; ++i) {
      packer << byte_swap_copy(*start++);
    }
  }
  return;
}

} // end namespace rtt_dsxx

#endif // rtt_ds_Packing_Utils_hh

//------------------------------------------------------------------------------------------------//
// end of ds++/Packing_Utils.hh
//------------------------------------------------------------------------------------------------//
