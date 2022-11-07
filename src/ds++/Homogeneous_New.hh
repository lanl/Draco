//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/Homogeneous_New.hh
 * \author Kent Budge
 * \date   Tue Nov 28 08:27:37 2006
 * \brief  Definition of class Homogeneous_New
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef dsxx_Homogeneous_New_hh
#define dsxx_Homogeneous_New_hh

#include "ds++/config.h"

namespace rtt_dsxx {

//================================================================================================//
/*!
 * \class Homogeneous_New
 * \brief Allocator for large number of individual objects of identical size
 *
 * This class optimizes allocation of memory under the following circumstances:
 *
 * - Memory is allocated or deallocated for one object at a time.
 * - The objects are all identical in size.
 * - Large numbers of these objects will be allocated and deallocated in an inpredictable manner
 * .
 *
 * This class does not support memory requests of any size other than the object size specified at
 * construction. However, requests of the correct size are handled more efficiently than is possible
 * for the C++ library default allocator, or for any other allocator that handles requests of
 * arbitrary size.
 *
 * Note that this class is \em not suitable for implementing STL allocator classes, which is why we
 * call it \c Homogeneous_New rather than \c Homogeneous_Allocator.  STL allocators generally must
 * be able to support memory requests of varying size and have other requirements that cannot be
 * supported well by this class.
 *
 * The implementation of this class is a linked list of blocks of memory whose size should be tuned
 * to minimize memory fragmentation without wasting excessive memory in the last, partially filled,
 * block in the list.  Each block is itself divided into segments of the anticipated object size,
 * and all unused segments from all the memory blocks are organized into their own linked list.
 *
 * When a memory request is received, the address of the first free segment is returned and the next
 * free segment is moved to the top of the list.  When a memory release is received, the returned
 * memory segment is placed at the top of the list.  These are fast operations unless all segments
 * are in use and a new memory block must be allocated.  Memory segments currently in use are not
 * part of the segment list and need no extraneous storage for management.  The only extraneous
 * storage is that required for the linked list of contiguous memory blocks in which the segments
 * live.  Even this extraneous storage could be eliminated if we did not care about destruction of
 * the allocator itself.
 *
 * Unlike the memory segments, the memory blocks from which the segments are drawn need not be
 * identical in size (though this is the default behavior.)  A user may at any time call the \c
 * reserve function to indicate that a number of objects are about to be allocated.  If the required
 * number of segments are not available to fill the anticipated requests, a new memory block will be
 * allocated that is at least large enough to accomodate all the anticipated requests.  If the user
 * anticipates that the maximum number of allocations will never exceed \c N, he should call <code>
 * reserve(N) </code> immediately after constructing the allocator.
 *
 * Example: A finite element code allocates cell objects one at a time off the heap, because each
 * cell type is represented by its own class (descended from a common type.)  There will likely be
 * large numbers of each cell type. The initial number of cells of each type is known, and the user
 * anticipates that mesh refinement and load balancing will not change these numbers by more than
 * about 20%. \c Homogeneous_New should then be ideal for implementing an <code> operator new
 * </code> for each cell type.  The allocators should initially \c reserve 1.2 times the number of
 * cells of each type.
 *
 * The concepts underlying this class are from a paper by Andy Koenig presented at a USENIX meeting
 * in the early 1990s.
 */
//================================================================================================//

class Homogeneous_New {
public:
  // NESTED CLASSES AND TYPEDEFS

  enum { DEFAULT_BLOCK_SIZE = 4096 };

  // CREATORS

  //! Create an allocator for objects of the specified size, using the specified block byte count
  //! (defaulting to a system-tuned value.)
  explicit Homogeneous_New(unsigned object_size, unsigned default_block_size = DEFAULT_BLOCK_SIZE);

  //! Copy constructor: not implemented
  Homogeneous_New(const Homogeneous_New &rhs) = delete;
  //! Move constructor: not implemented
  Homogeneous_New(const Homogeneous_New &&rhs) = delete;
  //! Assignment copy operator for Homogeneous_New: not implemented.
  Homogeneous_New &operator=(const Homogeneous_New &rhs) = delete;
  //! Assignment move operator for Homogeneous_New: not implemented.
  Homogeneous_New &operator=(const Homogeneous_New &&rhs) = delete;

  //! Destructor.
  ~Homogeneous_New();

  // MANIPULATORS

  // Allocate storage for a single object.
  void *allocate();

  // Release an object's storage.
  void deallocate(void *ptr);

  // Reserve storage for the specified number of objects.
  void reserve(unsigned object_count);

  // ACCESSORS

  bool check_class_invariants() const;

private:
  // NESTED CLASSES AND TYPEDEFS

  // IMPLEMENTATION

  void allocate_block_(unsigned const block_size);

  // DATA

  unsigned object_size_;
  unsigned default_block_size_;
  unsigned total_number_of_segments_;
  void *first_block_;
  void *first_segment_;
};

} // end namespace rtt_dsxx

#endif // dsxx_Homogeneous_New_hh

//------------------------------------------------------------------------------------------------//
// end of ds++/Homogeneous_New.hh
//------------------------------------------------------------------------------------------------//
