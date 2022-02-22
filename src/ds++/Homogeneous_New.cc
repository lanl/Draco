//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/Homogeneous_New.cc
 * \author Kent Budge
 * \brief  Implement methods of class Homogeneous_New
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Homogeneous_New.hh"
#include "Assert.hh"
#include <algorithm>
#include <cstdlib>

namespace rtt_dsxx {

//------------------------------------------------------------------------------------------------//
Homogeneous_New::Homogeneous_New(unsigned const object_size, unsigned const default_block_size)
    : object_size_(object_size), default_block_size_(default_block_size),
      total_number_of_segments_(0), first_block_(nullptr), first_segment_(nullptr) {

  Require(object_size >= sizeof(void *));
  Require(default_block_size > sizeof(void *));

  if ((default_block_size_ - sizeof(void *)) % object_size != 0) {
    unsigned count = (default_block_size_ - static_cast<unsigned>(sizeof(void *))) / object_size;
    ++count;
    default_block_size_ = count * object_size + static_cast<unsigned>(sizeof(void *));
  }

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
bool Homogeneous_New::check_class_invariants() const {
  // No checkable invariants at present
  return true;
}

//------------------------------------------------------------------------------------------------//
Homogeneous_New::~Homogeneous_New() {
  void *next = first_block_;
  while (next != nullptr) {
    void *next_next = *reinterpret_cast<void **>(next);
    delete[] reinterpret_cast<char *>(next);
    next = next_next;
  }
}

//------------------------------------------------------------------------------------------------//
void *Homogeneous_New::allocate() {
  if (first_segment_ != nullptr) {
    void *Result = first_segment_;
    first_segment_ = *reinterpret_cast<void **>(Result);

    Ensure(check_class_invariants());
    return Result;
  } else {
    // have to allocate a new block
    allocate_block_(default_block_size_);

    Ensure(check_class_invariants());
    return allocate();
  }
}

//------------------------------------------------------------------------------------------------//
void Homogeneous_New::allocate_block_(unsigned const /*block_size*/) {
  auto *const new_first_block = new char[default_block_size_];
  *reinterpret_cast<void **>(new_first_block) = first_block_;
  first_block_ = new_first_block;
  char *segptr = new_first_block + object_size_;
  char *new_segptr = segptr + object_size_;
  while (new_segptr < new_first_block + default_block_size_) {
    *reinterpret_cast<void **>(segptr) = first_segment_;
    first_segment_ = segptr;
    segptr = new_segptr;
    new_segptr += object_size_;
    ++total_number_of_segments_;
  }

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
void Homogeneous_New::deallocate(void *const ptr) {
  *reinterpret_cast<void **>(ptr) = first_segment_;
  first_segment_ = ptr;

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
void Homogeneous_New::reserve(unsigned const object_count) {
  if (object_count > total_number_of_segments_) {
    unsigned count = object_count - total_number_of_segments_;
    count = count * object_size_ + static_cast<unsigned>(sizeof(void *));
    count = std::max(count, default_block_size_);
    allocate_block_(count);
  }

  Ensure(check_class_invariants());
}

} // end namespace rtt_dsxx

//------------------------------------------------------------------------------------------------//
// end of Homogeneous_New.cc
//------------------------------------------------------------------------------------------------//
