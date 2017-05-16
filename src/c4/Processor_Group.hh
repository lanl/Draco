//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Processor_Group.hh
 * \author Kent Budge
 * \brief  Definition of class Processor_Group
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef c4_Processor_Group_hh
#define c4_Processor_Group_hh

#include "c4/config.h"

#ifdef C4_MPI

#include "c4_mpi.h"
#include <vector>

namespace rtt_c4 {

//===========================================================================//
/*!
 * \class Processor_Group
 * \brief Representation of subgroup of processors
 *
 * This class provides a parallel interface for a group of processors that is
 * a subset of the entire process set.  This allows communications libraries
 * to do things like sums over process groups efficiently.  In other words,
 * this class is essentially a representation of an MPI communicator.
 */
//===========================================================================//

class DLL_PUBLIC_c4 Processor_Group {
public:
  // NESTED CLASSES AND TYPEDEFS

  // CREATORS

  //! Create a Process_Group based on a stride through the ranks.
  explicit Processor_Group(unsigned const stride);

  //! Destructor.
  ~Processor_Group();

  // ACCESSORS

  //! Get the number of processors in the group.
  unsigned size() const { return size_; }

  bool check_class_invariants() const { return true; }

  // SERVICES

  //! Sum a set of values over the group, returning the sum to all
  //! processors.
  template <typename RandomAccessContainer>
  void sum(RandomAccessContainer &values);

  //! Assemble a set of local vectors into global vectors.
  template <typename T>
  void assemble_vector(std::vector<T> const &local_vector,
                       std::vector<T> &global_vector) const;

  //! Assemble a set of local vectors into global vectors.
  template <typename T>
  void assemble_vector(T const *local_vector, T *global_vector,
                       unsigned count) const;

private:
  // NESTED CLASSES AND TYPEDEFS

  // IMPLEMENTATION

  //! Not implemented
  Processor_Group(const Processor_Group &rhs);

  //! Not implemented
  Processor_Group &operator=(const Processor_Group &rhs);

  // DATA

  unsigned size_;
  MPI_Group group_;
  MPI_Comm comm_;
};

} // end namespace rtt_c4

#endif // C4_MPI

#endif // c4_Processor_Group_hh

//---------------------------------------------------------------------------//
// end of c4/Processor_Group.hh
//---------------------------------------------------------------------------//
