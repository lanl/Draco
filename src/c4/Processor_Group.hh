//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/Processor_Group.hh
 * \brief  Definition of class Processor_Group
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef c4_Processor_Group_hh
#define c4_Processor_Group_hh

#include "c4/config.h"
#include <vector>

#ifdef C4_MPI
#include "c4_mpi.h"
#endif // C4_MPI

namespace rtt_c4 {

//================================================================================================//
/*!
 * \class Processor_Group
 * \brief Representation of subgroup of processors
 *
 * This class provides a parallel interface for a group of processors that is a subset of the entire
 * process set.  This allows communications libraries to do things like sums over process groups
 * efficiently.  In other words, this class is essentially a representation of an MPI communicator.
 */
//================================================================================================//

class Processor_Group {
public:
  // NESTED CLASSES AND TYPEDEFS

  // CREATORS

  //! Create a Process_Group based on a stride through the ranks.
  explicit Processor_Group(unsigned const stride);

  //! Destructor.
  ~Processor_Group();

  //! Copy/Move assignment/construction Not implemented
  Processor_Group(const Processor_Group &rhs) = delete;
  Processor_Group(const Processor_Group &&rhs) = delete;
  Processor_Group &operator=(const Processor_Group &rhs) = delete;
  Processor_Group &operator=(const Processor_Group &&rhs) = delete;

  // ACCESSORS

  //! Get the number of processors in the group.
  unsigned size() const { return size_; }

  bool check_class_invariants() const { return (size() > 0); }

  // SERVICES

  //! Sum a set of values over the group, returning the sum to all processors.
  template <typename RandomAccessContainer> void sum(RandomAccessContainer &values);

  /*!
   * \brief Assemble a set of local vectors into global vectors (container-based version).
   */
  template <typename T>
  void assemble_vector(std::vector<T> const &local_vector, std::vector<T> &global_vector) const;

  /*!
   * \brief Assemble a set of local vectors into global vectors (pointer-based version).
   */
  template <typename T>
  void assemble_vector(T const *local_vector, T *global_vector, unsigned const N) const;

private:
  // DATA

  unsigned size_;

#ifdef C4_MPI
  MPI_Group group_;
  MPI_Comm comm_;
#endif // C4_MPI
};

} // end namespace rtt_c4

#include "Processor_Group.i.hh"

#endif // c4_Processor_Group_hh

//------------------------------------------------------------------------------------------------//
// end of c4/Processor_Group.hh
//------------------------------------------------------------------------------------------------//
