//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/Processor_Group.i.hh
 * \brief  Template method definitions of class Processor_Group
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef c4_Processor_Group_i_hh
#define c4_Processor_Group_i_hh

#ifdef C4_MPI
#include "MPI_Traits.hh"
#endif // C4_MPI

#include "Processor_Group.hh"
#include "ds++/Assert.hh"

namespace rtt_c4 {

#ifdef C4_MPI

//------------------------------------------------------------------------------------------------//
template <typename RandomAccessContainer> void Processor_Group::sum(RandomAccessContainer &x) {
  using T = typename RandomAccessContainer::value_type;

  // copy data into send buffer
  std::vector<T> y(x.begin(), x.end());

  // do global MPI reduction (result is on all processors) into x
  Check(y.size() < INT_MAX);
  int status = MPI_Allreduce(&y[0], &x[0], static_cast<int>(y.size()),
                             rtt_c4::MPI_Traits<T>::element_type(), MPI_SUM, comm_);

  Insist(status == 0, "MPI_Allreduce failed");
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Assemble a set of local vectors into global vectors (container-based version).
 *
 * \param[in]  local  Points to a region of storage of size N.
 * \param[out] global Points to a region of storage of size N*this->size()
 */
template <typename T>
void Processor_Group::assemble_vector(std::vector<T> const &local, std::vector<T> &global) const {
  global.resize(local.size() * size());

  Check(local.size() < INT_MAX);
  int status =
      MPI_Allgather(const_cast<T *>(&local[0]), static_cast<int>(local.size()),
                    rtt_c4::MPI_Traits<T>::element_type(), &global[0],
                    static_cast<int>(local.size()), rtt_c4::MPI_Traits<T>::element_type(), comm_);

  Insist(status == 0, "MPI_Gather failed");
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Assemble a set of local vectors into global vectors (pointer-based version).
 *
 * \param[in]  local  Points to a region of storage of size N.
 * \param[out] global Points to a region of storage of size N*this->size()
 * \param[in]  N      Number of local quantities to assemble.
 */
template <typename T>
void Processor_Group::assemble_vector(T const *local, T *global, unsigned const N) const {
  int status = MPI_Allgather(const_cast<T *>(local), N, rtt_c4::MPI_Traits<T>::element_type(),
                             global, N, rtt_c4::MPI_Traits<T>::element_type(), comm_);

  Insist(status == 0, "MPI_Gather failed");
}

#else // not C4_MPI

//------------------------------------------------------------------------------------------------//
template <typename RandomAccessContainer> void Processor_Group::sum(RandomAccessContainer & /*x*/) {
  // noop
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Assemble a set of local vectors into global vectors (container-based version).
 *
 * \param[in]  local  Points to a region of storage of size N.
 * \param[out] global Points to a region of storage of size N*this->size()
 */
template <typename T>
void Processor_Group::assemble_vector(std::vector<T> const &local, std::vector<T> &global) const {
  global = local;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Assemble a set of local vectors into global vectors (pointer-based version).
 *
 * \param[in]  local  Points to a region of storage of size N.
 * \param[out] global Points to a region of storage of size N*this->size()
 * \param[in]  N      Number of local quantities to assemble.
 */
template <typename T>
void Processor_Group::assemble_vector(T const *local, T *global, unsigned const N) const {
  std::copy(local, local + N, global);
}

#endif // C4_MPI
} // end namespace rtt_c4

#endif // c4_Processor_Group_i_hh

//------------------------------------------------------------------------------------------------//
// end of c4/Processor_Group.i.hh
//------------------------------------------------------------------------------------------------//
