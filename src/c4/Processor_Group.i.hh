//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Processor_Group.i.hh
 * \author Kent Budge
 * \date   Fri Oct 20 13:49:10 2006
 * \brief  Template method definitions of class Processor_Group
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef c4_Processor_Group_i_hh
#define c4_Processor_Group_i_hh

#include "MPI_Traits.hh"
#include "c4/config.h"
#include "ds++/Assert.hh"

#ifdef C4_MPI

namespace rtt_c4 {
//---------------------------------------------------------------------------//
template <typename T> void Processor_Group::sum(std::vector<T> &x) {
  // copy data into send buffer
  std::vector<T> y = x;

  // do global MPI reduction (result is on all processors) into x
  int status =
      MPI_Allreduce(&y[0], &x[0], y.size(),
                    rtt_c4::MPI_Traits<T>::element_type(), MPI_SUM, comm_);

  Insist(status == 0, "MPI_Allreduce failed");
}

//---------------------------------------------------------------------------//
template <typename T>
void Processor_Group::assemble_vector(std::vector<T> const &local,
                                      std::vector<T> &global) const {
  global.resize(local.size() * size());

  int status =
      MPI_Allgather(const_cast<T *>(&local[0]), local.size(),
                    rtt_c4::MPI_Traits<T>::element_type(), &global[0],
                    local.size(), rtt_c4::MPI_Traits<T>::element_type(), comm_);

  Insist(status == 0, "MPI_Gather failed");
}

//---------------------------------------------------------------------------//
/*!
 * \param local Points to a region of storage of size N.
 * \param global Points to a region of storage of size N*this->size()
 * \param N Number of local quantities to assemble.
 */

template <typename T>
void Processor_Group::assemble_vector(T const *local, T *global,
                                      unsigned const N) const {
  int status = MPI_Allgather(const_cast<T *>(local), N,
                             rtt_c4::MPI_Traits<T>::element_type(), global, N,
                             rtt_c4::MPI_Traits<T>::element_type(), comm_);

  Insist(status == 0, "MPI_Gather failed");
}

} // end namespace rtt_c4

#endif // C4_MPI

#endif // c4_Processor_Group_i_hh

//---------------------------------------------------------------------------//
// end of c4/Processor_Group.i.hh
//---------------------------------------------------------------------------//
