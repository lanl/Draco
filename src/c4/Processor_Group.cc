//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/Processor_Group.cc
 * \brief  Define non-template methods of class Processor_Group.
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/config.h"
#include <iostream>

#include "C4_Functions.hh"
#include "Processor_Group.hh"
#include "ds++/Assert.hh"

namespace rtt_c4 {

using namespace std;

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for Processor_Group
 * \param stride Number of processor groups to create. Each processor with rank \c pid will be
 *        assigned to processor group <code> pid mod stride</code>. For example, if the stride is 3,
 *        then processors 0, 3, 6, ... form the first processor group, processors 1, 4, 7, ...  form
 *        the second processor group, and processors 2, 5, 8, ... form the third processor group.
 */
//------------------------------------------------------------------------------------------------//

#ifdef C4_MPI

Processor_Group::Processor_Group(unsigned const stride) : size_(0), group_(), comm_() {
  int flag;
  MPI_Initialized(&flag);
  Insist(flag, "Processor_Group created before MPI is initialized");

  unsigned const number_of_processors = rtt_c4::nodes();
  unsigned const pid = rtt_c4::node();

  MPI_Group parent_group;
  flag = MPI_Comm_group(MPI_COMM_WORLD, &parent_group);

  vector<int> ranks;
  for (unsigned i = 0; i < number_of_processors; ++i) {
    if (i % stride == pid % stride) {
      ranks.push_back(i);
    }
  }
  Check(ranks.size() < UINT_MAX);
  size_ = static_cast<unsigned>(ranks.size());
  flag = MPI_Group_incl(parent_group, size_, &ranks[0], &group_);

  flag = MPI_Comm_create(MPI_COMM_WORLD, group_, &comm_);

  Ensure(check_class_invariants());
}

Processor_Group::~Processor_Group() {
  int flag;
  MPI_Finalized(&flag);
  if (!flag) {
    MPI_Comm_free(&comm_);
    MPI_Group_free(&group_);
  }
  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//

#else // not C4_MPI

Processor_Group::Processor_Group(unsigned const /*stride*/) : size_(1) {
  Check(rtt_c4::nodes() == 1);
  Check(rtt_c4::node() == 0);

  Ensure(check_class_invariants());
}

Processor_Group::~Processor_Group() { Ensure(check_class_invariants()); }

#endif //C4_MPI

//------------------------------------------------------------------------------------------------//

} // end namespace rtt_c4

//------------------------------------------------------------------------------------------------//
// end of Processor_Group.cc
//------------------------------------------------------------------------------------------------//
