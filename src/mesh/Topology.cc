//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/Topology.cc
 * \author H.Park
 * \brief  Topology class implementation file.
 * \note   Copyright (C) 2023-2026 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "mesh/Topology.hh"
#include "dsxx/Assert.hh"
#include <cstdint> // int64_t
#include <iostream>
#include <string>
#include <utility> // move

namespace rtt_mesh {

//------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Topology abstract base class constructor.
 *
 * This constructor leaves the number of global cells data member unset; thus, it must be manually
 * set by the derived class.
 *
 * \param[in] scheme string stating the parallel scheme, must be "replication", "DD", or
 *        "DD/replication"
 */
Topology::Topology(std::string scheme)
    : num_global_cells_d(0), parallel_scheme_d(std::move(scheme)) {
  Ensure(parallel_scheme_d == "replication" || parallel_scheme_d == "DD" ||
         parallel_scheme_d == "DD/replication");
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Topology abstract base class constructor.
 *
 * \param[in] num_global_cells the number of global cells in the problem.
 * \param[in] scheme string stating the parallel scheme, must be "replication", "DD", or
 *        "DD/replication"
 */
Topology::Topology(const int64_t num_global_cells, std::string scheme)
    : num_global_cells_d(num_global_cells), parallel_scheme_d(std::move(scheme)) {
  Ensure(parallel_scheme_d == "replication" || parallel_scheme_d == "DD" ||
         parallel_scheme_d == "DD/replication");
  Ensure(num_global_cells_d > 0);
}

//------------------------------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//------------------------------------------------------------------------------------------------//

/*!
 * \brief Enable operator<< overload for Topology objects.
 * \param[out] out The object's state will be printed to this output stream.
 * \param[in] topo The object whose data state will be printed.
 */
std::ostream &operator<<(std::ostream &out, const Topology &topo) {
  topo.print(out);
  return out;
}

} // namespace rtt_mesh

//------------------------------------------------------------------------------------------------//
// end of Topology.cc
//------------------------------------------------------------------------------------------------//
