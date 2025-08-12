//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/Draco_Mesh_Splitter.hh
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \date   Thursday, Jul 31, 2025, 15:00 pm
 * \brief  Draco_Mesh class header file.
 * \note   Copyright (C) 2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_mesh_Draco_Mesh_Splitter_hh
#define rtt_mesh_Draco_Mesh_Splitter_hh

#include "dsxx/Assert.hh"
#include "c4/C4_Functions.hh"
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace rtt_mesh {

// forward declaration
class Draco_Mesh;

//================================================================================================//
/*!
 * \class Draco_Mesh_Splitter
 *
 * \brief Decompose a global Draco_Mesh to rank-local Draco_Mesh objects, provided a partition
 *
 * The Draco_Mesh_Splitter class ...
 */
//================================================================================================//

class Draco_Mesh_Splitter {
private:
  // Partition of global cells over ranks
  std::vector<std::vector<int64_t>> global_cells_per_rank;

public:
  //! Constructor.
  explicit Draco_Mesh_Splitter(const std::vector<std::vector<int64_t>> &global_cells_per_rank_);

  //! Split the global mesh, producing a rank-local mesh (the main service of this class)
  std::shared_ptr<Draco_Mesh> split_mesh(const std::shared_ptr<Draco_Mesh> &global_mesh);

  //! Split an arbitrary global vector into local vectors according to global_cells_per_rank
  template <typename T> std::vector<T> split_vector(const std::vector<T> &global_vector);

  //! Reset the global cell partition over ranks
  void reset_global_cells_per_rank(const std::vector<std::vector<int64_t>> &gcells_rank_new) {
    global_cells_per_rank = gcells_rank_new;
  }
};

//------------------------------------------------------------------------------------------------//
// IMPLEMENTATION
//------------------------------------------------------------------------------------------------//
/*!
 * \brief split a global vector into rank-local vectors (assuming the global_vector is on each rank)
 *
 * \param[in] global_vector the global vector to partition
 *
 * \return an MPI rank-local vector for the partition of cells on the rank
 */
template <typename T>
std::vector<T> Draco_Mesh_Splitter::split_vector(const std::vector<T> &global_vector) {
  Ensure(global_vector.size() > 0);

  // get rank-local data from global cell partition
  const int my_rank = rtt_c4::node();
  const auto &my_gcells = global_cells_per_rank[my_rank];
  const size_t num_lcells = my_gcells.size();

  // check global vector is correct size for partition
  Remember(auto g_it = std::max_element(my_gcells.begin(), my_gcells.end()));
  Remember(auto gcell_max = *g_it);
  Remember(rtt_c4::global_max(gcell_max));
  Check(gcell_max == static_cast<int64_t>(global_vector.size() - 1));

  // intialized and populate local return vector
  std::vector<T> local_vector(my_gcells.size());
  for (size_t lcell = 0; lcell < num_lcells; ++lcell) {
    const auto &gcell = my_gcells[lcell];
    local_vector[lcell] = global_vector[gcell];
  }

  // return rank local vector
  return local_vector;
}

} // end namespace rtt_mesh

#endif // rtt_mesh_Draco_Mesh_Splitter_hh

//------------------------------------------------------------------------------------------------//
// end of mesh/Draco_Mesh_Splitter.hh
//------------------------------------------------------------------------------------------------//
