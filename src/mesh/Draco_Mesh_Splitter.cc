//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/Draco_Mesh_Splitter.cc
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \date   Thursday, Jul 31, 2025, 16:13 pm
 * \brief  Draco_Mesh class header file.
 * \note   Copyright (C) 2025 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Draco_Mesh_Splitter.hh"
#include "Draco_Mesh.hh"
#include "dsxx/Assert.hh"
#include "c4/C4_Functions.hh"
#include "mesh_element/Geometry.hh"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <set>
#include <vector>

namespace rtt_mesh {

//------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Draco_Mesh_Splitter constructor.
 *
 * \param[in] global_cells_per_rank_ vector of vector of global cell indexes per MPI rank
 */
Draco_Mesh_Splitter::Draco_Mesh_Splitter(
    const std::vector<std::vector<int64_t>> &global_cells_per_rank_)
    : global_cells_per_rank(global_cells_per_rank_) {

  // check cell partition
  Require(global_cells_per_rank.size() == static_cast<size_t>(rtt_c4::nodes()));
}

//------------------------------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief split a global mesh into rank-local meshes (assuming the global_mesh is on each rank)
 *
 * \param[in] global_mesh a Draco_Mesh object representing the global mesh (the full domain)
 *
 * \return an MPI rank-local mesh (Draco_Mesh object) for the partition of cells on the rank
 */
std::shared_ptr<Draco_Mesh>
Draco_Mesh_Splitter::split_mesh(const std::shared_ptr<Draco_Mesh> &global_mesh) {

  // extract interface data from global mesh ...
  const std::vector<std::vector<std::vector<unsigned>>> global_cell2node =
      global_mesh->get_cell_to_node_linkage();
  const std::vector<std::vector<double>> global_coord_vec = global_mesh->get_node_coord_vec();

  // use global cell partitioning to create rank-local interface data ...

  // ... get this (my) rank and corresponding global cell partition
  const int my_rank = rtt_c4::node();
  const auto &my_gcells = global_cells_per_rank[my_rank];

  Remember(auto g_it = std::max_element(my_gcells.begin(), my_gcells.end()));
  Remember(auto gcell_max = *g_it);
  Remember(rtt_c4::global_max(gcell_max));
  Check(gcell_max == static_cast<int64_t>(global_mesh->get_num_cells() - 1));

  // ... set number of local cells
  const size_t num_lcells = my_gcells.size();

  // ... set local data (aggregate serialized interface vectors)

  const unsigned dimension = global_mesh->get_dimension();
  const rtt_mesh_element::Geometry geometry = global_mesh->get_geometry();

  std::vector<unsigned> num_faces_per_cell(num_lcells, 0);
  unsigned tot_nface_per_cell = 0;
  for (size_t lcell = 0; lcell < num_lcells; ++lcell) {
    const auto &gcell = my_gcells[lcell];
    num_faces_per_cell[lcell] = static_cast<unsigned>(global_cell2node[gcell].size());
    tot_nface_per_cell += num_faces_per_cell[lcell];
  }
  Check(tot_nface_per_cell >= 2 * num_lcells);

  std::vector<unsigned> nvert_per_face(tot_nface_per_cell, 0);
  size_t fcntr = 0;
  size_t tot_nvert_per_face = 0;
  for (size_t lcell = 0; lcell < num_lcells; ++lcell) {
    const auto &gcell = my_gcells[lcell];
    const auto &nface = num_faces_per_cell[lcell];
    Check(nface == static_cast<unsigned>(global_cell2node[gcell].size()));
    for (unsigned face = 0; face < nface; ++face) {
      nvert_per_face[fcntr] = static_cast<unsigned>(global_cell2node[gcell][face].size());
      tot_nvert_per_face += nvert_per_face[fcntr];
      fcntr++;
    }
  }

  // get set of global node indices (so that entries are unique)
  std::set<unsigned> gnode_set;
  for (size_t lcell = 0; lcell < num_lcells; ++lcell) {
    const auto &gcell = my_gcells[lcell];
    const auto &nface = num_faces_per_cell[lcell];
    for (unsigned face = 0; face < nface; ++face) {
      const auto nvert = global_cell2node[gcell][face].size();
      for (size_t vert = 0; vert < nvert; ++vert) {
        gnode_set.insert(global_cell2node[gcell][face][vert]);
      }
    }
  }

  // ... set global_node_number (draco_mesh ctor arg)
  std::vector<unsigned> global_node_number(gnode_set.begin(), gnode_set.end());
  const auto num_nodes = static_cast<unsigned>(global_node_number.size());

  // ... use global node as key for conveniently setting cell_to_node_linkage
  std::map<unsigned, unsigned> global2local_node;
  for (unsigned lnode = 0; lnode < num_nodes; ++lnode)
    global2local_node[global_node_number[lnode]] = lnode;

  std::vector<unsigned> cell_to_node_linkage(tot_nvert_per_face);
  size_t ncntr = 0;
  for (size_t lcell = 0; lcell < num_lcells; ++lcell) {
    const auto &gcell = my_gcells[lcell];
    const auto &nface = num_faces_per_cell[lcell];
    for (unsigned face = 0; face < nface; ++face) {
      const auto nvert = global_cell2node[gcell][face].size();
      for (size_t vert = 0; vert < nvert; ++vert) {
        cell_to_node_linkage[ncntr] = global2local_node.at(global_cell2node[gcell][face][vert]);
        ncntr++;
      }
    }
  }

  // set coordinates for each local node index (draco_mesh ctor arg)
  std::vector<double> coordinates(num_nodes * dimension);
  for (unsigned lnode = 0; lnode < num_nodes; ++lnode) {
    for (unsigned idim = 0; idim < dimension; ++idim)
      coordinates[dimension * lnode + idim] = global_coord_vec[global_node_number[lnode]][idim];
  }

  // ... invert local cell -> global cell vector to global cell -> local cell map per rank
  const int num_rank = rtt_c4::nodes();
  std::vector<std::map<int64_t, unsigned>> g2lcell_map_per_rank(num_rank);
  for (int rank = 0; rank < num_rank; ++rank) {
    unsigned lcell = 0;
    for (auto gcell : global_cells_per_rank[rank]) {
      g2lcell_map_per_rank[rank][gcell] = lcell;
      lcell++;
    }
  }

  // ... use global cell-cell linkage to populate local ghost data
  const rtt_mesh::Draco_Mesh::Layout gcgc_linkage = global_mesh->get_cc_linkage();

  // count number of ghost cells and ghost face nodes
  size_t num_ghost = 0;
  size_t nvert_per_ghost = 0;
  for (size_t lcell = 0; lcell < num_lcells; ++lcell) {
    const auto &gcell = my_gcells[lcell];
    const auto &nbrs = gcgc_linkage.at(static_cast<unsigned>(gcell));
    for (const auto &nbr : nbrs) {
      if (g2lcell_map_per_rank[my_rank].count(nbr.first) == 0) {
        // increment number of ghost faces and nodes (neighbor must be off-rank)
        num_ghost++;
        nvert_per_ghost += nbr.second.size();
      }
    }
  }

  // ... ghost cell data (draco_mesh ctor args)
  std::vector<unsigned> ghost_cell_type;
  std::vector<int> ghost_cell_number;
  std::vector<int> ghost_cell_rank;
  std::vector<unsigned> ghost_cell_to_node_linkage;
  ghost_cell_type.reserve(num_ghost);
  ghost_cell_number.reserve(num_ghost);
  ghost_cell_rank.reserve(num_ghost);
  ghost_cell_to_node_linkage.reserve(nvert_per_ghost);

  for (size_t lcell = 0; lcell < num_lcells; ++lcell) {
    const auto &gcell = my_gcells[lcell];
    const auto &nbrs = gcgc_linkage.at(static_cast<unsigned>(gcell));
    for (const auto &nbr : nbrs) {
      if (g2lcell_map_per_rank[my_rank].count(nbr.first) == 0) {
        // neigbor global cell index is off-rank; find rank it is on
        for (int rank = 0; rank < num_rank; ++rank) {
          if (g2lcell_map_per_rank[rank].count(nbr.first) > 0) {
            ghost_cell_type.push_back(static_cast<unsigned>(nbr.second.size()));
            ghost_cell_number.push_back(g2lcell_map_per_rank[rank].at(nbr.first) + 1);
            ghost_cell_rank.push_back(rank);
            // converge global node indices to local node indices
            for (auto gnode : nbr.second) {
              const auto lnode = global2local_node.at(gnode);
              ghost_cell_to_node_linkage.push_back(lnode);
            }
            // break: only 1 rank should have this global cell (nbr.first)
            break;
          }
        }
      }
    }
  }

  // ... indicate the global cell partition is (in)valid, w.r.t. the global_mesh
  Check(ghost_cell_type.size() == num_ghost);
  Check(ghost_cell_number.size() == num_ghost);
  Check(ghost_cell_rank.size() == num_ghost);
  Check(ghost_cell_to_node_linkage.size() == nvert_per_ghost);

  // ... use global cell-side linkage
  const rtt_mesh::Draco_Mesh::Layout gcgs_linkage = global_mesh->get_cs_linkage();

  // ... count number of true boundaries (sides)
  size_t num_sides = 0;
  size_t nvert_per_side = 0;
  for (size_t lcell = 0; lcell < num_lcells; ++lcell) {
    const auto &gcell = my_gcells[lcell];
    if (gcgs_linkage.count(static_cast<unsigned>(gcell)) > 0) {
      const auto &bcs = gcgs_linkage.at(static_cast<unsigned>(gcell));
      for (const auto &bc : bcs) {
        num_sides++;
        nvert_per_side += bc.second.size();
      }
    }
  }

  // ... ghost cell data (draco_mesh ctor args)
  std::vector<unsigned> side_set_flag;
  std::vector<unsigned> side_node_count;
  std::vector<unsigned> side_to_node_linkage;
  side_set_flag.reserve(num_sides);
  side_node_count.reserve(num_sides);
  side_to_node_linkage.reserve(nvert_per_side);

  for (size_t lcell = 0; lcell < num_lcells; ++lcell) {
    const auto &gcell = my_gcells[lcell];
    if (gcgs_linkage.count(static_cast<unsigned>(gcell)) > 0) {
      const auto &bcs = gcgs_linkage.at(static_cast<unsigned>(gcell));
      for (const auto &bc : bcs) {
        side_set_flag.push_back(bc.first);
        side_node_count.push_back(static_cast<unsigned>(bc.second.size()));
        // converge global node indices to local node indices
        for (auto gnode : bc.second) {
          const auto lnode = global2local_node.at(gnode);
          side_to_node_linkage.push_back(lnode);
        }
      }
    }
  }

  // ... indicate the global cell partition is (in)valid, w.r.t. the global_mesh
  Check(side_set_flag.size() == num_sides);
  Check(side_node_count.size() == num_sides);
  Check(side_to_node_linkage.size() == nvert_per_side);

  // instantiate rank-local mesh
  const Draco_Mesh local_mesh(dimension, geometry, num_faces_per_cell, cell_to_node_linkage,
                              side_set_flag, side_node_count, side_to_node_linkage, coordinates,
                              global_node_number, nvert_per_face, ghost_cell_type,
                              ghost_cell_to_node_linkage, ghost_cell_number, ghost_cell_rank);

  // return shared ptr to on-rank mesh
  return std::make_shared<Draco_Mesh>(local_mesh);
}

} // namespace rtt_mesh

//------------------------------------------------------------------------------------------------//
// end of mesh/Draco_Mesh_Splitter.hh
//------------------------------------------------------------------------------------------------//
