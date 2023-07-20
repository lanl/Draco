//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/global_containers.t.hh
 * \author Kent Budge
 * \date   Mon Mar 24 09:26:31 2008
 * \brief  Member definitions of class global_containers
 * \note   Copyright (C) 2021-2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef c4_global_containers_t_hh
#define c4_global_containers_t_hh

#include "c4/config.h"

#ifdef C4_MPI

#include "C4_Functions.hh"
#include "C4_Req.hh"
#include "gatherv.hh"
#include "global_containers.hh"
#include "ds++/Assert.hh"

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
// Helper traits classes to reduce code duplication
//------------------------------------------------------------------------------------------------//
template <typename T> struct global_merge_traits { using ET = T; };
// Specialization for bool (Element type is 'int' instead of bool.)
template <> struct global_merge_traits<bool> { using ET = int; };

// void global_merge_setup(){}

//------------------------------------------------------------------------------------------------//
/*!
 * Merge a set across all processors.
 *
 * \param local_set On entry, contains a local set. On exit, contains a set consisting of the union
 *        of all the local sets that came into the function on all processors.
 */
template <typename ElementType> void global_merge(std::set<ElementType> &local_set) {
  using namespace std;

  // Break out promptly if not running in parallel.
  unsigned const number_of_processors = nodes();
  if (number_of_processors < 2)
    return;

  // Flatten the sets
  vector<ElementType> local_elements(local_set.size());
  copy(local_set.begin(), local_set.end(), local_elements.begin());

  // Gather the sets
  vector<vector<ElementType>> global_elements;
  indeterminate_gatherv(local_elements, global_elements);

  if (rtt_c4::node() == 0) {
    Check(global_elements.size() == number_of_processors);
    for (unsigned p = 1; p < number_of_processors; ++p) {
      Check(global_elements[p].size() < UINT_MAX);
      auto const count = static_cast<unsigned>(global_elements[p].size());
      for (unsigned i = 0; i < count; ++i) {
        local_set.insert(global_elements[p][i]);
      }
    }
  }

  Check(local_set.size() < UINT_MAX);
  auto number_of_elements = static_cast<unsigned>(local_set.size());
  broadcast(&number_of_elements, 1, 0);

  local_elements.resize(number_of_elements);
  if (node() == 0) {
    copy(local_set.begin(), local_set.end(), local_elements.begin());
  }
  broadcast(&local_elements[0], number_of_elements, 0);

  if (node() != 0) {
    for (unsigned i = 0; i < number_of_elements; ++i) {
      local_set.insert(local_elements[i]);
    }
  }
}

//------------------------------------------------------------------------------------------------//
template <typename IndexType, typename ElementType>
void global_merge(std::map<IndexType, ElementType> &local_map) {
  using namespace std;
  using ET = typename global_merge_traits<ElementType>::ET;

  size_t const number_of_processors = nodes();
  if (number_of_processors < 2)
    return;

  // Flatten the maps
  vector<IndexType> local_indices(local_map.size());
  vector<ET> local_elements(local_map.size());

  size_t j(0);
  for (auto it = local_map.begin(); it != local_map.end(); ++it, ++j) {
    local_indices[j] = it->first;
    local_elements[j] = static_cast<ET>(it->second); // for bool, cast to int
  }

  // Gather the indices
  vector<vector<IndexType>> global_indices;
  indeterminate_gatherv(local_indices, global_indices);

  // Gather the elements
  vector<vector<ET>> global_elements(number_of_processors);
  for (size_t ip = 0; ip < number_of_processors; ++ip) {
    global_elements[ip].resize(global_indices[ip].size());
  }
  determinate_gatherv(local_elements, global_elements);

  vector<IndexType> index;
  vector<ET> elements;
  unsigned number_of_elements;
  if (node() == 0) {
    for (size_t p = 1; p < number_of_processors; ++p) {
      vector<IndexType> const &other_index = global_indices[p];
      vector<ET> const &other_elements = global_elements[p];
      auto const number_of_other_elements = other_index.size();
      Check(other_index.size() == other_elements.size());
      for (size_t k = 0; k < number_of_other_elements; ++k) {
        // cast int back to bool when ElementType is bool
        local_map.insert(pair<IndexType, ElementType>(other_index[k],
                                                      static_cast<ElementType>(other_elements[k])));
      }
    }
    Check(local_map.size() < UINT_MAX);
    number_of_elements = static_cast<unsigned>(local_map.size());
    index.resize(number_of_elements);
    elements.resize(number_of_elements);
    j = 0;
    for (auto it = local_map.begin(); it != local_map.end(); ++it, ++j) {
      index[j] = it->first;
      elements[j] = static_cast<ET>(it->second); // for bool, cast to int.
    }
  }

  broadcast(&number_of_elements, 1, 0);

  index.resize(number_of_elements);
  elements.resize(number_of_elements);

  broadcast(number_of_elements ? &index[0] : nullptr, number_of_elements, 0);
  broadcast(number_of_elements ? &elements[0] : nullptr, number_of_elements, 0);

  if (node() != 0) {
    for (unsigned k = 0; k < number_of_elements; ++k) {
      // cast int back to bool when ElementType is bool
      local_map.insert(
          pair<IndexType, ElementType>(index[k], static_cast<ElementType>(elements[k])));
    }
  }
}

} // end namespace rtt_c4

#endif // C4_MPI

#endif // c4_global_containers_t_hh

//------------------------------------------------------------------------------------------------//
// end of c4/global_containers.t.hh
//------------------------------------------------------------------------------------------------//
