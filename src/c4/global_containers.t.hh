//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/global_containers.t.hh
 * \author Kent Budge
 * \date   Mon Mar 24 09:26:31 2008
 * \brief  Member definitions of class global_containers
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved. */
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
/*!
 * Merge a set across all processors.
 *
 * \param local_set On entry, contains a local set. On exit, contains a set consisting of the union
 *        of all the local sets that came into the function on all processors.
 */
template <typename ElementType> void global_merge(std::set<ElementType> &local_set) {
  using namespace std;

  // Break out promptly if not running in parallel.

  unsigned number_of_processors = nodes();
  if (number_of_processors < 2)
    return;

  // Flatten the sets
  Check(local_set.size() < UINT_MAX);
  auto const number_of_local_elements = static_cast<unsigned>(local_set.size());
  vector<ElementType> local_elements;
  local_elements.resize(number_of_local_elements);
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

  unsigned number_of_processors = nodes();
  if (number_of_processors < 2)
    return;

  // Flatten the maps
  Check(local_map.size() < UINT_MAX);
  auto const number_of_local_elements = static_cast<unsigned>(local_map.size());
  vector<IndexType> local_indices(number_of_local_elements);
  vector<ElementType> local_elements(number_of_local_elements);
  unsigned j;
  typename map<IndexType, ElementType>::const_iterator i;
  for (i = local_map.begin(), j = 0; i != local_map.end(); ++i, ++j) {
    local_indices[j] = i->first;
    local_elements[j] = i->second;
  }

  // Gather the indices
  vector<vector<IndexType>> global_indices;
  indeterminate_gatherv(local_indices, global_indices);

  // Gather the elements
  vector<vector<ElementType>> global_elements(number_of_processors);
  for (unsigned ip = 0; ip < number_of_processors; ++ip) {
    global_elements[ip].resize(global_indices[ip].size());
  }
  determinate_gatherv(local_elements, global_elements);

  unsigned number_of_elements;
  vector<IndexType> index;
  vector<ElementType> elements;
  if (node() == 0) {
    for (unsigned p = 1; p < number_of_processors; ++p) {
      vector<IndexType> const &other_index = global_indices[p];
      vector<ElementType> const &other_elements = global_elements[p];
      Check(other_index.size() < UINT_MAX);
      auto const number_of_other_elements = static_cast<unsigned>(other_index.size());
      Check(other_index.size() == other_elements.size());
      for (unsigned k = 0; k < number_of_other_elements; ++k) {
        local_map.insert(pair<IndexType, ElementType>(other_index[k], other_elements[k]));
      }
    }
    Check(local_map.size() < UINT_MAX);
    number_of_elements = static_cast<unsigned>(local_map.size());
    index.resize(number_of_elements);
    elements.resize(number_of_elements);

    for (i = local_map.begin(), j = 0; i != local_map.end(); ++i, ++j) {
      index[j] = i->first;
      elements[j] = i->second;
    }
  }

  broadcast(&number_of_elements, 1, 0);

  index.resize(number_of_elements);
  elements.resize(number_of_elements);

  broadcast(number_of_elements ? &index[0] : nullptr, number_of_elements, 0);
  broadcast(number_of_elements ? &elements[0] : nullptr, number_of_elements, 0);

  if (node() != 0) {
    for (unsigned k = 0; k < number_of_elements; ++k) {
      local_map.insert(pair<IndexType, ElementType>(index[k], elements[k]));
    }
  }
}

//------------------------------------------------------------------------------------------------//
/* We have specialized the case of bool map elements because the standard C++ STL library does
 * "clever" things with bool containers that don't play well with the generic implementation. In
 * particular, the communications steps promote the bool elements to int to ensure correct
 * communication. Char might work as well and be more efficient; we can experiment with this if this
 * code ever proves a computational bottleneck.
 */
template <typename IndexType> void global_merge(std::map<IndexType, bool> &local_map) {
  using namespace std;

  unsigned number_of_processors = nodes();
  if (number_of_processors < 2)
    return;

  // Flatten the maps, promoting the bool elements to int so they will play
  // well with C4.
  Check(local_map.size() < UINT_MAX);
  auto const number_of_local_elements = static_cast<unsigned>(local_map.size());
  vector<IndexType> local_indices(number_of_local_elements);
  vector<int> local_elements(number_of_local_elements);
  unsigned j;
  typename map<IndexType, bool>::const_iterator i;
  for (i = local_map.begin(), j = 0; i != local_map.end(); ++i, ++j) {
    local_indices[j] = i->first;
    local_elements[j] = i->second;
  }

  // Gather the indices
  vector<vector<IndexType>> global_indices;
  indeterminate_gatherv(local_indices, global_indices);

  // Gather the elements
  vector<vector<int>> global_elements(number_of_processors);
  for (unsigned ip = 0; ip < number_of_processors; ++ip) {
    global_elements[ip].resize(global_indices[ip].size());
  }
  determinate_gatherv(local_elements, global_elements);

  unsigned number_of_elements;
  vector<IndexType> index;
  vector<int> elements;
  if (node() == 0) {
    for (unsigned p = 1; p < number_of_processors; ++p) {
      vector<IndexType> const &other_index = global_indices[p];
      vector<int> const &other_elements = global_elements[p];
      Check(other_index.size() < UINT_MAX);
      auto const number_of_other_elements = static_cast<unsigned>(other_index.size());
      Check(other_index.size() == other_elements.size());
      for (unsigned k = 0; k < number_of_other_elements; ++k) {
        IndexType const &oindex = other_index[k];

        if (local_map.find(oindex) != local_map.end() &&
            local_map[oindex] != static_cast<bool>(other_elements[k])) {
          throw invalid_argument("inconsistent global map");
        }
        local_map[oindex] = other_elements[k];
      }
    }
    Check(local_map.size() < UINT_MAX);
    number_of_elements = static_cast<unsigned>(local_map.size());
    index.resize(number_of_elements);
    elements.resize(number_of_elements);

    for (i = local_map.begin(), j = 0; i != local_map.end(); ++i, ++j) {
      index[j] = i->first;
      elements[j] = i->second;
    }
  }

  broadcast(&number_of_elements, 1, 0);

  index.resize(number_of_elements);
  elements.resize(number_of_elements);
  broadcast(&index[0], number_of_elements, 0);
  broadcast(&elements[0], number_of_elements, 0);

  // Build the final map, converting the ints back to bool.
  if (node() != 0) {
    for (unsigned k = 0; k < number_of_elements; ++k) {
      local_map[index[k]] = static_cast<bool>(elements[k]);
    }
  }
}

} // end namespace rtt_c4

#endif // C4_MPI

#endif // c4_global_containers_t_hh

//------------------------------------------------------------------------------------------------//
// end of c4/global_containers.t.hh
//------------------------------------------------------------------------------------------------//
