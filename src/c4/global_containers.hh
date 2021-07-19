//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   c4/global_containers.hh
 * \author Kent Budge
 * \brief  Define class global_containers
 * \note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef c4_global_containers_hh
#define c4_global_containers_hh

#include <map>
#include <set>

namespace rtt_c4 {

//------------------------------------------------------------------------------------------------//
/*! Merge local sets into a global set
 *
 * This template takes a processor-local set and merges its elements globally. After this operation
 * is perfomed, each local container will contain all elements found on any processor.
 *
 * /param local_set Local set to be globally merged.
 */
template <typename ElementType> void global_merge(std::set<ElementType> &local_set);

//------------------------------------------------------------------------------------------------//
/*! Merge local maps into a global map
 *
 * This template takes a processor-local maps and merges its elements globally. After this operation
 * is perfomed, each local container will contain all elements found on any processor.
 *
 * /param local_set Local set to be globally merged.
 */
template <typename IndexType, typename ElementType>
void global_merge(std::map<IndexType, ElementType> &local_map);

//------------------------------------------------------------------------------------------------//
//! Specialization for bool
template <typename IndexType> void global_merge(std::map<IndexType, bool> &local_map);

} // end namespace rtt_c4

#endif // c4_global_containers_hh

//------------------------------------------------------------------------------------------------//
// end of c4/global_containers.hh
//------------------------------------------------------------------------------------------------//
