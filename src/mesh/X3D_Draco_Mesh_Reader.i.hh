//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/X3D_Draco_Mesh_Reader.i.hh
 * \author Ryan Wollaeger <wollaeger@lanl.gov>, Kendra Long
 * \date   Thursday, Jul 12, 2018, 08:46 am
 * \brief  X3D_Draco_Mesh_Reader class implementation header file.
 * \note   Copyright (C) 2018-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/DracoStrings.hh"

namespace rtt_mesh {

//------------------------------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Generate a map from a block in an x3d file
 *
 * \param[in] block_name name of parsed x3d block
 * \param[in] dist number of string-pairs to skip when looking for an x3d block
 *
 * \return map of mesh data with key of type KT and value of type VT
 */
template <typename KT, typename VT>
std::map<KT, std::vector<VT>> X3D_Draco_Mesh_Reader::map_x3d_block(const std::string &block_name,
                                                                   size_t &dist) {

  // parse x3d block and generate x3d_map
  auto label_first = find_iter_of_key(parsed_pairs, block_name, dist);
  auto label_last = find_iter_of_key(parsed_pairs, "end_" + block_name, dist);

  // add distance in map to exclude parsed file region
  dist += std::distance(label_first, label_last);

  Check(dist > 0);
  Check(label_first < label_last);

  // x3d map to return
  std::map<KT, std::vector<VT>> ret_x3d_map;

  for (auto it = label_first + 1; it < label_last; ++it) {

    std::vector<VT> tmp_vec((*it).second.size());
    size_t i = 0;

    // convert value types from string to VT
    for (auto const &j : (*it).second) {

      // try to convert value type to VT, throw if impossible
      try {
        tmp_vec[i] = rtt_dsxx::parse_number_impl<VT>(j);
      } catch (std::invalid_argument &err) {
        Insist(false, err.what());
      }

      // increment counter
      i++;
    }

    // try to convert to key type KT, throw assertion if impossible
    KT tmp_key = convert_key<KT>((*it).first);

    // insert the new entry
    ret_x3d_map.insert(std::pair<KT, std::vector<VT>>(tmp_key, tmp_vec));
  }

  Ensure(ret_x3d_map.size() > 0);
  return ret_x3d_map;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Convert key to template type.
 *
 * \param[in] skey string key
 *
 * \return numerical key of type "KT"
 */
template <typename KT> KT X3D_Draco_Mesh_Reader::convert_key(const std::string &skey) {

  // try to convert to key type KT, throw assertion if impossible
  KT ret_key;
  try {
    ret_key = rtt_dsxx::parse_number_impl<KT>(skey);
  } catch (std::invalid_argument &err) {
    Insist(false, err.what());
  }

  return ret_key;
}

} // namespace rtt_mesh

//------------------------------------------------------------------------------------------------//
// end of mesh/X3D_Draco_Mesh_Reader.i.hh
//------------------------------------------------------------------------------------------------//
