//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   predict/nnetwork_interface.cc
 * \author Mathew Cleveland
 * \date   Jan. 11t 2021
 * \brief  Explicit instansiation of the nnetwork_interface class
 * \note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "nnetwork_interface.hh"

namespace rtt_predict {

//================================================================================================//
/*!
 * \brief Load Neural Network state data
 *
 * \param[in] network_state_file (dependes on the NN backend)
 * 
 */
//================================================================================================//
nnetwork_interface::nnetwork_interface(const std::string &network_state_file) {
  nnetwork_wrapper::load_network(network_state_file);
}
std::vector<float> nnetwork_interface::predict(std::vector<float> &input, const size_t input_size,
                                               const size_t output_size) {
  return nnetwork_wrapper::predict(input, input_size, output_size);
}

} // namespace rtt_predict

//------------------------------------------------------------------------------------------------//
// end of nnetwork_interface.cc
//------------------------------------------------------------------------------------------------//
