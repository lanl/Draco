//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   predict/nnetwork_wrapper.hh
 * \author Mathew Cleveland
 * \date   Feb. 16th 2022
 * \brief  Simple draco neural network class implementation
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "draco_nn.hh"
#include "ds++/dbc.hh"
#include <fstream>
#include <numeric>
#include <string>
#include <vector>

namespace rtt_predict {

//================================================================================================//
/*!
 * \brief Load a draco neural network file (*.nnb)
 *
 * \param[in] net_file draco NN binary file
 */
//================================================================================================//
void draco_nn::load_network(const std::string &net_file) {
  std::ifstream in(net_file, std::ofstream::binary);
  Insist(in.good(), "Draco NN file -> " + net_file + " dose not exist!!!");
  int test_int;
  in.read(reinterpret_cast<char *>(&test_int), sizeof(int));
  Insist(test_int == 1,
         "Test integer in binary file should be equal to one check *.nnb binary file "
         "for endianness and correct binary format");
  int ibuffer;
  in.read(reinterpret_cast<char *>(&ibuffer), sizeof(int));
  n_layers = static_cast<size_t>(ibuffer);
  nn_layers.resize(n_layers);
  layer_shape.resize(n_layers);
  weights.resize(n_layers);
  bias.resize(n_layers);
  std::array<int, 5> layer_buffer{};
  for (size_t n = 0; n < n_layers; n++) {
    in.read(reinterpret_cast<char *>(&layer_buffer[0]), sizeof(layer_buffer));
    Insist(layer_buffer[3] == layer_buffer[4],
           "We currently only support bias sizes equal to the output size for linear layers");
    nn_layers[n] = std::pair<ACTIVATION_LAYER_TYPE, NN_LAYER_TYPE>(
        static_cast<ACTIVATION_LAYER_TYPE>(layer_buffer[0]),
        static_cast<NN_LAYER_TYPE>(layer_buffer[1]));
    layer_shape[n] = std::pair<size_t, size_t>(static_cast<size_t>(layer_buffer[2]),
                                               static_cast<size_t>(layer_buffer[3]));
    weights[n] = std::vector<std::vector<float>>(layer_shape[n].second,
                                                 std::vector<float>(layer_shape[n].first));
    bias[n] = std::vector<float>(layer_buffer[4]);
  }
  // fill up bias and weights
  for (size_t n = 0; n < n_layers; n++) {
    Check(bias[n].size() == layer_shape[n].second);
    const size_t weight_buffer_size = sizeof(float) * layer_shape[n].first;
    const size_t bias_buffer_size = sizeof(float) * layer_shape[n].second;
    for (auto &weight : weights[n]) {
      Check(weight.size() == layer_shape[n].first);
      in.read(reinterpret_cast<char *>(&weight[0]), weight_buffer_size);
    }
    in.read(reinterpret_cast<char *>(&bias[n][0]), bias_buffer_size);
  }
  valid_ = true;
}

//================================================================================================//
/*!
 * \brief Generate a prediction from an input signal using the simple draco NN.
 *
 * \param[in,out] signal strided input data used to generate a prediction (input_size*output_size)
 * \param[in] input_size number of input nodes per prediction
 * \param[in] output_size number of input signals included in the signal vector
 *
 * \return prediction for each input signal
 */
//================================================================================================//
std::vector<float> draco_nn::predict(std::vector<float> &signal, const size_t Remember(input_size),
                                     const size_t output_size) {
  Require(signal.size() == (input_size * output_size));
  std::vector<float> input = signal;
  std::vector<float> result;
  for (size_t n = 0; n < n_layers; n++) {
    Insist(nn_layers[n].second == NN_LAYER_TYPE::LINEAR, "Currently only supporting Linear Layers");
    Insist(nn_layers[n].first < ACTIVATION_LAYER_TYPE::ACTIVATION_END,
           "Specified unsported activation type");
    // apply activation to input if necessary
    if (nn_layers[n].first == ACTIVATION_LAYER_TYPE::RELU)
      std::transform(input.begin(), input.end(), input.begin(),
                     [](auto i) { return std::max(i, 0.0F); });
    // resize the result
    result.resize(layer_shape[n].second * output_size);
    // loop over all signals
    size_t r = 0;
    for (size_t si = 0; si < output_size * layer_shape[n].first; si += layer_shape[n].first) {
      size_t b = 0;
      for (auto &weight : weights[n]) {
        result[r] =
            std::inner_product(weight.begin(), weight.end(), input.begin() + si, bias[n][b]);
        b++;
        r++;
      }
    }
    input = result;
  }
  return result;
}

} // namespace rtt_predict

//------------------------------------------------------------------------------------------------//
// end of predict/draco_nn.cc
//------------------------------------------------------------------------------------------------//
