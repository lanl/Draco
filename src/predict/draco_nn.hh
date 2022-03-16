//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   predict/draco_nn.hh
 * \author Mathew Cleveland
 * \brief  Simple draco neural network back end class
 * specific implementations.
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef predict_draco_nn_hh
#define predict_draco_nn_hh

#include <string>
#include <vector>

namespace rtt_predict {

//! Supported Activation Types
enum class ACTIVATION_LAYER_TYPE { NONE, RELU, ACTIVATION_END };
//! Support NN Layers
enum class NN_LAYER_TYPE { LINEAR, LAYER_END };

//================================================================================================//
/*!
 * \brief draco_nn provides a simple neural network back end in place of more advanced backends like
 * py-torch. This class has limited support and is mostly intended to serve as a proof of concept.
 * This class reads a "hand-rolled" py_torch binary NN definition. There is a conversion tool
 * (tools/convert_nn_pt_to_nnb.py) to convert from the *.pt pytorch format to the draco format
 * *.nnb.
 *
 */
//================================================================================================//

class draco_nn {
public:
  //! Load the network
  void load_network(const std::string &net_file);

  //! Prediction call
  std::vector<float> predict(std::vector<float> &signal, const size_t input_size,
                             const size_t output_size);

  //! Check if the network is valid
  bool valid() { return valid_; }

private:
  //! Private network data used for predictions
  //! Valid network state
  bool valid_{false};
  //! Number of layers
  size_t n_layers{0};
  //! Layer activation and layer type (activation, type)
  std::vector<std::pair<ACTIVATION_LAYER_TYPE, NN_LAYER_TYPE>> nn_layers;
  //! Layer shape (input_size, output_size)
  std::vector<std::pair<size_t, size_t>> layer_shape;
  //! NN Weights
  std::vector<std::vector<std::vector<float>>> weights;
  //! NN bias
  std::vector<std::vector<float>> bias;
};

} // namespace rtt_predict

#endif // predict_draco_nn_hh

//------------------------------------------------------------------------------------------------//
// end of predict/draco_nn.hh
//------------------------------------------------------------------------------------------------//
