//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   predict/nnetwork_wrapper.cc
 * \author Mathew Cleveland
 * \brief Definition of the neural network backend to enable build specific implementations.
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "nnetwork_wrapper.hh"

namespace rtt_predict {

#ifdef LIBTORCH

//------------------------------------------------------------------------------------------------//
// Pytorch Network calls
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/*!
 * \brief Load a py-torch neural network file (*.pt)
 *
 * \param[in] net_file py-torch NN file name
 *
 */
//================================================================================================//
void nnetwork_wrapper::load_network(const std::string &net_file) {
  // store tensor format in torch objects
  net_ = torch::jit::load(net_file + ".pt");
  is_valid_ = true;
}

//================================================================================================//
/*!
 * \brief Generate a prediction from a py-torch NN
 *
 * \param[in,out] signal strided input data used to generate a prediction (input_size*output_size)
 * \param[in] input_size number of input nodes per prediction
 * \param[in] output_size number of input signals included in the signal vector
 *
 * \return prediction for each input signal
 */
//================================================================================================//
std::vector<float> nnetwork_wrapper::predict(std::vector<float> &signal, const size_t input_size,
                                             const size_t output_size) {
  Require(is_valid_);
  Require(signal.size() == input_size * output_size);
  // Assign input to tensor data type
  std::vector<torch::jit::IValue> inputs;
  at::Tensor T_input = torch::from_blob(
      signal.data(), {1L, static_cast<long int>(output_size), static_cast<long int>(input_size)});
  inputs.push_back(T_input);
  // Generate Prediction from pre-loaded network
  at::Tensor T_output = net_.forward(inputs).toTensor();
  // Assign output tensor data type to standard vector
  std::vector<float> output(T_output.data_ptr<float>(), T_output.data_ptr<float>() + output_size);
  return output;
}
//! Check if NN is valid
bool nnetwork_wrapper::valid() { return is_valid_; }

#else

//------------------------------------------------------------------------------------------------//
// Draco Neural Network
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/*!
 * \brief Load a draco neural network file (*.nnb)
 *
 * \param[in] net_file draco NN binary file
 *
 */
//================================================================================================//
void nnetwork_wrapper::load_network(const std::string &net_file) {
  net_.load_network(net_file + ".nnb");
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
std::vector<float> nnetwork_wrapper::predict(std::vector<float> &signal, const size_t input_size,
                                             const size_t output_size) {
  return net_.predict(signal, input_size, output_size);
}

//! Check if NN is valid
bool nnetwork_wrapper::valid() { return net_.valid(); }

#endif

} // namespace rtt_predict

//------------------------------------------------------------------------------------------------//
// end of predict/nnetwork_wrapper.cc
//------------------------------------------------------------------------------------------------//
