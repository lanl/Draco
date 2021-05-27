//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   predict/nnetwork_wrapper.hh
 * \author Mathew Cleveland
 * \brief  HEADER ONLY Definition of the neural network backend to enable build
 * specific implementations.
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef predict_nnetwork_wrapper_hh
#define predict_nnetwork_wrapper_hh

#include "ds++/dbc.hh"
#include <string>
#include <vector>
#ifdef LIBTORCH
#include <torch/script.h>
#include <torch/torch.h>
#endif

namespace rtt_predict {
#ifdef LIBTORCH
//================================================================================================//
/*!
 * \class nnetwork_wrapper
 * \brief
 *
 * A libtorch implementation of the neural network backend
 * 
 *
 */
//================================================================================================//

class nnetwork_wrapper {
public:
  void load_network(const std::string &net_file) {
    // store tensor format in torch objects
    net_ = torch::jit::load(net_file);
    is_valid_ = true;
  }

  std::vector<float> predict(std::vector<float> &signal, const long int input_size,
                             const long int output_size) {
    Require(is_valid_);
    Require(signal.size() == static_cast<size_t>(input_size * output_size));
    // Assign input to tensor data type
    std::vector<torch::jit::IValue> inputs;
    at::Tensor T_input = torch::from_blob(signal.data(), {1L, output_size, input_size});
    inputs.push_back(T_input);
    // Generate Prediction from pre-loaded network
    at::Tensor T_output = net_.forward(inputs).toTensor();
    // Assign output tensor data type to standard vector
    std::vector<float> output(T_output.data_ptr<float>(),
                              T_output.data_ptr<float>() + output_size);
    return output;
  }

  bool valid() { return is_valid_; }

private:
  //! Touch Network Data
  torch::jit::script::Module net_;
  bool is_valid_ = false;
};
#else // neural network backend stub

//================================================================================================//
/*!
 * \class nnetwork_wrapper
 * \brief
 *
 * A stub variation of a neural network backend
 *
 */
//================================================================================================//

class nnetwork_wrapper {
public:
  //! load_network stub
  void load_network(const std::string & /*net_file*/){};

  //! predict stub
  std::vector<float> predict(std::vector<float> & /*signal*/, const long int /*input_size*/,
                             const long int /*output_size*/) {
    Insist(false, "Neural Network Backend was not compiled into this build");
    return std::vector<float>{-1.0};
  }

  //! valid stub
  bool valid() { return false; }
};

#endif

} // namespace rtt_predict

#endif // predict_nnetwork_wrapper_hh

//------------------------------------------------------------------------------------------------//
// end of predict/nnetwork_wrapper.hh
//------------------------------------------------------------------------------------------------//
