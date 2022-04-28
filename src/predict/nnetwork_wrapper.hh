//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   predict/nnetwork_wrapper.hh
 * \author Mathew Cleveland
 * \brief  Definition of the neural network backend to enable build specific implementations.
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef predict_nnetwork_wrapper_hh
#define predict_nnetwork_wrapper_hh

#include "ds++/dbc.hh"
#include <string>
#include <vector>
#ifdef LIBTORCH
#include <torch/script.h>
#include <torch/torch.h>
#else
#include "draco_nn.hh"
#endif

namespace rtt_predict {

//================================================================================================//
/*!
 * \class nnetwork_wrapper
 * \brief This supports both a simple draco (when libtorch is not available) and a full libtorch
 *        implementation of the neural network backend
 */
//================================================================================================//

class nnetwork_wrapper {
public:
  //! Load network data
  void load_network(const std::string &net_file);

  //! Generate network prediction
  std::vector<float> predict(std::vector<float> &signal, const size_t input_size,
                             const size_t output_size);

  //! Check if network is valid
  bool valid();

private:
#ifdef LIBTORCH
  //! Touch Network Data
  torch::jit::script::Module net_;
  bool is_valid_ = false;
#else
  //! Draco Network Data
  draco_nn net_;
#endif
};

} // namespace rtt_predict

#endif // predict_nnetwork_wrapper_hh

//------------------------------------------------------------------------------------------------//
// end of predict/nnetwork_wrapper.hh
//------------------------------------------------------------------------------------------------//
