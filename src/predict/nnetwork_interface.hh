//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   predict/nnetwork_interface.hh
 * \author Mathew Cleveland
 * \brief  Provide a basic neural network interface. Given there are many open source NN libraries
 *         we decided to provide a general wrapper class to provide flexibility in the
 *         implementation details.
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef predict_nnetwork_interface_hh
#define predict_nnetwork_interface_hh

#include "nnetwork_wrapper.hh"
#include <string>
#include <vector>

namespace rtt_predict {

//================================================================================================//
/*!
 * \class nnetwork_interface
 * \brief
 *
 * Provide a flexible interface for external Neural Network software for generating inline
 * predictions from external Machine Learning libraries
 *
 * \sa nnetwork_interface.cc for detailed descriptions.
 */
//================================================================================================//

class nnetwork_interface : private nnetwork_wrapper {
public:
  //! Default constructors.
  explicit nnetwork_interface(const std::string &network_state_file);

  //! Predict outcome
  std::vector<float> predict(std::vector<float> &input, const size_t input_size,
                             const size_t output_size);

  //! Check that the backend is valid
  bool valid() { return nnetwork_wrapper::valid(); }
};

} // namespace rtt_predict

#endif // predict_nnetwork_interface_hh

//------------------------------------------------------------------------------------------------//
// end of predict/nnetwork_interface.hh
//------------------------------------------------------------------------------------------------//
