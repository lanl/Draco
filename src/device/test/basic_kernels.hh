///--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   device/test/basic_kernels.hh
 * \author Alex R. Long
 * \date   Mon Mar 25 2019
 * \brief  Simple kernels for basic GPU tests
 * \note   Copyright (C) 2019-2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_device_test_basic_kernels_hh
#define rtt_device_test_basic_kernels_hh

#include "device/device_gpu.h"
#include <algorithm>
#include <numeric>
#include <vector>

namespace rtt_device_test {

__global__ void vector_add(double const *A_dev, double const *B_dev, double *C_dev, int const N);
__global__ void sum(int *dest, int a, int b);

} // namespace rtt_device_test

#endif // rtt_device_test_basic_kernels_hh

//------------------------------------------------------------------------------------------------//
// end of device/test/basic_kernels.hh
//------------------------------------------------------------------------------------------------//
