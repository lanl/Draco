//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   device/test/basic_kernels.cu
 * \author Kelly Thompson
 * \brief  Small kernel code for testing GPU Device framework.
 * \note   Copyright (C) 2016-2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

// #include <hip/hip_runtime.h>
#include "basic_kernels.hh"

namespace rtt_device_test {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief CUDA kernel for adding two numbers
 *
 * \param[in,out] dest location to store sum
 * \param[in] a value to add
 * \param[in] b value to add
 */
__global__ void sum(int *dest, int a, int b) {
  // Assuming a single thread, 1x1x1 block, 1x1 grid
  *dest = a + b;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief CUDA kernel for adding two vectors
 *
 * \param[in] A_dev vector to add
 * \param[in] B_dev vector to add
 * \param[in,out] C_dev location to store solution vector
 * \param[in] N length of vectors
 */
__global__ void vector_add(double const *A_dev, double const *B_dev, double *C_dev, int const N) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  /* if(i%512==0)
         *   printf("index %d\n",i); */
  if (i < N)
    C_dev[i] = A_dev[i] + B_dev[i];
}

} // namespace rtt_device_test

//------------------------------------------------------------------------------------------------//
// end of basic_kernels.cu
//------------------------------------------------------------------------------------------------//
