//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   config/query_gpu.cu
 * \author Alex Long
 * \brief  Small CUDA code that prints the architecture version, used by CMake
 * \date   Thu Mat 21 15:53:51 2019
 * \note   Copyright (C) 2019-2023 Triad National Security, LLC., All rights reserved.
 *
 * This code is from wagonhelm.github.io/articles/2018-03/detecting-cuda-capability-with-cmake */
//------------------------------------------------------------------------------------------------//

#include <stdio.h>

// This ifdef block is related to the code found at src/device/device_gpu.h.  We need this code
// locally in this file because it is compled before cmake is done creating device/config.h.
#ifdef __HIP__
#include "hip/hip_runtime.h"
#define gpuDeviceProp hipDeviceProp_t
#define gpuGetDeviceProperties hipGetDeviceProperties
#define gpuGetLastError hipGetLastError
#define gpuSuccess hipSuccess
#define gpuGetErrorString hipGetErrorString
#define gpuError_t hipError_t
#else
#define gpuDeviceProp cudaDeviceProp
#define gpuGetDeviceProperties cudaGetDeviceProperties
#define gpuGetLastError cudaGetLastError
#define gpuSuccess cudaSuccess
#define gpuGetErrorString cudaGetErrorString
#define gpuError_t cudaError_t
#endif

int main(int argc, char **argv) {
  gpuDeviceProp dP;
  float min_cc = 3.0;

  int rc = gpuGetDeviceProperties(&dP, 0);
  if (rc != gpuSuccess) {
    gpuError_t error = gpuGetLastError();
    printf("CUDA/HIP error: %s", gpuGetErrorString(error));
    return rc; /* Failure */
  }
  if ((dP.major + (dP.minor / 10)) < min_cc) {
    printf("Min Compute Capability of %2.1f required:  %d.%d found\n Not building CUDA/HIP code",
           min_cc, dP.major, dP.minor);
    return 1; /* Failure */
  } else {
    // printf("-arch=sm_%d%d", dP.major, dP.minor);
    printf("%d%d", dP.major, dP.minor);
    return 0; /* Success */
  }
}

//------------------------------------------------------------------------------------------------//
// end of query_gpu.cu
//------------------------------------------------------------------------------------------------//
