///--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   device/test/gpu_hello_rt_api.cu
 * \author Kelly (KT) Thompson
 * \date   Thu Oct 25 15:28:48 2011
 * \brief  Simple test of the CUDA Runtime API.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "dsxx/Soft_Equivalence.hh"
#include "device/GPU_CheckError.hh"
#include <iostream>
#include <string>
#include <vector>

// GPU kernels
#include "device/test/basic_kernels.hh"

using namespace rtt_device_test;

//------------------------------------------------------------------------------------------------//
// Tests
//------------------------------------------------------------------------------------------------//

void simple_add(rtt_dsxx::ScalarUnitTest &ut) {
  // C = A+B

  std::cout << "\n==> Running test: gpu_hello_rt_api::simple_add(ut)\n" << std::endl;

  // create some data
  int const N(2053);
  size_t const BLOCK_SIZE(512);
  std::vector<double> A(N, 1.0);
  std::vector<double> B(N, 0.1);
  std::vector<double> C(N, 0.0);

  // device info
  dim3 blockSize;
  dim3 gridSize;
  blockSize.x = BLOCK_SIZE;
  blockSize.y = 1;
  blockSize.z = 1;
  gridSize.x = ((N + BLOCK_SIZE - 1) / BLOCK_SIZE);
  gridSize.y = 1;
  gridSize.z = 1;

  // allocate memory on the device and move data to the GPU
  double *A_dev(0), *B_dev(0), *C_dev(0);
  int N_bytes = N * sizeof(double);
  DBS_CHECK_ERROR(cudaMalloc((void **)&A_dev, N_bytes));
  DBS_CHECK_ERROR(cudaMalloc((void **)&B_dev, N_bytes));
  DBS_CHECK_ERROR(cudaMalloc((void **)&C_dev, N_bytes));

  DBS_CHECK_ERROR(cudaMemcpy(A_dev, &A[0], N_bytes, cudaMemcpyHostToDevice));
  DBS_CHECK_ERROR(cudaMemcpy(B_dev, &B[0], N_bytes, cudaMemcpyHostToDevice));

  // Launch the cuda kernel
  vector_add<<<gridSize, blockSize>>>(A_dev, B_dev, C_dev, N);

  DBS_CHECK_ERROR(cudaDeviceSynchronize());
  DBS_CHECK_ERROR(cudaGetLastError());

  // retrieve the solution
  DBS_CHECK_ERROR(cudaMemcpy(&C[0], C_dev, N_bytes, cudaMemcpyDeviceToHost));

  // free device memory
  DBS_CHECK_ERROR(cudaFree(A_dev));
  DBS_CHECK_ERROR(cudaFree(B_dev));
  DBS_CHECK_ERROR(cudaFree(C_dev));
  DBS_CHECK_ERROR(cudaDeviceReset());

  // Check solution
  std::vector<double> const Ref(N, 1.1);
  if (rtt_dsxx::soft_equiv(C.begin(), C.end(), Ref.begin(), Ref.end(), 1.0e-15))
    PASSMSG("CUDA kernel vector add works!");
  else
    FAILMSG("CUDA kernel vector add failed!");

  return;
}

//------------------------------------------------------------------------------------------------//
void copy_and_retrieve(rtt_dsxx::ScalarUnitTest &ut) {

  using namespace std;

  std::cout << "\n==> Running test: gpu_hello_rt_api::copy_and_retrieve(ut)\n" << std::endl;

  // Size of array
  size_t const count(100);

  // Host arrays
  vector<float> inH(count, -1.0);
  vector<float> outH(count, -1.0);

  // Allocate device array
  float *inD = NULL; // inD = NULL;
  cudaError_t ret = cudaMalloc((void **)&inD, count * sizeof(float));
  if (ret != cudaSuccess)
    FAILMSG(string("cudaMalloc returned \"") + cudaGetErrorString(ret) + "\"");
  else
    PASSMSG("Allocated inD array on device.");

  // Initialize host data
  for (size_t i = 0; i < count; ++i)
    inH[i] = i;

  // Copy to device
  ret = cudaMemcpy(inD, &inH[0], count * sizeof(float), cudaMemcpyHostToDevice);
  if (ret != cudaSuccess)
    FAILMSG(string("cudaMemcpy returned \"") + cudaGetErrorString(ret) + "\"");
  else
    PASSMSG("Copied array from host to device.");

  // Copy back to host
  ret = cudaMemcpy(&outH[0], inD, count * sizeof(float), cudaMemcpyDeviceToHost);
  if (ret != cudaSuccess)
    FAILMSG(string("cudaMemcpy returned \"") + cudaGetErrorString(ret) + "\"");
  else
    PASSMSG("Copied array from device to host.");

  // Compare original and result arrays.
  size_t diff(0);
  for (size_t i = 0; i < count; ++i) {
    if (!rtt_dsxx::soft_equiv(inH[i], outH[i], 1.0e-8f)) {
      diff = 1;
      cout << "diff: inH[" << i << "] = " << inH[i] << "; outH[" << i << "] = " << outH[i] << endl;
    }
  }
  if (diff != 0)
    FAILMSG("inH != outH");
  else
    PASSMSG("inH == outH");

  // Clean up
  ret = cudaFree(inD);
  if (ret != cudaSuccess)
    FAILMSG(string("cudaFree returned \"") + cudaGetErrorString(ret) + "\"");
  else
    PASSMSG("Freed allocated memory on device.");

  return;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    copy_and_retrieve(ut);
    simple_add(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of gpu_hello_rt_api.cu
//------------------------------------------------------------------------------------------------//
