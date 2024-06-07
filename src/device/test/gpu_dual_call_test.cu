///--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   device/test/gpu_dual_call_test.cc
 * \author Alex R. Long
 * \date   Mon Mar 25 2019
 * \brief  Show how code can be called from GPU and host
 * \note   Copyright (C) 2019-2020 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "dsxx/Release.hh"
#include "dsxx/ScalarUnitTest.hh"
#include "device/GPU_Device.hh"
#include "device/test/Dual_Call.hh"

#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using namespace rtt_device;
using namespace rtt_device_test;

int dual_call_test(rtt_dsxx::ScalarUnitTest &ut) {

  int n_cells = 1029;
  vector<double> src_cell_bias(n_cells, 1.0);
  vector<double> e_field(n_cells, 1.0);
  vector<int> n_field(n_cells, 0);
  unsigned long long device_n_tot = 0;
  unsigned long long host_n_tot = 0;

  const double part_per_e = 1.0;
  const unsigned max_particles_pspc = 100;

  constexpr int threads_per_block = 512;
  int n_blocks = (n_cells + threads_per_block - 1) / threads_per_block;
  // setup and copy all fields
  vector<unsigned long long> n_tot_block(n_blocks, 0);
  int *D_n_field = NULL;
  double *D_e_field = NULL;
  double *D_src_cell_bias = NULL;
  unsigned long long *D_n_tot = NULL;
  rtt_device::GPU_Device gpu;
  cudaError_t err = cudaMalloc((void **)&D_n_field, n_cells * sizeof(int));
  std::cout << gpu.getErrorMessage(cudaGetLastError()) << std::endl;
  err = cudaMalloc((void **)&D_e_field, n_cells * sizeof(double));
  std::cout << gpu.getErrorMessage(cudaGetLastError()) << std::endl;
  err = cudaMalloc((void **)&D_src_cell_bias, n_cells * sizeof(double));
  std::cout << gpu.getErrorMessage(cudaGetLastError()) << std::endl;
  err = cudaMalloc((void **)&D_n_tot, n_blocks * sizeof(unsigned long long));
  std::cout << gpu.getErrorMessage(cudaGetLastError()) << std::endl;
  err = cudaMemcpy(D_n_field, &n_field[0], n_cells * sizeof(int), cudaMemcpyHostToDevice);
  std::cout << gpu.getErrorMessage(cudaGetLastError()) << std::endl;
  err = cudaMemcpy(D_e_field, &e_field[0], n_cells * sizeof(double), cudaMemcpyHostToDevice);
  std::cout << gpu.getErrorMessage(cudaGetLastError()) << std::endl;
  err = cudaMemcpy(D_src_cell_bias, &src_cell_bias[0], n_cells * sizeof(double),
                   cudaMemcpyHostToDevice);
  std::cout << gpu.getErrorMessage(cudaGetLastError()) << std::endl;

  dim3 blockSize;
  dim3 gridSize;
  blockSize.x = n_blocks;
  blockSize.y = 1;
  blockSize.z = 1;
  gridSize.x = threads_per_block;
  gridSize.y = 1;
  gridSize.z = 1;
  cuda_conserve_calc_num_src_particles<<<blockSize, gridSize>>>(
      part_per_e, max_particles_pspc, n_cells, D_e_field, D_src_cell_bias, D_n_field, D_n_tot);
  std::cout << gpu.getErrorMessage(cudaGetLastError()) << std::endl;
  err = cudaDeviceSynchronize();
  std::cout << gpu.getErrorMessage(cudaGetLastError()) << std::endl;

  vector<double> e_field_out(n_cells, 0.0);
  err = cudaMemcpy(&n_field[0], D_n_field, n_cells * sizeof(int), cudaMemcpyDeviceToHost);
  std::cout << gpu.getErrorMessage(err) << std::endl;

  err = cudaMemcpy(&n_tot_block[0], D_n_tot, n_blocks * sizeof(unsigned long long),
                   cudaMemcpyDeviceToHost);
  std::cout << gpu.getErrorMessage(err) << std::endl;

  err = cudaDeviceReset();
  std::cout << gpu.getErrorMessage(err) << std::endl;

  err = cudaFree(D_n_field);
  std::cout << gpu.getErrorMessage(cudaGetLastError()) << std::endl;
  err = cudaFree(D_e_field);
  std::cout << gpu.getErrorMessage(cudaGetLastError()) << std::endl;
  err = cudaFree(D_src_cell_bias);
  std::cout << gpu.getErrorMessage(cudaGetLastError()) << std::endl;
  err = cudaFree(D_n_tot);
  std::cout << gpu.getErrorMessage(cudaGetLastError()) << std::endl;
  // reduce the n_tot over all thread blocks
  for (int i = 0; i < n_blocks; ++i)
    device_n_tot += n_tot_block[i];

  cout << "N total: " << device_n_tot << endl;

  host_n_tot = sub_conserve_calc_num_src_particles(part_per_e, max_particles_pspc, 0, n_cells,
                                                   &e_field[0], &src_cell_bias[0], &n_field[0]);

  cout << "Host N total: " << host_n_tot << endl;
  if (host_n_tot != device_n_tot)
    FAILMSG(string("Host and device totals don't match!"));
  else
    PASSMSG("Host and device calls agree on value.");

  return 0;
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    dual_call_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of device/Dual_Call.cc
//------------------------------------------------------------------------------------------------//
