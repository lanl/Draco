///--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   device/test/Dual_Call.cu
 * \author Alex R. Long
 * \date   Mon Mar 25 2019
 * \brief  Show how code can be called from GPU and host
 * \note   Copyright (C) 2019-2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Dual_Call.hh"
#include "ds++/Assert.hh"

namespace rtt_device_test {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Calculate the number of source particles for a range of cells
 *
 * \param[in] part_per_e energy for this source
 * \param[in] max_particles_pspc max partices per species
 * \param[in] cell_start starting cell index
 * \param[in] cell_end ending cell index
 * \param[in] e_field energy in a cell
 * \param[in] src_cell_bias bias in a cell
 * \param[in,out] n_field destination for particles in a cell
 * \param[out] return number of particles over this cell range
 */
__host__ __device__ unsigned long long sub_conserve_calc_num_src_particles(
    const double part_per_e, unsigned max_particles_pspc, const size_t cell_start,
    const size_t cell_end, const double *e_field, const double *src_cell_bias, int *n_field) {
  unsigned long long ntot = 0;

  ntot = 0;

  Insist_device(part_per_e > 0, "Use the Insist_device function");

  // sweep through cells and calculate number of particles per cell
  for (size_t cell = cell_start; cell < cell_end; cell++) {
    // if the cell has any energy try to put some particles in it
    if (e_field[cell] > 0.0) {
      // get estimate of number of particles per cell to nearest integer per species, a cell-based
      // bias can be added that simply multiplies the expected number by a user defined bias; the
      // energy balance will still be correct because particles will simply be subtracted from other
      // cells to compensate
      const double d_num = e_field[cell] * part_per_e * src_cell_bias[cell];

      // We are about to cast d_num back to int.  Ensure that the conversion is valid.  If not, set
      // the number of particles to the ceiling value provided in Source.hh.
      if (d_num < static_cast<double>(max_particles_pspc - 1)) {
        n_field[cell] = static_cast<int>(d_num + 0.5);

        // try to get at least one particle per cell per species
        if (n_field[cell] == 0)
          n_field[cell] = 1;

      } else {
        n_field[cell] = max_particles_pspc;
      }

      // increment particle counter (uint64_t += int)
      ntot += n_field[cell];
    } else
      n_field[cell] = 0;
  }
  return ntot;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Launch a kernel to calculate the number of source particles
 *
 * \param[in] part_per_e energy for this source
 * \param[in] max_particles_pspc max partices per species
 * \param[in] cont_size size of all fields
 * \param[in] e_field energy in a cell
 * \param[in] src_cell_bias bias in a cell
 * \param[in,out] n_field destination for particles in a cell
 * \param[in,out] ntot total particles per thread block
 */
__global__ void cuda_conserve_calc_num_src_particles(const double part_per_e,
                                                     unsigned max_particles_pspc, int cont_size,
                                                     const double *e_field,
                                                     const double *src_cell_bias, int *n_field,
                                                     unsigned long long *ntot) {

  __shared__ unsigned long long shared_data[512];
  size_t cell_start = threadIdx.x + blockIdx.x * blockDim.x;
  size_t cell_end = cell_start + 1;
  if (cell_start < cont_size) {
    shared_data[threadIdx.x] = sub_conserve_calc_num_src_particles(
        part_per_e, max_particles_pspc, cell_start, cell_end, e_field, src_cell_bias, n_field);
  } else
    shared_data[threadIdx.x] = 0;
  __syncthreads();
  for (unsigned int s = 1; s < blockDim.x; ++s) {
    if (s == threadIdx.x)
      shared_data[0] += shared_data[s];
    __syncthreads();
  }
  __syncthreads();
  ntot[blockIdx.x] = shared_data[0];
}

} // namespace rtt_device_test

//------------------------------------------------------------------------------------------------//
// end of device/test/Dual_Call.cc
//------------------------------------------------------------------------------------------------//
