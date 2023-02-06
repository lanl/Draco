///--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   device/test/Dual_Call.hh
 * \author Alex R. Long
 * \date   Mon Mar 25 2019
 * \brief  Show how code can be called from GPU and host
 * \note   Copyright (C) 2019-2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_device_test_Dual_Call_hh
#define rtt_device_test_Dual_Call_hh

#include "device/device_gpu.h"

namespace rtt_device_test {

__host__ __device__ unsigned long long sub_conserve_calc_num_src_particles(
    const double part_per_e, unsigned max_particles_pspc, const size_t cell_start,
    const size_t cell_end, const double *e_field, const double *src_cell_bias, int *n_field);

__global__ void cuda_conserve_calc_num_src_particles(const double part_per_e,
                                                     unsigned max_particles_pspc, int cont_size,
                                                     const double *e_field,
                                                     const double *src_cell_bias, int *n_field,
                                                     unsigned long long *ntot);

} // namespace rtt_device_test

#endif // rtt_device_test_Dual_Call_hh

//------------------------------------------------------------------------------------------------//
// end of device/test/Dual_Call.hh
//------------------------------------------------------------------------------------------------//
