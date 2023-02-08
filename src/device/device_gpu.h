//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   device/device_gpu.h
 * \author Kelly (KT) Thompson
 * \brief  Wrap the cuda.h or hipruntim.h headers while preventing comiler warnings about vendor
 *         code.
 * \note   Copyright (C) 2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef device_device_gpu_h
#define device_device_gpu_h

#include <device/config.h>

// Suppresses warnings found in "cuda.h" and "hip_runtime.h".
// http://wiki.services.openoffice.org/wiki/Writing_warning-free_code#When_all_else_fails
#if defined __GNUC__
#pragma GCC system_header
// Intel defines __GNUC__ by default
#ifdef __INTEL_COMPILER
#pragma warning push
#endif
#elif defined __SUNPRO_CC
#pragma disable_warn
#elif defined _MSC_VER
#pragma warning(push, 1)
#endif

#ifdef USE_HIP
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
#define cudaError_t hipError_t
#define cudaGetDevice hipGetDevice
#define cudaSetDevice hipSetDevice
#define cudaGetDeviceCount hipGetDeviceCount
#define cudaGetDeviceProperties hipGetDeviceProperties
#define cudaDeviceProp hipDeviceProp_t
#define cudaDevAttrTotalConstantMemory hipDeviceAttributeTotalConstantMemory
#define cudaDeviceGetAttribute hipDeviceGetAttribute
#define cudaDevAttrMaxRegistersPerBlock hipDeviceAttributeMaxRegistersPerBlock
#define cudaDevAttrClockRate hipDeviceAttributeClockRate
#define cudaDevAttrTextureAlignment hipDeviceAttributeTextureAlignment
#define cudaSuccess hipSuccess
#define cudaGetErrorString hipGetErrorString
#define cudaMalloc hipMalloc
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#define cudaFree hipFree
#define cudaMemcpy hipMemcpy
#define cudaDeviceSynchronize hipDeviceSynchronize
#define cudaGetLastError hipGetLastError
#define cudaDeviceReset hipDeviceReset
#endif

#ifdef USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif

#if defined __GNUC__
#pragma GCC system_header
#ifdef __INTEL_COMPILER
#pragma warning pop
#endif
#elif defined __SUNPRO_CC
#pragma enable_warn
#elif defined _MSC_VER
#pragma warning(pop)
#endif

#endif // device_device_gpu_h

//------------------------------------------------------------------------------------------------//
// end of device/device_gpu.h
//------------------------------------------------------------------------------------------------//
