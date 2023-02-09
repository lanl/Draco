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
#define cudaDevAttrClockRate hipDeviceAttributeClockRate
#define cudaDevAttrMaxRegistersPerBlock hipDeviceAttributeMaxRegistersPerBlock
#define cudaDevAttrPciBusId hipDeviceAttributePciBusId
#define cudaDevAttrTextureAlignment hipDeviceAttributeTextureAlignment
#define cudaDevAttrTotalConstantMemory hipDeviceAttributeTotalConstantMemory
#define cudaDeviceGetAttribute hipDeviceGetAttribute
#define cudaDeviceProp hipDeviceProp_t
#define cudaDeviceReset hipDeviceReset
#define cudaDeviceSynchronize hipDeviceSynchronize
#define cudaError_t hipError_t
#define cudaFree hipFree
#define cudaGetDevice hipGetDevice
#define cudaGetDeviceCount hipGetDeviceCount
#define cudaGetDeviceProperties hipGetDeviceProperties
#define cudaGetErrorString hipGetErrorString
#define cudaGetLastError hipGetLastError
#define cudaMalloc hipMalloc
#define cudaMemcpy hipMemcpy
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#define cudaSetDevice hipSetDevice
#define cudaSuccess hipSuccess
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
