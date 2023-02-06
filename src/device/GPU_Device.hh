//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   device/GPU_Device.hh
 * \author Kelly (KT) Thompson
 * \brief  Define class GPU_Device
 * \note   Copyright (C) 2011-2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef device_GPU_Device_hh
#define device_GPU_Device_hh

#include "device_gpu.h"
#include "ds++/Assert.hh"
#include <array>
#include <iostream>
#include <string>
#include <vector>

namespace rtt_device {

//================================================================================================//
/*!
 * \class GPU_Device
 * \brief A wrapper for the CUDA environment.
 *
 * This class encapsulates many features provided by the CUDA programming environment. It can query
 * the hardware for available features, find and load .cubin files for GPU kernel execution, etc.
 *
 * \sa GPU_Device.cc for detailed descriptions.
 *
 * \par Code Sample:
 * \code
 *     rtt_device::GPU_Device gpu;
 *
 *     // Create and then print a summary of the devices found.
 *     std::ostringstream out;
 *     size_t const numDev( gpu.numDevicesAvailable() );
 *     out << "GPU device summary:\n\nNumber of devices found: " << numDev << "\n" << endl;
 *     for( size_t device=0; device<numDev; ++device )
 *         gpu.printDeviceSummary( device, out );
 *
 *   // Print the message to stdout
 *   cout << out.str();
 * \endcode
 *
 * \example device/test/gpu_hello_rt_api.cu
 * Test of GPU_Device for with CUDA runtime API.
 */
//================================================================================================//
class GPU_Device {
public:
  // CREATORS

  //! Default constructors.
  GPU_Device();

  //! Destructor.
  ~GPU_Device();

  // ACCESSORS

  //! How many GPU devices are found on the bus?
  size_t numDevicesAvailable() const { return deviceCount; }
  //! maximum number of threads per block
  int maxThreadsPerBlock(int devId = 0) const { return m_maxthreadsperblock[devId]; }
  //! maxThreadsDim[3] is the maximum sizes of each dimension of a block.
  int maxThreadsDim(int dim, int devId = 0) const { return m_maxthreadsdim[devId][dim]; }
  //! maxGridSize[3] is the maximum sizes of each dimension of a grid;
  int maxGridSize(int dim, int devId = 0) const { return m_maxgridsize[devId][dim]; }
  //! the total amount of shared memory available per block in bytes
  size_t sharedMemPerBlock(size_t devId = 0) const { return m_sharedmemperblock[devId]; }
  //! the total amount of constant memory available on the device in bytes;
  int totalConstantMemory(int devId = 0) const { return m_totalconstantmemory[devId]; }
  //! the warp size
  int SIMDWidth(int devId = 0) const { return m_simdwidth[devId]; }
  /*! the maximum pitch allowed by the memory copy functions that involve memory
   *  regions allocated through cuMemAllocPitch() */
  size_t memPitch(size_t devId = 0) const { return m_mempitch[devId]; }
  //! the total number of registers available per block
  int regsPerBlock(int devId = 0) const { return m_regsperblock[devId]; }
  //! the clock frequency in kilohertz
  int clockRate(int devId = 0) const { return m_clockrate[devId]; }
  /*! the alignment requirement; texture base addresses that are aligned to textureAlign bytes do
   *  not need an offset applied to texture fetches */
  int textureAlign(int devId = 0) const { return m_texturealign[devId]; }

  // SERVICES
  //! Print a summary of idevice's features to ostream out.
  void printDeviceSummary(int const idevice, std::ostream &out = std::cout) const;

  // STATICS
  inline static int align(int offset, int alignment) {
    return (offset + alignment - 1) & ~(alignment - 1);
  }
  //! Check cuda return code and throw an Insist on error.
  static std::string getErrorMessage(cudaError_t const err);
  //! Check the value of the return code for CUDA calls.
  static void checkForCudaError(cudaError_t const errorCode);
  //! Wrap the cuMemAlloc call to include error checking
  static void *MemAlloc(unsigned const nbytes);
  //! Wrap cuMemcpyHtoD() to include error checking.
  static void MemcpyHtoD(void *ptr, void const *loc, unsigned nbytes);
  //! Wrap cuMemcpyDtoH() to include error checking.
  static void MemcpyDtoH(void *loc, void *ptr, unsigned nbytes);
  //! Wrap cuMemFree() to include error checking.
  static void MemFree(void *ptr);

protected:
  // IMPLEMENTATION

private:
  // NESTED CLASSES AND TYPEDEFS

  // IMPLEMENTATION

  // DATA
  int deviceCount;
  std::vector<std::vector<size_t>> computeCapability;
  std::vector<std::string> deviceName;

  // Device Properties: 1st dim is the device number.
  std::vector<int> m_maxthreadsperblock;
  std::vector<std::array<int, 3>> m_maxthreadsdim;
  std::vector<std::array<int, 3>> m_maxgridsize;
  std::vector<size_t> m_sharedmemperblock;
  std::vector<int> m_totalconstantmemory;
  std::vector<int> m_simdwidth;
  std::vector<size_t> m_mempitch;
  std::vector<int> m_regsperblock;
  std::vector<int> m_clockrate;
  std::vector<int> m_texturealign;
};

} // end namespace rtt_device

#endif // device_GPU_Device_hh

//------------------------------------------------------------------------------------------------//
// end of device/GPU_Device.hh
//------------------------------------------------------------------------------------------------//
