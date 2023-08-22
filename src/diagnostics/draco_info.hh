//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/draco_info.hh
 * \author Kelly Thompson
 * \date   Wednesday, Nov 07, 2012, 18:49 pm
 * \brief  Small executable that prints the version and copyright strings.
 * \note   Copyright (C) 2014-2023 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_diagnostics_draco_info_hh
#define rtt_diagnostics_draco_info_hh

#include "ds++/config.h"
#include <string>
#include <vector>

namespace rtt_diagnostics {

//================================================================================================//
/*!
 * \class DracoInfo
 * \brief Store and present basic information about the current draco build.
 *
 * The constructed string will take this form:
 *
 * \verbatim
Draco-7_11_20210915, build date 2021 Sep 15, build type: Debug, DBC: 7, DRACO_DIAGNOSTICS: 0

CCS-2 Draco Team: Kelly G. Thompson, Kent G. Budge, Matt A. Cleveland, Ryan T. Wollaeger,
     Ben R. Ryan, Alex R. Long, Kendra P. Long, James S. Warsa, Jae H. Chang, Andrew T. Till,
     and David A. Dixon.

Prior Contributors: Gabriel M. Rockefeller, Allan B. Wollaber, Tim Kelley, Rob B. Lowrie,
     Paul W. Talbot, Katherine J. Wang, Seth D. Cook, Ondrej Certik, Peter Ahrens,
     Massimiliano Rosa, Todd J. Urbatsch, Daniel Holladay, Jeff D. Densmore, Howard Pritchard,
     Jeff Furnish, John McGhee, Kris C. Garrett, Mike Buksas, Nick Myers, Paul Henning,
     Randy Roberts, Seth Johnson, Todd Adams, Tom Evans, and Lori Pritchett-Sheats.

For information, send e-mail to draco@lanl.gov.

Build information:
    Build type        : DEBUG
    Library type      : Shared
    System type       : Linux
    Site name         : ccscs3
    CUDA support      : disabled
    MPI support       : enabled
      mpirun cmd      : /ccs/opt/.../bin/mpiexec -n <N> -bind-to none
    OpenMP support    : enabled
    Design-by-Contract: 7, features = Insist, Require, Check, Ensure
    Diagnostics       : 0
    Diagnostics Timing: disabled
    CXX Compiler      : /scratch/vendors/spack.20180425/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/
                        gcc-8.1.0-3c5hjkqndywdp3w2l5vts62xlllrsbtq/bin/g++
    CXX_FLAGS         : -Wcast-align -Wpointer-arith -Wall -pedantic -Wno-expansion-to-defined
                        -Wnarrowing -march=native -fopenmp -Werror
    C Compiler        : /scratch/vendors/spack.20180425/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/
                        gcc-8.1.0-3c5hjkqndywdp3w2l5vts62xlllrsbtq/bin/gcc
    C_FLAGS           : -Wcast-align -Wpointer-arith -Wall -pedantic -Wno-expansion-to-defined
                        -Wnarrowing -march=native -fopenmp -Werror
    Fortran Compiler  : /scratch/vendors/spack.20180425/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/
                        gcc-8.1.0-3c5hjkqndywdp3w2l5vts62xlllrsbtq/bin/gfortran
    Fortran_FLAGS     : -ffree-line-length-none -cpp -march=native -fopenmp

 * \endverbatim
 */
//================================================================================================//
class DracoInfo {
public:
  // IMPLELEMENTATION
  // ================

  // Constructors
  // ------------
  DracoInfo();

  // Actions
  // -------

  /*! \brief Construct an information message that includes Draco's version, copyright and basic
   *         build parameters. */
  std::string fullReport() const;

  //! Version and (C)
  std::string briefReport() const;

  //! Version only
  std::string versionReport() const;

private:
  // DATA
  // ----

  std::string const release;
  std::string const copyright;
  std::string const contact;
  std::string const build_type;
  std::string library_type;
  std::string system_type;
  std::string site_name;
  bool gpu = {false};
  bool mpi = {false};
  std::string mpirun_cmd;
  bool openmp = {false};
  std::string diagnostics_level;
  bool diagnostics_timing = {false};
  std::string cxx;
  std::string cxx_flags;
  std::string cc;
  std::string cc_flags;
  std::string fc;
  std::string fc_flags;
  std::string gpu_compiler; //! < full path to the cuda compiler (if any)
  std::string gpu_flags;    //!< Compiler flags for cuda
};

} // end namespace rtt_diagnostics

#endif // rtt_diagnostics_draco_info_hh

//------------------------------------------------------------------------------------------------//
// end of draco_info.hh
//------------------------------------------------------------------------------------------------//
