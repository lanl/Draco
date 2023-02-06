# Draco

## GitHub Status

[![GNU Build Status](https://github.com/lanl/Draco/actions/workflows/ci-gcc-mpi.yml/badge.svg)](https://github.com/lanl/Draco/actions/workflows/ci-gcc-mpi.yml)
[![LLVM Build Status](https://github.com/lanl/Draco/actions/workflows/ci-llvm-mpi.yml/badge.svg)](https://github.com/lanl/Draco/actions/workflows/ci-llvm-mpi.yml)
[![Windows Build Status](https://ci.appveyor.com/api/projects/status/yp8r9jxl2gc9n1fs/branch/develop?svg=true)](https://ci.appveyor.com/project/lanl/Draco)
[![codecov.io](https://codecov.io/github/lanl/Draco/coverage.svg?branch=develop)](https://codecov.io/github/lanl/Draco/branch/develop)
[![Latest Version](https://img.shields.io/github/release/lanl/draco.svg?style=flat-square)](https://github.com/lanl/Draco/releases)
[![PyPI](https://img.shields.io/pypi/l/Django.svg)](https://github.com/lanl/Draco/blob/develop/LICENSE.md)
[![GitHub Super-Linter](https://github.com/lanl/Draco/workflows/Lint%20Code%20Base/badge.svg)](https://github.com/marketplace/actions/super-linter)

## Summary

Draco is an object-oriented component library geared towards numerically intensive, radiation
(particle) transport applications built for parallel computing hardware.  It consists of
semi-independent packages and a robust build system.  The packages in draco provide a set of
components that can be used by multiple clients to build transport codes.  The build system can also
be extracted for use in clients.

To clone draco:

    git clone https://github.com/lanl/Draco.git

To get started, please see
[Development - Quick Start Guide](https://github.com/lanl/Draco/wiki/Development---Quick-Start)
in the wiki. Pull requests must satisfy the requirements listed in the
[Style Guide](https://github.com/lanl/Draco/wiki/Style-Guide).

Spack - Draco can be installed via the [spack package manager](https:://github.com/spack/spack)

## Synopsis of Active Draco Packages

* c4 - A communications library for message passing interfaces (MPI).
  * For builds without MPI, all communication commands will be no-op functions
    (`DRACO_C4={MPI;SCALAR}`). Also provides a wrapper to [libquo](https://github.com/lanl/libquo)
* cdi - Access to material data. The Common Data Interface (CDI) specifies a common abstraction for
  objects and libraries that return material data (opacities, atomic cross sections,
  equation-of-state data, etc.)
* cdi_analytic - Analytic models for physical data
* cdi_CPEloss - Stopping power models for charged particles.
* cdi_eospac - Equation-of-State data
  * These classes wrap the [EOSPAC6 libraries](https://laws.lanl.gov/projects/data/eos.html) that
    read sesame files; Commonly used to access gray opacity data and heat capacities.
    (`EOSPAC_LIB_DIR=<path>`).
* cdi_ipcress - Gray and multigroup opacities
  * The classes in this component will read and parse opacity values from an IPCRESS file produced by
    TOPS.
* cdi_ndi - Wrapper around LANL's Nuclear Data Interface (NDI) library for accessing [nuclear data
  tables](https://xweb.lanl.gov/projects/data/nuclear/ndi/ndi.html).
* compton_interface - Provides access to Compton scattering models and data as provided by the [CSK
  library](https://gitlab.lanl.gov/CSK/CSK)..
* compton_tools - Provides performant converter, reader, and interpolation tools for Compton data
* device - Wrapper for heterogeneous device communication
  * The classes in this component provide access to DaCS (deprecated) and CUDA calls for use on
    heterogeneous architecture platforms (GPU machines).
* diagnostics - CPP macros that are activated at compile time that can provide
  additional diagnostic verbosity during calculations.
* ds++ - Basic services and data structures library.
  * Array containers, assertion and Design-by-Contract, access to low level OS functions, file and
    path manipulation, unit test system, etc.
* experimental - Implementations of C++ features that are not yet available to all compilers that
    are used when installing Draco libraries.  Currently, this component provides
    [mdspan](https://github.com/kokkos/mdspan).
* fc4 - Fortran wrapper around portions of c4
* fit - Least squares fitting routines.
* fpe_trap - Catch IEEE floating point exceptions
* FortranCheck - Test Fortran compatibility and interoperability
  * The examples in this component will demonstrate if the Fortran compiler is working; if Fortran/C
    inter-language linking/running is working and sample `ISO_C_BINDING` calls.
* kde - KDE are a natural extension of quadrature. They reconstruct a smooth representation of a
  sampled field based on a specified Kernel shape. This populates the initial KDE infrastructure and
  some simple 1D reconstructions.
* linear - direct solvers for small linear systems of equations.
* mesh - Encapsulate mesh definition and accessors.
* mesh_element - defines fundamental mesh element types used by mesh, meshReaders and
  RTT_Format_Reader.
* meshReaders - Read RTT format mesh files.
* min - Optimization routines. Find the minimum of a function.
* norms - Calculate norms for data sets.
* ode - Ordinary differential equation solvers (e.g.: Runge-Kutta).
* [parser](https://github.com/lanl/Draco/tree/develop/src/parser) - Lexical file/input parser.
* quadrature - access to angular integration functions and related data.
* rng - A random number generator component
  * The primary set of functions provided by this component were derived from
    [Random123](https://www.deshawresearch.com/downloads/download_random123.cgi) random number
    library.  A few additional random number generators are also provided.
* roots - Root finding algorithms
* RTT_Format_Reader - meshReaders implementation for RTT format mesh files or input-streams.
* special_functions - Specialized math functions like factorial and Dirac delta.
* units - Provides encapsulated unit systems, functions to convert between unit systems and physical
  constants.
* VendorChecks - A testing component to ensure that discovered and used vendor libraries behave as
  expected.
* viz - Generates [Ensight](https://www.ansys.com/products/fluids/ansys-ensight) files for data
  visualization.

## Deprecated packages

* lapack_wrap - C++ wrapper for BLAS and LAPACK (removed 2022 July 12).
* plot2d - C++ interface to XMGrace 2-dimensional plotting (removed 2020 Dec 3).
* shared_lib - Dynamically load/unload shared object libraries via dl load (removed 2020 Sep 7).
* traits - A traits class used by viz.

## Authors

Many thanks go to Draco's [contributors](https://github.com/lanl/Draco/graphs/contributors).

Draco was originally written by staff from Los Alamos's [CCS-2 Computational Physics and Methods
Group](http://www.lanl.gov/org/padste/adtsc/computer-computational-statistical-sciences/computational-physics-methods/index.php):

> *CCS-2 Draco Team:* Kelly G. Thompson, Matt A. Cleveland, Ryan T. Wollaeger, Ben R. Ryan,
> Alex R. Long, Kendra P. Long, James S. Warsa, and Jae H. Chang.
>
> *Prior Contributors:* Kent. G. Budge, Gabriel M. Rockefeller, Allan B. Wollaber, Tim Kelley,
> Rob B. Lowrie, Andrew T. Till, Paul W. Talbot, Katherine J. Wang, Ondrej Certik, Peter Ahrens,
> David A. Dixon, Massimiliano Rosa, Todd J. Urbatsch, Daniel Holladay, Howard Pritchard,
> Jeff D. Densmore, Jeff Furnish, John McGhee, Kris C. Garrett, Mike Buksas, Nick Myers,
> Paul Henning, Randy Roberts, Seth Johnson, Todd Adams, Tom Evans, Lori Pritchett-Sheats,
> and Seth D. Cook.

## Release

Draco is released under the BSD 3-Clause License. For more details see the
[LICENSE file](https://github.com/lanl/Draco/blob/develop/LICENSE.md).
