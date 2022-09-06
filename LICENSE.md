# DRACO Software

## Terms and Conditions of Use

Please read these terms and conditions of use carefully before using this software. Version 7_X_X of
this code has been reviewed by the Los Alamos National Laboratory (LANL) Export Control and the
Classification Teams and deemed to be Open Source Software under a 3-Clause BSD license.  The
copyright disclosure record is C19028 for "Draco-7_X_X".

Version 6_X_X was also released as Open Source Software and given the copyright disclosure record
LA-CC-16-016.

## Copyright

Copyright (C) 2019. Triad National Security, LLC.  All rights reserved.

This software was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National
Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of
Energy/National Nuclear Security Administration.

All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of
Energy/National Nuclear Security Administration. The U.S. Government is granted for itself and
others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material
to reproduce, prepare derivative works, distribute copies to the public, perform publicly and
display publicly, and to permit others to do so.

This is open source software; you can redistribute it and/or modify it under the terms of the
3-Clause BSD License. If software is modified to produce derivative works, such modified software
should be clearly marked, so as not to confuse it with the version available from LANL. Full text of
the 3-Clause BSD License can be found in the LICENSE file in the main development branch of the
repository.

## 3-Clause BSD License

This program is open source under the BSD-3 License.  Redistribution and use in source and binary
forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and
  the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice, this list of conditions
  and the following disclaimer in the documentation and/or other materials provided with the
  distribution.

- Neither the name of Triad National Security, LLC, Los Alamos National Laboratory, LANL, the
  U.S. Government, nor the names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY TRIAD NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL TRIAD NATIONAL SECURITY, LLC OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Additional License Information

The Draco software suite links against external libraries. Copyrights to these associated software
packages may are listed below:

1. Various MPI implementations ([OpenMPI](https://www.open-mpi.org/community/license.php), MPICH,
   Microsoft MPI, Intel MPI, etc.);

   An open source Message Passing Interface implementation. The Open MPI Project is an open source
   Message Passing Interface implementation that is developed and maintained by a consortium of
   academic, research, and industry partners. Open MPI is therefore able to combine the expertise,
   technologies, and resources from all across the High Performance Computing community in order to
   build the best MPI library available.  Open MPI offers advantages for system and software
   vendors, application developers and computer science researchers.

2. [Gnu Scientific Library](http://www.gnu.org/software/gsl/) - GSL is distributed under the terms
   of the GNU General Public License (GPL);

   The GNU Scientific Library (GSL) is a numerical library for C and C++ programmers. It is free
   software under the GNU General Public License.  The library provides a wide range of mathematical
   routines such as random number generators, special functions and least-squares fitting.  There
   are over 1000 functions in total with an extensive test suite.

3. [LAPACK](http://www.netlib.org/lapack/LICENSE.txt) or equivalent

   LAPACK version 3.X is a comprehensive FORTRAN library that does linear algebra operations
   including matrix inversions, least squared solutions to linear sets of equations, eigenvector
   analysis, singular value decomposition, etc. It is a very comprehensive and reputable package
   that has found extensive use in the scientific community.

4. [Random123](https://www.deshawresearch.com/downloads/download_random123.cgi/Random123_License.txt)

   Random123 is a library of 'counter-based' random number generators (CBRNGs), in which the Nth
   random number can be obtained by applying a stateless mixing function to N instead of the
   conventional approach of using N iterations of a stateful transformation.

5. [CUDA](http://docs.nvidia.com/cuda/eula/index.html#nvidia-cuda-toolkit-license-agreement)

6. [Metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)

   METIS is a set of serial programs for partitioning graphs, partitioning finite element meshes,
   and producing fill reducing orderings for sparse matrices. The algorithms implemented in METIS
   are based on the multilevel recursive-bisection, multilevel k-way, and multi-constraint
   partitioning schemes.

7. experimental/mdspan - [link](https://github.com/kokkos/mdspan) Kokkos v. 2.0, Copyright (2019)
   Sandia Corporation, BSD-3 License Git version hash: a7990884f090365787a90cdc12e689822d642c65

8. [ParMetis](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview)

   ParMETIS is an MPI-based parallel library that implements a variety of algorithms for
   partitioning unstructured graphs, meshes, and for computing fill-reducing orderings of sparse
   matrices.

9. `ds++/terminal.h` and `ds++/terminal_base.h` developed by Ondřej Čertík,
   [link](https://github.com/certik/terminal) MIT License
