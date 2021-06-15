#!/bin/bash -l
#--------------------------------------------------------------------------------------------------#
# File  : ./travis-run-tests.sh
# Date  : Tuesday, Jan 17, 2017, 15:55 pm
# Author: Kelly Thompson
# Note  : Copyright (C) 2017-2020, Triad National Security, LLC., All rights are reserved.
#
# .travis.yml calls this script to build draco and run the tests.
#--------------------------------------------------------------------------------------------------#

# preliminaries and environment
set -e

cd "${SOURCE_DIR:-/home/travis/Draco}"
source tools/common.sh

#--------------------------------------------------------------------------------------------------#
# Style Check
#--------------------------------------------------------------------------------------------------#
if [[ ${STYLE} ]]; then
  echo "checking style conformance..."

  run "git branch -a"

  # Return integer > 0 if 'develop' branch is found.
  function find_dev_branch
  {
    set -f
    git branch -a | grep -v origin | grep -c develop
    set +f
  }

  # Ensure the 'develop' branch is available.  In some cases (merge a branch that lives at
  # github.com/lanl), the develop branch is missing in the Travis checkout. Since we only test files
  # that are modified when compared to the 'develop' branch, the develop branch must be available
  # locally.
  num_dev_branches_found=$(find_dev_branch)
  if [[ $num_dev_branches_found == 0 ]]; then
    echo "no develop branches found."
    # Register the develop branch in draco/.git/config
    run "git config --local remote.origin.fetch +refs/heads/develop:refs/remotes/origin/develop"
    # Download the meta-data for the 'develop' branch
    run "git fetch"
    # Create a local tracking branch
    run "git branch -t develop origin/develop"
  fi

  # clang-format is installed at /usr/bin.
  export PATH=$PATH:/usr/bin
  # extract the TPL list from the Dockerfile
  CLANG_FORMAT_VER="$(grep CLANG_FORMAT_VER tools/Dockerfile-style-checks.txt | head -n 1 | sed -e 's/.*=//')"
  export CLANG_FORMAT_VER
  tools/check_style.sh -t

#--------------------------------------------------------------------------------------------------#
# 1. GCC based build and test, or
# 2. Autodoc
#--------------------------------------------------------------------------------------------------#
elif [[ "${COMPILER:=GCC}" == "GCC" ]]; then
  if [[ "${AUTODOC}" == "ON" ]]; then
    echo "checking autodoc compliance..."
  else
    echo "checking build and test (GCC)..."
  fi

  # extract the TPL list from the Dockerfile
  DRACO_TPL="$(grep DRACO_TPL tools/Dockerfile-spack-gcc.txt | sed -e 's/.*=//' | sed -e 's/\"//g')"
  DRACO_GCC_TPL="$(grep DRACO_GCC_TPL tools/Dockerfile-spack-gcc.txt | head -n 1 | sed -e 's/.*=//' | sed -e 's/\"//g')"
  DRACO_DOC_TPL="$(grep DRACO_DOC_TPL tools/Dockerfile-spack-gcc.txt | head -n 1 | sed -e 's/.*=//' | sed -e 's/\"//g')"
  export DRACO_TPL DRACO_GCC_TPL DRACO_DOC_TPL

  # Environment setup for the build...
  dmodules="${DRACO_TPL//@//} ${DRACO_GCC_TPL}"
  if [[ "${AUTODOC}" == "ON" ]]; then
    dmodules+=" ${DRACO_DOC_TPL}"
  fi
  for item in ${dmodules}; do
    run "module load $item"
  done
  run "module list"
  run "module avail"

  # Provide a newer lcov that is compatible with gcc-8.
  #run "cp -r ${SOURCE_DIR}/tools/spack/lcov ${SPACK_ROOT}/var/spack/repos/builtin/packages/."
  #run "spack install lcov@1.14"
  #run "spack load lcov"

  echo "GCCVER = ${GCCVER}"
  if [[ "${GCCVER:=notset}" == "notset" ]]; then
    CXX=$(which g++)
    CC=$(which gcc)
    FC=$(which gfortran)
    GCOV=$(which gcov)
  else
    CXX=$(which "g++-${GCCVER}")
    CC=$(which "gcc-${GCCVER}")
    FC=$(which "gfortran-${GCCVER}")
    GCOV=$(which "gcov-${GCCVER}")
  fi
  export CXX CC FC GCOV
  echo "CXX    = ${CXX}"
  echo "FC     = ${FC}"
  echo "GCOV   = ${GCOV}"

  # Extra options for mpirun. Travis requires these extra mpirun options
  export MPIEXEC_PREFLAGS="--allow-run-as-root --oversubscribe"

  export OMP_NUM_THREADS=2
  if [[ ${WERROR} ]]; then
    for i in C CXX Fortran; do
      eval export ${i}_FLAGS+=\" -Werror\"
    done
    eval export CUDA_FLAGS+=\" -Werror all-warnings\"
  fi
  if [[ ${COVERAGE:-OFF} == "ON" ]]; then
    CMAKE_OPTS="-DCODE_COVERAGE=ON"
  fi

  echo -e "\n========== printenv =========="
  #printenv
  printenv | grep _FLAGS
  echo " "

  if [[ ${BUILD_DIR:-notset} == "notset" ]]; then die "BUILD_DIR not set by environment."; fi
  run "mkdir -p ${BUILD_DIR}"
  run "cd ${BUILD_DIR}"

  echo " "
  if [[ -f CMakeCache.txt ]]; then
    echo "===== CMakeCache.txt ====="
    run "cat CMakeCache.txt"
  fi
  echo "========"
  run "cmake -DDRACO_C4=${DRACO_C4} ${CMAKE_OPTS} ${SOURCE_DIR}"
  echo "========"
  if [[ "${AUTODOC}" == "ON" ]]; then
    run "make autodoc"
    echo "========"
  else
    run "make -j 2"
    echo "========"
    # tstOMP_2 needs too many ppr (threads * cores) for Travis.
    run "ctest -j 2 -E \(c4_tstOMP_2\|c4_tstTermination_Detector_2\) --output-on-failure"
  fi
  if [[ ${COVERAGE} == "ON" ]]; then
    echo "========"
    run "make VERBOSE=1 covrep"
    # Uploading report to CodeCov
    bash <(curl -s https://codecov.io/bash) -f coverage.info || \
      echo "Codecov did not collect coverage reports"
  fi
  cd -

#--------------------------------------------------------------------------------------------------#
# 1. LLVM based build and test, or
# 2. clang-tidy
#--------------------------------------------------------------------------------------------------#
elif [[ "${COMPILER}" == "LLVM" ]]; then
  echo "checking build and test (LLVM)..."

  # extract the TPL list from the Dockerfile
  DRACO_TPL="$(grep DRACO_TPL tools/Dockerfile-spack-llvm.txt | sed -e 's/.*=//' | sed -e 's/\"//g')"
  DRACO_LLVM_TPL="$(grep DRACO_LLVM_TPL tools/Dockerfile-spack-llvm.txt | head -n 1 | sed -e 's/.*=//' | sed -e 's/\"//g')"
  export DRACO_TPL DRACO_LLVM_TPL

  # Environment setup for the build...
  dmodules="llvm ${DRACO_TPL//@//} ${DRACO_LLVM_TPL}"
  for item in ${dmodules}; do
    run "module load $item"
  done
  run "module list"
  run "module avail"

  CXX=$(which clang++)
  CC=$(which clang)
  FC=$(which gfortran)
  export CXX CC FC
  echo "CXX    = ${CXX}"
  echo "FC     = ${FC}"
  echo "GCOV   = ${GCOV}"

  # Extra options for mpirun. Travis requires these extra mpirun options
  export MPIEXEC_PREFLAGS="--allow-run-as-root --oversubscribe"

  export OMP_NUM_THREADS=2
  if [[ ${WERROR} ]]; then
    for i in C CXX Fortran; do
      eval export ${i}_FLAGS+=\" -Werror\"
    done
    eval export CUDA_FLAGS+=\" -Werror all-warnings\"
  fi

  echo -e "\n========== printenv =========="
  #printenv
  printenv | grep _FLAGS
  echo " "

  if [[ ${BUILD_DIR:-notset} == "notset" ]]; then die "BUILD_DIR not set by environment."; fi
  run "mkdir -p ${BUILD_DIR}"
  run "cd ${BUILD_DIR}"

  echo " "
  if [[ -f CMakeCache.txt ]]; then
    echo "===== CMakeCache.txt ====="
    run "cat CMakeCache.txt"
  fi
  if [[ "${CLANGTIDY:=OFF}" == "ON" ]]; then
    echo "==> enable clang-tidy CI mode"
    # -DCI_CLANG_TIDY=ON
    CMAKE_OPTS+=" -DDRACO_STATIC_ANALYZER=clang-tidy "
    echo "==> CMAKE_OPTS = ${CMAKE_OPTS}"
  fi
  echo "========"
  run "cmake -DDRACO_C4=${DRACO_C4} ${CMAKE_OPTS} ${SOURCE_DIR}"
  echo "========"
  run "make -j 2"
  echo "========"
  # tstOMP_2 needs too many ppr (threads * cores) for Travis.
  run "ctest -j 2 -E \(c4_tstOMP_2\|c4_tstTermination_Detector_2\) --output-on-failure"
  cd -
fi

# Finish up and report
echo "======== end .travis-run-tests.sh =========="

#--------------------------------------------------------------------------------------------------#
# End .travis-run-tests.sh
#--------------------------------------------------------------------------------------------------#
