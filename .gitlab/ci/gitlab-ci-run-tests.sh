#!/bin/bash -l
#--------------------------------------------------------------------------------------------------#
# File  : ./.gitlab/ci/gitlab-ci-run-tests.sh
# Date  : Tuesday, Jun 02, 2020, 12:28 pm
# Author: Kelly Thompson <kgt@lanl.gov>
# Note  : Copyright (C) 2021-2026 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#

# Staging directory should be protected. Deployed artifacts will be set with chmod at the end of
# this script.
umask 77
umask
if [[ "${DEPLOY_RELEASE}" == "TRUE" ]]; then
  chmod_params="g+rX,o+rX" # release directory
else
  chmod_params="g+rwX,o+rX" # deploy directory
fi

# preliminaries and environment
set -e
# shellcheck source=.gitlab/ci/common.sh
source "${DRACO_SOURCE_DIR}/.gitlab/ci/common.sh"
# shellcheck source=.gitlab/ci/environments.sh
source "${DRACO_SOURCE_DIR}/.gitlab/ci/environments.sh"

if [[ $(which lscpu | grep -c lscpu) -gt 0 ]]; then
  # run "lscpu"
  # NPROC=`lscpu | grep CPU\(s\) | head -n 1 | awk '{ print $2 }'`
  NPROC=$(grep -c processor < /proc/cpuinfo)
fi

[[ -z "${CTEST_NPROC}" ]] && CTEST_NPROC=$NPROC || echo "limiting CTEST_NPROC = $CTEST_NPROC"
[[ -z "${MAXLOAD}" ]] && MAXLOAD=$NPROC || echo "limiting MAXLOAD = $MAXLOAD"
[[ ${EXTRA_CMAKE_ARGS} =~ "Ninja" ]] && MAKE_COMMAND="ninja" || MAKE_COMMAND="make -j -l ${MAXLOAD}"
export MAKE_COMMAND
[[ "${DEPLOY}" == "TRUE" ]] && EXTRA_CMAKE_ARGS="-DBUILD_TESTING=NO ${EXTRA_CMAKE_ARGS}"
if [[ -n "${DRACO_INSTALL_DIR}" ]]; then
  EXTRA_CMAKE_ARGS=" ${EXTRA_CMAKE_ARGS} -DCMAKE_INSTALL_PREFIX=${DRACO_INSTALL_DIR}"
fi

if [[ $(groups | grep -c ccsrad) == 1 ]]; then
  deploy_group=ccsrad
elif [[ $(groups | grep -c dacodes) == 1 ]]; then
  deploy_group=dacodes
else
  die "Neither ccsrad or dacodes can be used as Unix group for deployment."
fi

#if [[ -d /ccs/opt/texlive/2018/texmf-dist/tex/latex/newunicodechar ]]; then
#  TEXINPUTS="${TEXINPUTS:+${TEXINPUTS}:}/ccs/opt/texlive/2018/texmf-dist/tex/latex/newunicodechar"
#fi

# Report number of cores/threads, etc.
# run "lstopo --no-io --no-bridges --no-icaches"

echo -e "\n========== printenv ==========\n"
echo "HOSTNAME       = ${HOSTNAME}"
echo -e "SITE_ID     = ${SITE_ID}\n"
[[ -z "${SLURM_JOBID}" ]] || echo "SLURM_NODELIST = ${SLURM_JOBID}"
[[ -z "${SLURM_NODELIST}" ]] || echo "SLURM_NODELIST = ${SLURM_NODELIST}"
echo -e "NPROC       = ${NPROC}\n"
echo -e "CTEST_NPROC = ${CTEST_NPROC}\n"
echo -e "MAXLOAD     = ${MAXLOAD}\n"
echo -e "DEPLOY      = ${DEPLOY}\n"
echo -e "EXTRA_CMAKE_ARGS = ${EXTRA_CMAKE_ARGS}\n"
echo -e "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}\n"
echo -e "EXTRA_CTEST_ARGS = ${EXTRA_CTEST_ARGS}\n"
echo -e "MAKE_COMMAND      = ${MAKE_COMMAND}\n"
if [[ "${AUTODOC}" == "ON" ]]; then
  run "which doxygen"
  run "which latex"
  run "which pdflatex"
  echo "TEXINPUTS = ${TEXINPUTS}"
fi
if [[ "${EXTRA_CTEST_ARGS}" =~ memcheck ]]; then
  run "which valgrind"
fi
printenv >& environment.log
run "pwd"

[[ "${EXTRA_CMAKE_ARGS}" =~ "CODE_COVERAGE" ]] && CODECOV=ON || CODECOV=OFF

#--------------------------------------------------------------------------------------------------#
# Setup compiler flags
#--------------------------------------------------------------------------------------------------#

if [[ "${DEPLOY_RELEASE}" != "TRUE" ]]; then
  for i in C CXX Fortran CUDA; do

    # Turn build warnings into fatal errors for non-release-deploys
    if [[ "${i}" == "Fortran" ]]; then
      if [[ "${FC}" =~ "ifort" ]] || [[ "${FC}" =~ "ifx" ]] || \
            [[ "${INTEL_COMPILER_TYPE}" == "ONEAPI" ]] || \
            [[ "${PE_ENV}" == "INTEL" ]] ; then
        # do not add the flag
        echo "Skip adding -Werror to ${i}_FLAGS (FC = ${FC})"
      elif [[ "${FC}" =~ "crayftn" ]]; then
        eval export ${i}_FLAGS+=\" -h error_on_warning\"
      else
        eval export ${i}_FLAGS+=\" -Werror\"
      fi
    elif [[ "${i}" == "CUDA" ]]; then
      eval export ${i}_FLAGS+=\" -Werror all-warnings\"
    else
      # Works for gcc, llvm, and icpc (not ifort or ifx).
      eval export ${i}_FLAGS+=\" -Werror\"
    fi

    # Enable coverage for Debug builds
    if [[ ${CMAKE_BUILD_TYPE} == Debug ]] && [[ "${CODECOV}" == "ON" ]] ; then
      eval export ${i}_FLAGS+=\" --coverage\"
    fi
  done
fi
echo "C_FLAGS       = $C_FLAGS"
echo "CXX_FLAGS     = $CXX_FLAGS"
# shellcheck disable=2154
echo "Fortran_FLAGS = $Fortran_FLAGS"
echo "CUDA_FLAGS = $CUDA_FLAGS"

# Tweak MPI settings
# export OMPI_MCA_btl=self,sm
export OMPI_MCA_btl=^openib

#--------------------------------------------------------------------------------------------------#
# Draco
#--------------------------------------------------------------------------------------------------#

echo -e "\n========== Draco ==========\n"
if ! [[ -d ${DRACO_BINARY_DIR} ]]; then run "mkdir -p ${DRACO_BINARY_DIR}"; fi
run "cd ${DRACO_BINARY_DIR}"
set > ci-env.log
echo -e "\n========"
run "cmake ${EXTRA_CMAKE_ARGS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} ${DRACO_SOURCE_DIR}"
echo " "
echo -e "\n========\n"
if [[ "${AUTODOC}" == "ON" ]]; then
  # run "make -j -l $MAXLOAD autodoc"
  run "cmake --build . --target autodoc -- -j $MAXLOAD -l $MAXLOAD"
  # Run one test so that the CDash report isn't empty.
  run "ctest -R draco_release_tIMCInfo_a_d"
else
  # run "make -j -l $MAXLOAD" #  VERBOSE=1
  run "cmake --build . --target all -- -j $MAXLOAD -l $MAXLOAD"
  echo -e "\n========\n"
  if [[ "${TEST_EXCLUSIONS:-no}" != "no"  ]]; then
    EXTRA_CTEST_ARGS+=" -E ${TEST_EXCLUSIONS}"
  fi
  if [[ "${DEPLOY}" == "TRUE" ]]; then
    # Remove existing deployment directory (but only if 'make' w/o install works)
    # - This ensures removed files are also removed from the deployment.
    # - This ensures permissions are consistent and correct.
    if [[ -n "${DRACO_INSTALL_DIR}" ]]; then
      run "rm -rf ${DRACO_INSTALL_DIR}"
      run "mkdir -p ${DRACO_INSTALL_DIR}"
      run "chgrp ${deploy_group} ${DRACO_INSTALL_DIR}" || echo "Unable to chgrp"
      run "chmod ${chmod_params} ${DRACO_INSTALL_DIR}" || echo "Unable to chmod"
      run "chmod g+s ${DRACO_INSTALL_DIR}" || echo "Unable to chmod"
    fi
    # run "make -j -l $MAXLOAD install"
    run "cmake --build . --target install --  -j $MAXLOAD -l $MAXLOAD"

    # Double check permissions for installed (deployed) files
    run "chgrp -R ${deploy_group} ${DRACO_INSTALL_DIR}" || echo "Unable to chgrp"
    run "chmod -R ${chmod_params} ${DRACO_INSTALL_DIR}" || echo "Unable to chmod"

  else
    EXTRA_CTEST_ARGS+=" --output-on-failure --repeat until-pass:2"
    run "ctest ${EXTRA_CTEST_ARGS} -j ${CTEST_NPROC} --test-load ${MAXLOAD}"
  fi
fi

# Generate a coverage report (text and html)
if [[ "${CODECOV}" == "ON" && "${DEPLOY}" != "TRUE" ]]; then
  # run "make covrep"
  run "cmake --build . --target covrep"
fi

if [[ ${CTEST_MODE} == "Nightly" ]]; then
  echo -e "\n========\n"
  echo "Nighlty regressions should run a ctest in script mode and post to CDash"
fi

echo -e "\n======== end .gitlab-ci-run-tests.sh ==========\n"

#--------------------------------------------------------------------------------------------------#
# End .gitlab/ci/gitlab-ci-run-tests.sh
#--------------------------------------------------------------------------------------------------#
