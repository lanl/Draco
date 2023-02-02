#!/bin/bash -l
# ------------------------------------------------------------------------------------------------ #
# File  : ./.gitlab/ci/gitlab-nightly-regress.sh
# Date  : Tuesday, Jun 02, 2020, 10:31 am
# Author: Kelly Thompson <kgt@lanl.gov>
# Note  : Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

print_use()
{
  echo " "
  echo "Usage: ${0#*/} -h -m [Configure,Build,Test,Submit]"
  echo " "
}

# defaults:
modes="Configure,Build,Test,Submit"

while getopts "hm:" opt; do
case $opt in
h)  print_use; exit 0 ;;
m)  modes="$OPTARG" ;;
\?) echo "" ;echo "invalid option: -$OPTARG"; print_use; exit 1 ;;
:)  echo "" ;echo "option -$OPTARG requires an argument."; print_use; exit 1 ;;
esac
done

echo "==> Setting modes to OPTARG = $OPTARG."

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
[[ "${EXTRA_CMAKE_ARGS}" =~ "CODE_COVERAGE" ]] && CODECOV=ON || CODECOV=OFF
if [[ -z "${MEMCHECK_CONFIGURATION}" ]] && [[ "${EXTRA_CMAKE_ARGS}" =~ "ENABLE_MEMORYCHECK" ]]; then
  MEMCHECK_CONFIGURATION=ON
fi
export CODECOV MEMCHECK_CONFIGURATION MAXLOAD

#echo -e "\n========== printenv ==========\n"
#[[ -z "${SLURM_NODELIST}" ]] || echo "SLURM_NODELIST = ${SLURM_NODELIST}"
#echo "HOSTNAME       = ${HOSTNAME}"
# echo -e "NPROC       = ${NPROC}\n"
# echo -e "CTEST_NPROC = ${CTEST_NPROC}\n"
# echo -e "MAXLOAD     = ${MAXLOAD}\n"
# echo -e "EXTRA_CMAKE_ARGS = ${EXTRA_CMAKE_ARGS}\n"
# echo -e "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}\n"
# echo -e "EXTRA_CTEST_ARGS = ${EXTRA_CTEST_ARGS}\n"
# echo -e "CODECOV          = ${CODECOV}\n"
# echo -e "MEMCHECK_CONFIGURATION = ${MEMCHECK_CONFIGURATION}\n"
# if [[ "${AUTODOC}" == "ON" ]]; then
#   run "which doxygen"
#   run "which latex"
#   run "which pdflatex"
#   echo "TEXINPUTS = ${TEXINPUTS}"
# fi
# if [[ "${EXTRA_CTEST_ARGS}" =~ memcheck ]]; then
#   run "which valgrind"
# fi
printenv >& environment.log
# run "pwd"

#--------------------------------------------------------------------------------------------------#
# Setup compiler flags
#--------------------------------------------------------------------------------------------------#

for i in C CXX Fortran CUDA; do

  # Turn build warnings into fatal errors
  if [[ "${DRACO_ENV}" =~ "xl16" ]] || [[ "${DRACO_ENV}" =~ "power9-xl" ]]; then
    eval export ${i}_FLAGS+=\" -qhalt=w\"
  else
    if [[ "${i}" == "Fortran" ]] && [[ "${FC}" =~ "ifort" ]]; then
      # do not add the flag
      echo "Skip adding -Werror to ${i}_FLAGS (FC = ${FC})"
    elif [[ "${i}" == "CUDA" ]]; then
      eval export ${i}_FLAGS+=\" -Werror all-warnings\"
    else
      # Works for gcc, llvm, and icpc (not ifort).
      eval export ${i}_FLAGS+=\" -Werror\"
    fi
  fi

  # Enable coverage for Debug builds
  if [[ ${CMAKE_BUILD_TYPE} == Debug ]] && [[ "${CODECOV}" == "ON" ]] ; then
    eval export ${i}_FLAGS+=\" --coverage\"
  fi
done
echo "C_FLAGS       = $C_FLAGS"
echo "CXX_FLAGS     = $CXX_FLAGS"
# shellcheck disable=2154
echo "Fortran_FLAGS = $Fortran_FLAGS"

# Tweak MPI settings
# export OMPI_MCA_btl=self,sm
export OMPI_MCA_btl=^openib

export CTEST_NPROC

#--------------------------------------------------------------------------------------------------#
# Build and run the tests for draco; post results to CDash.
#--------------------------------------------------------------------------------------------------#

# echo "To rerun manually, cd to CI directory (pwd), set these variables, run ctest."
# echo " "
# run "pwd"
# echo " "
# echo "ARCH               = ${ARCH}"
# echo "AUTODOCDIR         = ${AUTODOCDIR}"
# echo "BUILD_FLAGS        = ${BUILD_FLAGS}"
# echo "CI_PROJECT_DIR     = ${CI_PROJECT_DIR}"
# echo "CMAKE_BUILD_TYPE   = ${CMAKE_BUILD_TYPE}"
# echo "CTEST_BUILD_NAME   = ${CTEST_BUILD_NAME}"
# echo "CTEST_MODE         = ${CTEST_MODE}"
# echo "CTEST_NPROC        = ${CTEST_NPROC}"
# echo "EXTRA_CMAKE_ARGS   = ${EXTRA_CMAKE_ARGS}"
# echo "DRACO_BINARY_DIR   = ${DRACO_BINARY_DIR}"
# echo "DRACO_SOURCE_DIR   = ${DRACO_SOURCE_DIR}"
# echo "PROJECT            = ${PROJECT}"
# echo "SITE_ID            = ${SITE_ID}"
# echo "TEST_EXCLUSIONS    = ${TEST_EXCLUSIONS}"
# echo " "
# echo "CMAKE_TOOLCHAIN_FILE     = ${CMAKE_TOOLCHAINFILE}"
# echo "CODECOV                  = ${CODECOV}"
# #echo "MEMCHECK_COMMAND_OPTIONS = ${MEMCHECK_COMMAND_OPTIONS}"
# echo "MEMCHECK_CONFIGURATION   = ${MEMCHECK_CONFIGURATION}"
# #echo "MEMORYCHECK_TYPE         = ${MEMORYCHECK_TYPE}"
# echo " "
# echo "modes = ${modes}"
echo " "

run "ctest -V -S ${DRACO_SOURCE_DIR}/.gitlab/ci/draco-nightly.cmake,${modes}"

[[ "${AUTODOC}" == "ON" ]] && cp "${DRACO_SOURCE_DIR}/.gitlab/ci/index.html" "${AUTODOCDIR}/."

echo -e "\n======== end .gitlab-nightly-regress.sh ==========\n"

#--------------------------------------------------------------------------------------------------#
# End .gitlab/ci/gitlab-nightly-regress.sh
#--------------------------------------------------------------------------------------------------#
