#!/bin/bash -l
# ------------------------------------------------------------------------------------------------ #
# File  : ./.gitlab/ci/gitlab-nightly-regress.sh
# Date  : Tuesday, Jun 02, 2020, 10:31 am
# Author: Kelly Thompson <kgt@lanl.gov>
# Note  : Copyright (C) 2021-2023 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

print_use()
{
  echo " "
  echo "Usage: ${0#*/} -h"
  echo "               -a [EXTRA_CMAKE_ARGS]"
  echo "               -c Use current environment."
  echo "               -e Do not source the environments.sh file used by nightly regressions."
  echo "               -g [GENERATOR]"
  echo "               -i [DRACO_BINARY_DIR]"
  echo "               -m [Configure,Build,Test,Submit]"
  echo "               -n [numproc]"
  echo "               -p [PROJECT]"
  echo "               -s [DRACO_SOURCE_DIR]"
  echo "               -t [CTEST_BUILD_TYPE]"
  echo "               -v VERBOSITY = {-V, -VV}"
  echo "               -x [MAXLOAD]"
  echo "               -z cdash upload token"
  echo " "
}

# defaults:
modes="Configure,Build,Test,Submit"
export VERBOSITY="-V"

while getopts "a:eg:hi:m:n:p:s:t:v:x:z:" opt; do
case $opt in
a)  EXTRA_CMAKE_ARGS="$OPTARG"
    export EXTRA_CMAKE_ARGS ;;
e)  source_environment=0 ;;
g)  CTEST_CMAKE_GENERATOR="${OPTARG}"
    export CTEST_CMAKE_GENERATOR ;;
h)  print_use; exit 0 ;;
i)  DRACO_BINARY_DIR="$OPTARG"
    export DRACO_BINARY_DIR ;;
m)  modes="$OPTARG" ;;
n)  CTEST_NPROC="$OPTARG"
    export CTEST_NPROC;;
p)  PROJECT="$OPTARG"
    export PROJECT ;;
s)  DRACO_SOURCE_DIR="$OPTARG"
    export DRACO_SOURCE_DIR ;;
t)  CMAKE_BUILD_TYPE="$OPTARG"
    export CMAKE_BUILD_TYPE ;;
v)  VERBOSITY="$OPTARG"
    export VERBOSITY ;;
x)  MAXLOAD="$OPTARG"
    export MAXLOAD ;;
z)  cdash_upload_auth_token="$OPTARG"
    export cdash_upload_auth_token ;;
\?) echo "" ;echo "invalid option: -$OPTARG"; print_use; exit 1 ;;
:)  echo "" ;echo "option -$OPTARG requires an argument."; print_use; exit 1 ;;
esac
done

[[ -n $CI_PROJECT_DIR ]] || CI_PROJECT_DIR=${DRACO_SOURCE_DIR}
[[ -n $CTEST_MODE ]] || CTEST_MODE="Experimental"
[[ -n $CTEST_SITE ]] || CTEST_SITE=$(uname -n)
CTEST_SITE=${CTEST_SITE//.lanl.gov/}
if [[ "${CTEST_BUILD_NAME}x" == "x" ]]; then
  # shellcheck disable=2001
  compiler=$(echo "$CXX" | sed -e 's%.*/%%')
  CTEST_BUILD_NAME="$CMAKE_BUILD_TYPE-$compiler"
fi
export CI_PROJECT_DIR CTEST_MODE CTEST_SITE CTEST_BUILD_NAME

echo "==> Setting modes to $modes."

# preliminaries and environment
set -e
# shellcheck source=.gitlab/ci/common.sh
source "${DRACO_SOURCE_DIR}/.gitlab/ci/common.sh"

if [[ $source_environment == "0" ]]; then
  # do nothing
  echo "Using current build environment..."
else
  # shellcheck source=.gitlab/ci/environments.sh
  source "${DRACO_SOURCE_DIR}/.gitlab/ci/environments.sh"
fi

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
[[ -n $MAKEFILE_FLAGS ]] || export MAKEFILE_FLAGS="-j $MAXLOAD"
export CTEST_NPROC MAXLOAD CODECOV MEMCHECK_CONFIGURATION

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

echo "To rerun manually, cd to CI directory (pwd), set these variables, run ctest."
echo " "
run "pwd"
echo " "
echo "export ARCH=${ARCH}"
echo "export AUTODOCDIR=${AUTODOCDIR}"
echo "export BUILD_FLAGS=${BUILD_FLAGS}"
echo "export CI_PROJECT_DIR=${CI_PROJECT_DIR}"
echo "export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
echo "export CTEST_BUILD_NAME=${CTEST_BUILD_NAME}"
echo "export CTEST_MODE=${CTEST_MODE}"
echo "export CTEST_NPROC=${CTEST_NPROC}"
echo "export EXTRA_CMAKE_ARGS=${EXTRA_CMAKE_ARGS}"
echo "export DRACO_BINARY_DIR=${DRACO_BINARY_DIR}"
echo "export DRACO_SOURCE_DIR=${DRACO_SOURCE_DIR}"
echo "export PROJECT=${PROJECT}"
echo "export SITE_ID=${SITE_ID}"
echo "export TEST_EXCLUSIONS=${TEST_EXCLUSIONS}"
echo " "
echo "export CMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAINFILE}"
echo "export CODECOV=${CODECOV}"
echo "export MEMCHECK_COMMAND_OPTIONS=${MEMCHECK_COMMAND_OPTIONS}"
echo "export MEMCHECK_CONFIGURATION=${MEMCHECK_CONFIGURATION}"
echo " "
echo "modes       = ${modes}"
echo " "

if [[ "${cdash_upload_auth_token}x" == "x" ]]; then
  # shellcheck disable=2154
  die "environment variable cdash_upload_auth_token was not set. Aborting"
fi

run "ctest ${VERBOSITY} -S ${DRACO_SOURCE_DIR}/.gitlab/ci/draco-nightly.cmake,${modes}"

[[ "${AUTODOC}" == "ON" ]] && cp "${DRACO_SOURCE_DIR}/.gitlab/ci/index.html" "${AUTODOCDIR}/."

echo -e "\n======== end .gitlab-nightly-regress.sh ==========\n"

#--------------------------------------------------------------------------------------------------#
# End .gitlab/ci/gitlab-nightly-regress.sh
#--------------------------------------------------------------------------------------------------#
