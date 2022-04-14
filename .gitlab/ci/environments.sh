#!/bin/bash -l
# -*- Mode: sh -*-
#--------------------------------------------------------------------------------------------------#
# File  : .gitlab/ci/environments.sh
# Date  : Monday, Jun 01, 2020, 15:43 pm
# Author: Kelly Thompson
# Note  : Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#

echo "==> Setting up CI environment..."
echo "    SITE_ID = ${SITE_ID}"
[[ -n "${SCHEDULER_PARAMETERS}" ]] && echo "    Using: salloc ${SCHEDULER_PARAMETERS}"

case ${SITE_ID} in
  darwin | ccscs* | sn* | trinitite*) ;;
  *) die ".gitlab/ci/environments.sh :: SITE_ID not recognized, SITE_ID = ${SITE_ID}" ;;
esac

echo "    DRACO_ENV = ${DRACO_ENV}"
echo "    ARCH     = ${ARCH}"

#------------------------------------------------------------------------------#
# Darwin
#------------------------------------------------------------------------------#
if [[ "${SITE_ID}" == "darwin" ]]; then
  DRACO_ARCH=$(/usr/projects/draco/vendors/bin/target_arch)
  run "module use --append /projects/draco/Modules"
  case ${DRACO_ENV} in
    arm-gcc930 | power9-gcc930-smpi | power9-xl16117 | x64-gcc930* | x64-intel1905) ;;
    *) die ".gitlab/ci/environments.sh :: DRACO_ENV not recognized, DRACO_ENV = ${DRACO_ENV}" ;;
  esac
  if [[ "${MPIARCH:-notset}" == "openmpi" ]]; then
    disable_openib=$(sinfo -N -n "${HOSTNAME}" -o %all | grep -c ib:none)
    if [[ ${disable_openib} != 0 ]]; then
      export MPI_PREFLAGS="--mca btl ^openib"
    fi
  fi
  export DRACO_ARCH
  run "module load draco/${DRACO_ENV}"
  if [[ "${SLURM_JOB_PARTITION}" =~ "volta" || "${SLURM_JOB_PARTITION}" =~ "gpu" ]]; then
    module load cuda/11.0
  fi

#------------------------------------------------------------------------------#
# CCS-NET
#------------------------------------------------------------------------------#
elif [[ "${SITE_ID}" =~ "ccscs" ]]; then
  run "module use --append /ccs/codes/radtran/Modules"
  export PATH=/scratch/vendors/bin:$PATH # clang-format
  if ! [[ -d "/ccs/codes/radtran/Modules/draco/${DRACO_ENV}" ]]; then
    die ".gitlab/ci/environments.sh :: DRACO_ENV not recognized, DRACO_ENV = ${DRACO_ENV}" ;;
  fi
  run "module load draco/${DRACO_ENV}"

#------------------------------------------------------------------------------#
# Snow
#------------------------------------------------------------------------------#
elif [[ "${SITE_ID}" =~ "snow" ]]; then
  run "module use --append /usr/projects/draco/Modules/cts1"
  # export PATH=/scratch/vendors/bin:$PATH # clang-format
  case ${DRACO_ENV} in
    lapse* | draco* ) ;;
    *) die ".gitlab/ci/environments.sh :: DRACO_ENV not recognized, DRACO_ENV = ${DRACO_ENV}" ;;
  esac
  run "module load ${DRACO_ENV}"

#------------------------------------------------------------------------------#
# Trinitite
#------------------------------------------------------------------------------#
elif [[ "${SITE_ID}" =~ "trinitite" ]]; then
  run "module use --append /usr/projects/draco/Modules/trinitite"
  # export PATH=/scratch/vendors/bin:$PATH # clang-format
  case ${DRACO_ENV} in
    lapse* | draco* ) ;;
    *) die ".gitlab/ci/environments.sh :: DRACO_ENV not recognized, DRACO_ENV = ${DRACO_ENV}" ;;
  esac
  run "module unload draco lapse"
  run "module load ${DRACO_ENV}"

fi

#--------------------------------------------------------------------------------------------------#
# Report the environment...
#--------------------------------------------------------------------------------------------------#

run "module list"

#--------------------------------------------------------------------------------------------------#
# End environments.sh
#--------------------------------------------------------------------------------------------------#
