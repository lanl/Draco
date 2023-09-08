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
  darwin | ccscs* | roci* | rz* | sn*) ;;
  *) die ".gitlab/ci/environments.sh :: SITE_ID not recognized, SITE_ID = ${SITE_ID}" ;;
esac

echo "    DRACO_ENV = ${DRACO_ENV}"

#------------------------------------------------------------------------------#
# Darwin
#------------------------------------------------------------------------------#
if [[ "${SITE_ID}" == "darwin" ]]; then
  [[ "${SLURM_JOB_PARTITION}x" == "x" ]] && DRACO_ARCH="x86_64" || \
    DRACO_ARCH="${SLURM_JOB_PARTITION}"
  export DRACO_ARCH
  run "module use --append /projects/draco/Modules"
  if ! [[ -f "/projects/draco/Modules/${DRACO_ENV}.lua" ]]; then
    # Look for TCL version of the module.
    if ! [[ -f "/projects/draco/Modules/${DRACO_ENV}" ]]; then
      die ".gitlab/ci/environments.sh :: DRACO_ENV not recognized, DRACO_ENV = ${DRACO_ENV}"
    fi
    # if [[ "${MPIARCH:-notset}" == "openmpi" ]]; then
    #   disable_openib=$(sinfo -N -n "${HOSTNAME}" -o %all | grep -c ib:none)
    #   if [[ ${disable_openib} != 0 ]]; then
    #     export MPI_PREFLAGS="--mca btl ^openib"
    #   fi
    # fi
  fi
  run "module load ${DRACO_ENV}"
  if [[ "${SLURM_JOB_PARTITION}" =~ "volta" || "${SLURM_JOB_PARTITION}" =~ "gpu" ]]; then
    module load cuda
  fi

#------------------------------------------------------------------------------#
# CCS-NET
#------------------------------------------------------------------------------#
elif [[ "${SITE_ID}" =~ "ccscs" ]]; then
  run "module use --append /ccs/codes/radtran/Modules"
  export PATH=/scratch/vendors/bin:$PATH # clang-format
  if ! [[ -f "/ccs/codes/radtran/Modules/draco/${DRACO_ENV}.lua" ]]; then
    die ".gitlab/ci/environments.sh :: DRACO_ENV not recognized, DRACO_ENV = ${DRACO_ENV}"
  fi
  run "module load draco/${DRACO_ENV}"

#------------------------------------------------------------------------------#
# RZ systems (LLNL)
#------------------------------------------------------------------------------#
elif [[ "${SITE_ID}" =~ "rzvernal" ]] || [[ "${SITE_ID}" =~ "rzansel" ]]; then
  run "module use --append /usr/gapps/jayenne/Modules/${SITE_ID}"
  for sm in draco lapse; do
    [[ $(module list 2>&1 | grep -c $sm) -gt 0 ]] && module unload $sm
  done
  case ${DRACO_ENV} in
    lapse* | draco* ) ;;
    *) die ".gitlab/ci/environments.sh :: DRACO_ENV not recognized, DRACO_ENV = ${DRACO_ENV}" ;;
  esac
  run "module load ${DRACO_ENV}"

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
# Rocinante
#------------------------------------------------------------------------------#
elif [[ "${SITE_ID}" =~ "rocinante" ]]; then
  run "module use --append /usr/projects/draco/Modules/ats3"
  case ${DRACO_ENV} in
    lapse* | draco* ) ;;
    *) die ".gitlab/ci/environments.sh :: DRACO_ENV not recognized, DRACO_ENV = ${DRACO_ENV}" ;;
  esac
  run "module load ${DRACO_ENV}"

fi

#--------------------------------------------------------------------------------------------------#
# Report the environment...
#--------------------------------------------------------------------------------------------------#

run "module list"

#--------------------------------------------------------------------------------------------------#
# End environments.sh
#--------------------------------------------------------------------------------------------------#
