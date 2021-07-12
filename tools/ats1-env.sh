#!/bin/bash
#--------------------------------------------------------------------------------------------------#
# Cray (ATS-1) Environment setups
#--------------------------------------------------------------------------------------------------#

# verbose=true
export KMP_AFFINITY=disabled
export VENDOR_DIR=/usr/projects/draco/vendors
export CRAYPE_LINK_TYPE=dynamic

# Sanity Check (Cray machines have very fragile module systems!)
if [[ -d "${ParMETIS_ROOT_DIR:-false}" ]]; then
  echo "ERROR: This script should be run from a clean environment."
  echo "       Try running 'rmdracoenv'."
  exit 1
fi

# symlinks will be generated for each machine that point to the correct
# installation directory.
if [[ $(df | grep yellnfs2 | grep -c jayenne) -gt 0 ]]; then
  export siblings="trinitite"
else
  export siblings="trinity"
fi

# The following toolchains will be used when releasing code
environments="cce11env intel1904env intel1904env-knl"

# Extra cmake options
export CONFIG_BASE+=" -DCMAKE_VERBOSE_MAKEFILE=ON"

# SLURM
avail_queues=$(sacctmgr -np list assoc user="$LOGNAME" | sed -e 's/.*|\(.*dev.*\|.*access.*\)|.*/\1/' | sed -e 's/|.*//')
case "$avail_queues" in
  *access*) access_queue="-A access --qos=access" ;;
  *dev*) access_queue="--qos=dev" ;;
esac
export access_queue

#--------------------------------------------------------------------------------------------------#
# Specify environments (modules)
#--------------------------------------------------------------------------------------------------#

if [[ "${ddir:-false}" == false ]] ;then
  echo "FATAL ERROR: Expected ddir to be set in the environment. (ats1-env.sh)"
  exit 1
fi

case "$ddir" in

  #------------------------------------------------------------------------------------------------#
  draco-7_9* | draco-7_10* | draco-7_11*)
    function cce11env()
    {
      unset partition
      unset jobnameext

      sysname=$(/usr/projects/hpcsoft/utilities/bin/sys_name)
      module use --append "/usr/projects/draco/Modules/${sysname}"

      run "module unload draco"
      run "module unload PrgEnv-intel"
      run "module load PrgEnv-cray"
      run "module load draco/cce1101"
      run "module list"
      CC=$(which cc)
      CXX=$(which CC)
      FC=$(which ftn)
      # export LD_LIBRARY_PATH="$CRAY_LD_LIBRARY_PATH":"$LD_LIBRARY_PATH"
      export CC CXX FC
    }

    function intel1904env()
    {
      unset partition
      unset jobnameext

      sysname=$(/usr/projects/hpcsoft/utilities/bin/sys_name)
      module use --append "/usr/projects/draco/Modules/${sysname}"

      run "module unload draco"
      run "module unload PrgEnv-cray"
      run "module load PrgEnv-intel"
      run "module load draco/intel19"
      run "module list"
      CC=$(which cc)
      CXX=$(which CC)
      FC=$(which ftn)
      # export LD_LIBRARY_PATH="$CRAY_LD_LIBRARY_PATH":"$LD_LIBRARY_PATH"
      export CC CXX FC
    }

    function intel1904env-knl()
    {
      export partition="-p knl --exclude=nid00192"
      export jobnameext="-knl"

      sysname=$(/usr/projects/hpcsoft/utilities/bin/sys_name)
      module use --append "/usr/projects/draco/Modules/${sysname}"

      run "module unload draco"
      run "module unload PrgEnv-cray"
      run "module load PrgEnv-intel"
      run "module load draco/intel19"
      if ! [[ "${CRAY_CPU_TARGET}" == mic-knl ]]; then
        run "module swap craype-haswell craype-mic-knl"
      fi
      run "module list"
      CC=$(which cc)
      CXX=$(which CC)
      FC=$(which ftn)
      export LD_LIBRARY_PATH="$CRAY_LD_LIBRARY_PATH":"$LD_LIBRARY_PATH"
      export CC CXX FC
    }
    ;;

  #------------------------------------------------------------------------------------------------#

  *)
    die "ats1-env.sh:: did not set any build environments, ddir = $ddir."
    ;;


  #------------------------------------------------------------------------------------------------#
esac

#--------------------------------------------------------------------------------------------------#
# Sanity check
#--------------------------------------------------------------------------------------------------#

for env in $environments; do
  if [[ $(fn_exists "$env") -gt 0 ]]; then
    if [[ ${verbose:-false} != false ]]; then echo "export -f $env"; fi
    # exporting these environment is required on CTS-1, but breaks the module system on trinitite.
    # Ugh.
    #
    # export -f $env
  else
    die "Requested environment $env is not defined."
  fi
done

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
