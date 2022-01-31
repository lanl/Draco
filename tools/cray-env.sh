#!/bin/bash
#--------------------------------------------------------------------------------------------------#
# Cray (ARM) Environment setups (capulin, thunder)
#--------------------------------------------------------------------------------------------------#

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
if [[ $(df | grep yellow | grep -c lustre) -gt 0 ]]; then
  export siblings="capulin"
else
  export siblings="thunder"
fi

# The following toolchains will be used when releasing code
environments="cce11env"

# Extra cmake options
export CONFIG_BASE+=" -DCMAKE_VERBOSE_MAKEFILE=ON"

#--------------------------------------------------------------------------------------------------#
# Specify environments (modules)
#--------------------------------------------------------------------------------------------------#

if [[ "${ddir:-false}" == false ]] ;then
  echo "FATAL ERROR: Expected ddir to be set in the environment. (cray-env.sh)"
  exit 1
fi

case "$ddir" in

  #--------------------------------------------------------------------------------------------------#
  draco-7_10* | draco-7_11* | draco-7_12* | draco-7_13*)
    function cce11env()
    {
      unset partition
      unset jobnameext

      sysname=$(/usr/projects/hpcsoft/utilities/bin/sys_name)
      module use --append "/usr/projects/draco/Modules/${sysname}"

      run "module unload draco"
      run "module unload PrgEnv-intel"
      run "module load PrgEnv-cray"
      run "module load draco/cce11"
      run "module list"
      CC=$(which cc)
      CXX=$(which CC)
      FC=$(which ftn)
      # export LD_LIBRARY_PATH="$CRAY_LD_LIBRARY_PATH":"$LD_LIBRARY_PATH"
      export CC CXX FC
    }

    ;;

  *)
    die "cray-env.sh:: did not set any build environments, ddir = $ddir."
    ;;

  #--------------------------------------------------------------------------------------------------#
esac

#--------------------------------------------------------------------------------------------------#
# Sanity check
#--------------------------------------------------------------------------------------------------#

for env in $environments; do
  if [[ $(fn_exists "$env") -gt 0 ]]; then
    if [[ ${verbose:-false} != false ]]; then echo "export -f $env"; fi
    export -f "${env?}"
  else
    die "Requested environment $env is not defined."
  fi
done

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
