#!/bin/bash
#-------------------------------------------------------------------------------------------------#
# Darwin Environment setups (ARM)
#-------------------------------------------------------------------------------------------------#

# shellcheck source=tools/darwin-env.sh
source "${draco_script_dir:-unknown}/darwin-env.sh"

# symlinks will be generated for each machine that point to the correct installation directory.
export siblings="darwin-arm"

# The following toolchains will be used when releasing code:
environments="armgcc930env"

#--------------------------------------------------------------------------------------------------#
# Specify environments (modules)
#--------------------------------------------------------------------------------------------------#

case "${ddir}" in

  #---------------------------------------------------------------------------#
  draco-7_10* | draco-7_11*)
    function armgcc930env()
    {
      export darwin_queue="-p arm"
      run "module purge"
      echo "VENDOR_DIR = ${VENDOR_DIR}"
      echo "DRACO_ARCH = ${DRACO_ARCH}"
      run "module use --append /projects/draco/Modules"
      run "module load draco/arm-gcc930"
      run "module list"

      CXX=$(which g++)
      CC=$(which gcc)
      FC=$(which gfortran)
      MPIEXEC_EXECUTABLE=$(which mpirun)
      unset MPI_ROOT
      # work around for known openmpi issues: https://rtt.lanl.gov/redmine/issues/1229
      OMPI_MCA_btl=^openib
      UCX_NET_DEVICES=mlx5_0:1
      export CXX
      export CC
      export FC
      export MPIEXEC_EXECUTABLE
      export OMPI_MCA_btl
      export UCX_NET_DEVICES
    }
    ;;

  *)
    die "darwin-arm-env.sh:: did not set any build environments, ddir = $ddir."
    ;;

  #---------------------------------------------------------------------------#

esac

#--------------------------------------------------------------------------------------------------#
# Sanity check
#--------------------------------------------------------------------------------------------------#

for env in $environments; do
  if [[ $(fn_exists $env) -gt 0 ]]; then
    if [[ "${verbose:-false}" != "false" ]]; then "export -f $env"; fi
    export -f "${env?}"
  else
    die "Requested environment $env is not defined."
  fi
done

#--------------------------------------------------------------------------------------------------#
# End darwin-arm-env.sh
#--------------------------------------------------------------------------------------------------#
