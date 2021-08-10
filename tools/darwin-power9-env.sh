#!/bin/bash
#--------------------------------------------------------------------------------------------------#
# Darwin Environment setups (Power9)
#--------------------------------------------------------------------------------------------------#

# shellcheck source=tools/darwin-env.sh
source "${draco_script_dir:-unknown}/darwin-env.sh"

# symlinks will be generated for each machine that point to the correct installation directory.
export siblings="darwin-power9"

# The following toolchains will be used when releasing code:
environments="p9gcc930env p9xl16117env"

#--------------------------------------------------------------------------------------------------#
# Specify environments (modules)
#--------------------------------------------------------------------------------------------------#

case "${ddir}" in

  #---------------------------------------------------------------------------#
  draco-7_10* | draco-7_11*)
    function p9gcc930env()
    {
      export darwin_queue="-p power9-asc -A asc-priority"
      run "module purge"
      run "module use --append /projects/draco/Modules"
      run "module load draco/power9-gcc930-smpi"
      run "module list"

      # work around for known openmpi issues: https://re-git.lanl.gov/draco/draco/-/issues/938
      # eliminates warnings: "there are more than one active ports on host"

      # export OMPI_MCA_btl=^openib
      UCX_NET_DEVICES=mlx5_0:1
      UCX_WARN_UNUSED_ENV_VARS=n
      OMPI_MCA_pml=ob1
      OMPI_MCA_btl=self,vader

      CC=$(which gcc)
      FC=$(which gfortran)
      MPIEXEC_EXECUTABLE=$(which mpirun)
      unset MPI_ROOT
      export CXX
      export CC
      export FC
      export MPIEXEC_EXECUTABLE
      export UCX_NET_DEVICES
      export UCX_WARN_UNUSED_ENV_VARS
      export OMPI_MCA_pml
      export OMPI_MCA_btl
    }

    function p9xl16117env()
    {
      export darwin_queue="-p power9-asc -A asc-priority"
      run "module purge"
      echo "VENDOR_DIR = $VENDOR_DIR"
      echo "DRACO_ARCH = $DRACO_ARCH"
      run "module use --append /projects/draco/Modules"
      run "module load draco/power9-xl16117-smpi"
      run "module list"

      # work around for known openmpi issues: https://re-git.lanl.gov/draco/draco/-/issues/938
      # warnings: "there are more than one active ports on host"

      # export OMPI_MCA_btl=^openib
      UCX_NET_DEVICES=mlx5_0:1
      UCX_WARN_UNUSED_ENV_VARS=n
      OMPI_MCA_pml=ob1
      OMPI_MCA_btl=self,vader
      MPIEXEC_EXECUTABLE=$(which mpirun)
      unset MPI_ROOT
      export UCX_NET_DEVICES
      export UCX_WARN_UNUSED_ENV_VARS
      export OMPI_MCA_pml
      export OMPI_MCA_btl
      export MPIEXEC_EXECUTABLE
    }
    ;;

  *) die "darwin-power9-env.sh:: did not set any build environments, ddir = $ddir." ;;
esac

#--------------------------------------------------------------------------------------------------#
# Sanity check
#--------------------------------------------------------------------------------------------------#

for env in ${environments}; do
  if [[ $(fn_exists "$env") -gt 0 ]]; then
    if [[ "${verbose:-false}" != false ]]; then echo "export -f $env"; fi
    export -f "${env?}"
  else
    die "Requested environment $env is not defined."
  fi
done

#--------------------------------------------------------------------------------------------------#
# End darwin-power9-env.sh
#--------------------------------------------------------------------------------------------------#
