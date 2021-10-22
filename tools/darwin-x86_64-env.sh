#!/bin/bash
#--------------------------------------------------------------------------------------------------#
# Darwin Environment setups (x86_64 + gpu)
#--------------------------------------------------------------------------------------------------#

# shellcheck source=tools/darwin-env.sh
source "${draco_script_dir:-unknown}/darwin-env.sh"

# symlinks will be generated for each machine that point to the correct installation directory.
export siblings="darwin-x86_64"

# The following toolchains will be used when releasing code:
environments="x86gcc930env x86intel1905env"

# Special setup for CTS-1: replace the 'draco-latest' symlink
pushd "/usr/projects/${package:-notset}" || exit
if [[ -L draco-latest ]]; then rm draco-latest; fi
if [[ -d "${source_prefix:=notset}" ]]; then ln -s "${source_prefix}" draco-latest; fi
popd || exit

#--------------------------------------------------------------------------------------------------#
# Specify environments (modules)
#--------------------------------------------------------------------------------------------------#

case "${ddir:=notset}" in

  #---------------------------------------------------------------------------#
  draco-7_10* | draco-7_11* | draco-7_12*)
    function x86gcc930env()
    {
      export darwin_queue="-p volta-x86"
      run "module purge"
      run "module use --append /projects/draco/Modules"
      run "module load draco/x64-gcc930"
      if [[ ${SLURM_JOB_PARTITION} =~ "volta-" || ${SLURM_JOB_PARTITION} =~ "gpu" ]]; then
        run "module load cuda/11.0"
      fi
      run "module list"

      CXX=$(which g++)
      CC=$(which gcc)
      FC=$(which gfortran)
      MPIEXEC_EXECUTABLE=$(which mpirun)
      unset MPI_ROOT
      export CXX
      export CC
      export FC
      export MPIEXEC_EXECUTABLE
    }

    function x86intel1905env()
    {
      export darwin_queue="-p volta-x86"
      run "module purge"
      run "module use --append /projects/draco/Modules"
      run "module load draco/x64-intel1905"
      #if [[ ${SLURM_JOB_PARTITION} =~ "volta-" || ${SLURM_JOB_PARTITION} =~ "gpu" ]]; then
      #  run "module load cuda/11.0"
      #fi
      run "module list"
      MPIEXEC_EXECUTABLE=$(which mpirun)
      export MPIEXEC_EXECUTABLE
      unset MPI_ROOT
    }
    ;;

  *) die "darwin-x86_64-env.sh:: did not set any build environments, ddir = $ddir." ;;
esac

#--------------------------------------------------------------------------------------------------#
# Sanity check
#--------------------------------------------------------------------------------------------------#

for env in $environments; do
  if [[ $(fn_exists "$env") -gt 0 ]]; then
    [[ "${verbose:-notset}" != "notset" ]] && echo "export -f $env"
    export -f "${env?}"
  else
    die "Requested environment $env is not defined."
  fi
done

#--------------------------------------------------------------------------------------------------#
# End darwin-x86_64-env.sh
#--------------------------------------------------------------------------------------------------#
