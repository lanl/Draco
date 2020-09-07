#!/bin/bash
#--------------------------------------------------------------------------------------------------#
# ATS-2 Environment setups
#--------------------------------------------------------------------------------------------------#

export VENDOR_DIR=/usr/gapps/jayenne/vendors

# Sanity Check (Cray machines have very fragile module systems!)
if [[ -d $ParMETIS_ROOT_DIR ]]; then
  echo "ERROR: This script should be run from a clean environment."
  echo "       Try running 'rmdracoenv'."
  exit 1
fi

# symlinks will be generated for each machine that point to the correct
# installation directory.
if [[ `df | grep -c rz_gapps` -gt 0 ]]; then
  export siblings="rzansel"
else
  export siblings="sierra"
fi

# The following toolchains will be used when releasing code
environments="gcc831env xl20200318env"

# Extra cmake options
export CONFIG_BASE+=" -DCMAKE_VERBOSE_MAKEFILE=ON"

# job launch options
case $siblings in
  rzansel) job_launch_options="-q pdebug" ;;
  sierra) job_launch_options="-q pbatch" ;;
  *) die "*FATAL ERROR* (ats2-env.sh) I only know how to set job_launch_options for rzansel and sierra." ;;
esac
export job_launch_options

# Special setup for ATS-1: replace the 'latest' symlink
(cd /usr/gapps/jayenne; if [[ -L draco-latest ]]; then rm draco-latest; fi; ln -s $source_prefix draco-latest)

#--------------------------------------------------------------------------------------------------#
# Specify environments (modules)
#--------------------------------------------------------------------------------------------------#

if ! [[ $ddir ]] ;then
  echo "FATAL ERROR: Expected ddir to be set in the environment. (ats2-env.sh)"
  exit 1
fi

case $ddir in

  #--------------------------------------------------------------------------------------------------#
  draco-7_7*)
    function gcc831env()
    {
      run "module purge"
      run "module use /usr/gapps/user_contrib/spack.20200402/share/spack/lmod/linux-rhel7-ppc64le/Core"
      run "module load StdEnv"
      run "module unload cuda spectrum-mpi xl"
      unset CMAKE_PREFIX_PATH
      unset CPATH
      unset LD_LIBRARY_PATH
      unset LIBRARY_PATH
      run "module load gcc/8.3.1 spectrum-mpi/2019.06.24 cuda"
      run "module load python/3.7.2 cmake/3.17.0 git gsl numdiff random123 metis netlib-lapack"
      run "module load parmetis superlu-dist trilinos csk"
      run "module load eospac/6.4.0 libquo/1.3.1 ndi"
      run "module list"
      run "module avail"
      unset MPI_ROOT
      CXX=`which g++`
      CC=`which gcc`
      FC=`which gfortran`
    }
    function xl20200318env()
    {
      run "module purge"
      run "module use /usr/gapps/user_contrib/spack.20200402/share/spack/lmod/linux-rhel7-ppc64le/Core"
      run "module load StdEnv"
      run "module unload cuda spectrum-mpi xl"
      unset CMAKE_PREFIX_PATH
      unset CPATH
      unset LD_LIBRARY_PATH
      unset LIBRARY_PATH
      run "module load xl/2020.03.18 spectrum-mpi/2019.06.24"
      run "module load python/3.7.2 cmake/3.17.0 git gsl numdiff random123 metis netlib-lapack"
      run "module load parmetis superlu-dist trilinos csk"
      run "module load eospac/6.4.0 libquo/1.3.1 ndi"
      # trilinos possible non-spack solution at:
      # /usr/gapps/jayenne/vendors/trilinos-12.18.1/xl-2019.12.23-spectrum-mpi-2019.06.24
      # export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/usr/gapps/jayenne/vendors/trilinos-12.18.1/xl-2019.12.23-spectrum-mpi-2019.06.24
      run "module list"
      run "module avail"
      unset MPI_ROOT
    }
    ;;

  #--------------------------------------------------------------------------------------------------#
  draco-7_4* | draco-7_5* | draco-7_6*)
    function gcc731env()
    {
      run "module purge"
      run "module use /usr/gapps/jayenne/vendors-ec/spack.20190616/share/spack/lmod/linux-rhel7-ppc64le/Core"
      run "module load StdEnv"
      run "module unload cuda spectrum-mpi xl"
      unset CMAKE_PREFIX_PATH
      unset CPATH
      unset LD_LIBRARY_PATH
      unset LIBRARY_PATH
      run "module load gcc/7.3.1 spectrum-mpi/2019.06.24 cuda"
      run "module load python/3.7.2 cmake/3.14.5 git gsl numdiff random123 metis netlib-lapack"
      run "module load parmetis superlu-dist trilinos csk"
      run "module load eospac/6.4.0 libquo/1.3.1 ndi"
      run "module list"
      run "module avail"
      unset MPI_ROOT
      CXX=`which g++`
      CC=`which gcc`
      FC=`which gfortran`
    }
    function xl20191223()
    {
      run "module purge"
      run "module use /usr/gapps/jayenne/vendors-ec/spack.20190616/share/spack/lmod/linux-rhel7-ppc64le/Core"
      run "module load StdEnv"
      run "module unload cuda spectrum-mpi xl"
      unset CMAKE_PREFIX_PATH
      unset CPATH
      unset LD_LIBRARY_PATH
      unset LIBRARY_PATH
      run "module load xl/2019.12.23 spectrum-mpi/2019.06.24 cuda"
      run "module load python/3.7.2 cmake/3.14.5 git gsl numdiff random123 metis netlib-lapack"
      run "module load parmetis superlu-dist trilinos csk"
      run "module load eospac/6.4.0 libquo/1.3.1 ndi"
      # trilinos possible non-spack solution at:
      # /usr/gapps/jayenne/vendors/trilinos-12.18.1/xl-2019.12.23-spectrum-mpi-2019.06.24
      # export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/usr/gapps/jayenne/vendors/trilinos-12.18.1/xl-2019.12.23-spectrum-mpi-2019.06.24
      run "module list"
      run "module avail"
      unset MPI_ROOT
    }
    ;;

  *)
    die "ats2-env.sh:: did not set any build environments, ddir = $ddir."
    ;;

  #--------------------------------------------------------------------------------------------------#
esac

#--------------------------------------------------------------------------------------------------#
# Sanity check
#--------------------------------------------------------------------------------------------------#

for env in $environments; do
  if [[ `fn_exists $env` -gt 0 ]]; then
    if [[ $verbose ]]; then echo "export -f $env"; fi
    export -f $env
  else
    die "Requested environment $env is not defined."
  fi
done


##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
