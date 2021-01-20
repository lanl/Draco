#!/bin/bash
#-----------------------------------------------------------------------------#
# Darwin Environment setups (Power9)
#-----------------------------------------------------------------------------#

source $draco_script_dir/darwin-env.sh

# symlinks will be generated for each machine that point to the correct
# installation directory.
export siblings="darwin-power9"

# The following toolchains will be used when releasing code:
environments="p9gcc730env p9xl16117env"

#--------------------------------------------------------------------------------------------------#
# Specify environments (modules)
#--------------------------------------------------------------------------------------------------#

case $ddir in

  #---------------------------------------------------------------------------#
  draco-7_7* | draco-7_8*)
    function p9gcc730env()
    {
      export darwin_queue="-p power9-asc -A asc-priority"
      run "module purge"
      module use --append ${VENDOR_DIR}/user_contrib
      module use --append ${VENDOR_DIR}-ec/Modules/$DRACO_ARCH
      module load user_contrib

      cflavor="gcc-7.3.0"
      mflavor="$cflavor-openmpi-3.1.3"
      lflavor="lapack-3.8.0"
      noflavor="git gcc/7.3.0 cuda/11.0 cmake/3.17.3 random123 numdiff"
      compflavor="gsl/2.5-$cflavor netlib-lapack/3.8.0-$cflavor metis/5.1.0-$cflavor eospac/6.4.0-$cflavor openmpi/p9/3.1.3-gcc_7.3.0"
      mpiflavor="parmetis/4.0.3-$mflavor superlu-dist/5.4.0-${mflavor}-lapack-3.8.0 trilinos/12.14.1-${mflavor}-lapack-3.8.0 libquo/1.3-$mflavor caliper/2.0.1-$mflavor"
      ec_mf="csk/0.5.0-$cflavor ndi"
      # work around for known openmpi issues:
      # https://rtt.lanl.gov/redmine/issues/1229
      # eliminates warnings: "there are more than one active ports on host"
      # export OMPI_MCA_btl=^openib
      export UCX_NET_DEVICES=mlx5_0:1
      export UCX_WARN_UNUSED_ENV_VARS=n
      export OMPI_MCA_pml=ob1
      export OMPI_MCA_btl=self,vader

      export dracomodules="$noflavor $compflavor $mpiflavor $ec_mf"

      for m in $dracomodules; do
        module load $m
      done
      export CXX=`which g++`
      export CC=`which gcc`
      export FC=`which gfortran`
      export MPIEXEC_EXECUTABLE=`which mpirun`
      unset MPI_ROOT
      run "module list"
    }

    function p9xl16117env()
    {
      export darwin_queue="-p power9-asc -A asc-priority"
      run "module list"
      run "module purge"
      echo "VENDOR_DIR = $VENDOR_DIR"
      echo "DRACO_ARCH = $DRACO_ARCH"
      module use --append ${VENDOR_DIR}/user_contrib
      module use --append ${VENDOR_DIR}-ec/Modules/$DRACO_ARCH
      module load user_contrib

      cflavor="xl-16.1.1.7"
      mflavor="$cflavor-spectrum-mpi-10.3.0.1"
      noflavor="git ibm/xlc-16.1.1.7-xlf-16.1.1.7-gcc-7.4.0-cuda-11.0 cmake/3.17.3 random123 numdiff"
      compflavor="gsl/2.5-$cflavor netlib-lapack/3.8.0-$cflavor metis/5.1.0-$cflavor eospac/6.4.0-$cflavor ibm/smpi/p9/10.3.0.1-xlc-16.1.1.7-xlf-16.1.1.7-gcc-7.4.0-cuda-11.0"
      mpiflavor="parmetis/4.0.3-$mflavor libquo/1.3.1-$mflavor superlu-dist/5.4.0-${mflavor}-lapack-3.8.0 trilinos/12.14.1-${mflavor}-lapack-3.8.0 caliper/2.0.1-$mflavor"
      ec_mf="csk/0.5.0-$cflavor"
      # ndi

      # work around for known openmpi issues:
      # https://rtt.lanl.gov/redmine/issues/1229
      # eliminates warnings: "there are more than one active ports on host"
      # export OMPI_MCA_btl=^openib
      export UCX_NET_DEVICES=mlx5_0:1
      export UCX_WARN_UNUSED_ENV_VARS=n
      export OMPI_MCA_pml=ob1
      export OMPI_MCA_btl=self,vader

      export dracomodules="$noflavor $compflavor $mpiflavor $ec_mf"

      for m in $dracomodules; do
        module load $m
      done
      # export CXX=`which g++`
      # export CC=`which gcc`
      # export FC=`which gfortran`
      export MPIEXEC_EXECUTABLE=`which mpirun`
      unset MPI_ROOT
      run "module list"
    }
    ;;

  #---------------------------------------------------------------------------#
  draco-7_2* | draco-7_3* | draco-7_4* | draco-7_5* | draco-7_6* )
    function p9gcc730env()
    {
      export darwin_queue="-p power9-asc -A asc-priority"
      run "module purge"
      module use --append ${VENDOR_DIR}/user_contrib
      module use --append ${VENDOR_DIR}-ec/Modules/$DRACO_ARCH
      module load user_contrib

      cflavor="gcc-7.3.0"
      mflavor="$cflavor-openmpi-3.1.3"
      lflavor="lapack-3.8.0"
      noflavor="git gcc/7.3.0 cuda/10.1"
      compflavor="eospac/6.4.0-$cflavor cmake/3.17.0-$cflavor random123 numdiff/5.9.0-$cflavor gsl/2.5-$cflavor netlib-lapack/3.8.0-$cflavor metis/5.1.0-$cflavor openmpi/p9/3.1.3-gcc_7.3.0"
      mpiflavor="libquo/1.3-$mflavor parmetis/4.0.3-$mflavor superlu-dist/5.2.2-${mflavor}-$lflavor trilinos/12.14.1-cuda-10.1-${mflavor}-$lflavor"

      ec_mf="csk/0.5.0-$cflavor ndi"

      # work around for known openmpi issues:
      # https://rtt.lanl.gov/redmine/issues/1229
      # eliminates warnings: "there are more than one active ports on host"
      # export OMPI_MCA_btl=^openib
      export UCX_NET_DEVICES=mlx5_0:1
      export UCX_WARN_UNUSED_ENV_VARS=n
      export OMPI_MCA_pml=ob1
      export OMPI_MCA_btl=self,vader

      export dracomodules="$noflavor $compflavor $mpiflavor $ec_mf"

      for m in $dracomodules; do
        module load $m
      done
      export CXX=`which g++`
      export CC=`which gcc`
      export FC=`which gfortran`
      export MPIEXEC_EXECUTABLE=`which mpirun`
      unset MPI_ROOT
      run "module list"
    }
    ;;

  *)
    die "darwin-power9-env.sh:: did not set any build environments, ddir = $ddir."
    ;;


  #---------------------------------------------------------------------------#

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
