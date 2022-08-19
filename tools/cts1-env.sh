#!/bin/bash
#--------------------------------------------------------------------------------------------------#
# CTS-1 Environment setups
#--------------------------------------------------------------------------------------------------#

export VENDOR_DIR=/usr/projects/draco/vendors

# symlinks will be generated for each machine that point to the correct
# installation directory.
if [[ $(df | grep yellow | grep -c jayenne) -gt 0 ]]; then
  export siblings="snow badger kodiak"
else
  export siblings="fire ice cyclone"
fi

# The following toolchains will be used when releasing code
environments="lapse20intelenv lapse20gccenv eapintelenv"

# Extra cmake options
export CONFIG_BASE+=" -DCMAKE_VERBOSE_MAKEFILE=ON"

# SLURM
avail_queues=$(sacctmgr -np list assoc user="$LOGNAME" | sed -e 's/.*|\(.*dev.*\|.*access.*\)|.*/\1/' | sed -e 's/|.*//')
case "$avail_queues" in
  *access*) access_queue="-A access --qos=access" ;;
  *dev*) access_queue="--qos=dev" ;;
esac
export access_queue

# Special setup for CTS-1: replace the 'latest' symlink
if [[ ${package:-false} == false ]] ; then die "package not defined"; fi
if [[ ${source_prefix:-false} == false ]] ; then die "source_prefix not defined"; fi
(cd "/usr/projects/${package:=draco}" || exit; if [[ -L latest ]]; then rm latest; fi; ln -s "${source_prefix:-source_prefix}" latest)

#--------------------------------------------------------------------------------------------------#
# Specify environments (modules)
#--------------------------------------------------------------------------------------------------#

if [[ ${ddir:=false} == false ]] ;then
  echo "FATAL ERROR: Expected ddir to be set in the environment. (cts1-env.sh)"
  exit 1
fi

case "${ddir}" in

  #--------------------------------------------------------------------------------------------------#
  draco-7_14*)
    function lapse20intelenv()
    {
      # intel-19.0.4 + openmpi-3.1.6
      run "module purge"
      run "module use --append /usr/projects/draco/Modules/cts1"
      run "module load lapse/2.0-intel"
      run "module list"
    }
    function lapse20gccenv()
    {
      # gcc-9.3.0 + openmpi-3.1.6
      run "module purge"
      run "module use --append /usr/projects/draco/Modules/cts1"
      run "module load lapse/2.0-gnu"
      run "module list"
    }
    function eapintelenv()
    {
      # intel-19.0.4 + openmpi-4.1.1
      run "module purge"
      run "module use --append /usr/projects/draco/Modules/cts1"
      run "module load eapse/beta-intel-ompi4"
      run "module list"
    }
    ;;

  #--------------------------------------------------------------------------------------------------#
  draco-7_9* | draco-7_10* | draco-7_11* | draco-7_12* | draco-7_13*)
    function intel1904env
    {
      run "module purge"
      run "module use --append /usr/projects/draco/Modules/cts1"
      run "module load draco/intel19"
      run "module list"
    }
    function gcc930env()
    {
      run "module purge"
      run "module use --append /usr/projects/draco/Modules/cts1"
      run "module load draco/gcc9"
      run "module list"
    }
    function lapse18intelenv()
    {
      run "module purge"
      run "module use --append /usr/projects/draco/Modules/cts1"
      run "module load lapse/1.8-intel"
      run "module list"
    }
    ;;

#--------------------------------------------------------------------------------------------------#

  *)
    die "cts1-env.sh:: did not set any build environments, ddir = $ddir."
    ;;

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

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
