#!/bin/bash -l
# -*- Mode: sh -*-
#--------------------------------------------------------------------------------------------------#
# File  : regression/sripts/common.sh
# Date  : Tuesday, May 31, 2016, 14:48 pm
# Author: Kelly Thompson
# Note  : Copyright (C) 2016-2020, Triad National Security, LLC., All rights are reserved.
#
# Summary: Misc bash functions useful during development of code.
#--------------------------------------------------------------------------------------------------#

shopt -s extglob # Allow variable as case condition

#--------------------------------------------------------------------------------------------------#
# Helpful functions
#--------------------------------------------------------------------------------------------------#

function dracohelp()
{
  echo -e "Bash functions defined by Draco:\n\n"
  echo -e "Also try 'slurmhelp'\n"
  echo "allow_file_to_age - pause a program until a file is 'old'"
  echo "canonicalize_filename - standardize paths to files"
  echo "cleanemacs    - remove ~ and .elc files."
  echo "die           - exit with a message"
  echo "dracoenv/rmdracoenv - load/unload the draco environment"
  echo "establish_permissions - Change group to othello, dacodes or draco and change permissions to g+rwX,o-rwX"
  echo "flavor        - build a string that looks like fire-openmpi-2.0.2-intel-17.0.1"
  echo "fn_exists     - return true if named bash function is defined"
  echo "install_verions - helper for doing releases (see release_toss2.sh)"
  echo "lookupppn     - return PE's per node."
  echo "machineName   - return a string to represent the current machine."
  echo "matches_extension - Does the provided filename have a matching extension?"
  echo "npes_build    - return PE's to be used for compiling."
  echo "npes_test     - return PE's to be used for testing."
  echo "osName        - return a string to represent the current machine's OS."
  echo "proxy         - toggle the status of the LANL proxy variables."
  echo "publish_release - helper for doing releases (see release_toss2.sh)"
  echo "qrm           - quick rm for directories located in lustre scratch spaces."
  echo "rdde          - more agressive reset of the draco environment."
  echo "run           - echo a command and then run it."
  echo "selectscratchdir - find a scratch drive"
  echo "version_gt    - compare versions"
  echo "whichall      - Find all matchs."
  echo "xfstatus      - print status 'transfered' files."
  echo -e "\nUse 'type <function>' to print the full content of any function.\n"
}

#--------------------------------------------------------------------------------------------------#
# Function Definitions
#--------------------------------------------------------------------------------------------------#

# Print an error message and exit.
# e.g.: cd $dir || die "can't change dir to $dir".
function die () { echo " "; echo "FATAL ERROR: $1"; exit 1;}

# Echo a command and then run it.
function run ()
{
  echo "==> $1"; if [[ "${dry_run:-no}" == "no" ]]; then eval "$1"; fi
}

# Return 0 if provided name is a bash function.
function fn_exists ()
{
  type "$1" 2>/dev/null | grep -c 'is a function'
}

#--------------------------------------------------------------------------------------------------#
# The script starts here
#--------------------------------------------------------------------------------------------------#

# Logic taken from /usr/projects/hpcsoft/templates/header
function machineName
{
  sysName=${sysName="unknown"}
  if [[ -f /usr/projects/hpcsoft/utilities/bin/sys_name ]]; then
    sysName=$(/usr/projects/hpcsoft/utilities/bin/sys_name)
  elif [[ -d /projects/darwin ]]; then
    sysName=darwin
  elif [[ -d /usr/gapps/jayenne ]]; then
    case $(uname -n) in
      sq* | sequoia*) sysName=sq ;;
      rzansel*) sysName=rzansel ;;
      rzmanta*) sysName=rzmanta ;;
      sierra*) sysName=sierra ;;
    esac
  else
    sysName=$(uname -n | sed -e 's/[.].*//')
  fi
  if [[ "$sysName" == "unknown" ]]; then
    echo "Unable to determine machine name, please edit tools/common.sh."
    return 1
  fi
  echo "$sysName"
}

# Logic taken from /usr/projects/hpcsoft/templates/header
function osName
{
  osName=${osName="unknown"}
  if [[ -f /usr/projects/hpcsoft/utilities/bin/sys_os ]]; then
    if ! osName=$(/usr/projects/hpcsoft/utilities/bin/sys_os); then
      osName="unknown"
    fi
  fi
  if [[ ${osName} == "unknown" ]]; then
    if [[ -d /projects/darwin ]] ; then
      osName=darwin
      if [[ -f /projects/draco/vendors/bin/target_arch ]]; then
        osName=darwin-$(/projects/draco/vendors/bin/target_arch)
      fi
    elif [[ -d /usr/gapps/jayenne ]]; then
      osName=$(uname -p)
    fi
  fi
  if [[ ${osName} == "unknown" ]]; then
    osName=$(uname -p)
  fi
  if [[ "$osName" == "unknown" ]]; then
    echo "Unable to determine system OS, please edit scripts/common.sh."
    return 1
  fi
  echo "$osName"
}

function machineFamily
{
  machfam=$(uname -p)
  case $(osName) in
    darwin*) machfam=$(osName) ;;
    cle*)
      case $(uname -n | sed -e 's/[.].*//') in
        tt* | tr* ) machfam="ats1" ;;
        cp* | th* ) machfam="cray" ;;
      esac
      ;;
    toss3) machfam="cts1" ;;
    ppc64le)
      case $(uname -n | sed -e 's/[.].*//') in
        rzans* | sierra* ) machfam="ats2" ;;
      esac
      ;;
  esac
  echo $machfam
}

#--------------------------------------------------------------------------------------------------#
# Generates a string of the form <platform>-<mpi+ver>-<compiler+ver>
function flavor
{
  platform=$(machineName)
  os=$(osName)
  case $os in
    toss*)
      if [[ $LMPI ]]; then
        mpiflavor="$LMPI-$LMPIVER"
      else
        mpiflavor="unknown"
      fi
      if [[ $LCOMPILER ]]; then
        compilerflavor="$LCOMPILER-$LCOMPILERVER"
      else
        compilerflavor="unknown"
      fi
      ;;
    cle*)
      if [[ $CRAY_MPICH2_VER ]]; then
        mpiflavor=mpt-$CRAY_MPICH2_VER
      else
        mpiflavor="unknown"
      fi
      # Try to determine the loaded compiler
      loadedmodules="${LOADEDMODULES:-notset}"
      # select modules that look like compilers
      unset compilermodules
      for module in ${loadedmodules//:/ }; do
        case $module in
          PrgEnv*)
          # Ingore PrgEnv matches.
          ;;
          cce/* | cray/* | gcc/* | gnu/* | intel/* | pgi/* )
            tmp="${module////-}"
            compilermodules="$tmp $compilermodules"
            ;;
        esac
      done
      # pick the first compiler in the list
      compilerflavor=$(echo "${compilermodules}" | awk '{print $1;}')
      # append target if KNL
      if [[ $(echo "${CRAY_CPU_TARGET}" | grep -c knl) == 1 ]]; then
        compilerflavor+='-knl'
      fi
      ;;
    darwin*)
      platform=$os # e.g.: darwin-power9
      if [[ $(which mpirun) =~ openmpi ]]; then
        if [[ $MPIRUN ]]; then
          LMPIVER=$(echo "${MPIRUN}" | sed -r 's%.*/([0-9]+)[.]([0-9]+)[.]([0-9]+).*%\1.\2.\3%')
        else
          LMPIVER='unknown'
        fi
        mpiflavor="openmpi-$LMPIVER"
      elif [[ $(which mpirun) =~ smpi ]]; then
        if [[ $MPIRUN ]]; then
          LMPIVER=$(echo "${MPIRUN}" | sed -r 's%.*-([0-9]+)[.]([0-9]+)[.]([0-9]+)[.]([0-9]+).*%\1.\2.\3.\4%')
        else
          LMPIVER='unknown'
        fi
        mpiflavor="spectrum-mpi-$LMPIVER"
      else
        mpiflavor="unknown"
      fi
      if [[ $CC ]]; then
        if [[ $CC =~ "gcc" ]]; then
          compilerflavor=gnu-$("${CC}" --version | head -n 1 | sed -r 's%.* %%')
        elif [[ $CC =~ "icc" ]]; then
          compilerflavor=intel-$("${CC}" --version | head -n 1 | awk '{ print $3 }' | sed -r 's%([0-9]+).([0-9]+).([0-9]+).*%\1.\2.\3%')
        elif [[ $CC =~ "IBM" ]] || [[ $CC =~ "xlc" ]]; then
          compilerflavor=xl-$("${CC}" --version | grep Version | sed -r 's%.* %%' | sed -r 's%0+([1-9])%\1%g')
        fi
      else
        compilerflavor="unknown"
      fi
      ;;
    ppc64)
      # more /bgsys/drivers/V1R2M3/ppc64/comm/include/mpi.h
      # | grep MPI_VERSION    ==> 2 ==> (mpich2)
      # | grep MPICH2_VERSION ==> 1.5
      mpiflavor="mpich2-1.5"

      case $CC in
        *gcc*)
          LCOMPILER=gnu
          LCOMPILERVER=$("${CC}" --version | head -n 1 | sed -e 's/.*\([0-9][.][0-9][.][0-9]\)/\1/')
          compilerflavor="$LCOMPILER-$LCOMPILERVER"
          ;;
        *xlc*)
          LCOMPILER=ibm
          LCOMPILERVER=$("$CC" -V | head -n 1 | sed -e 's/.*[/]\([0-9]\+\.[0-9]\).*/\1/')
          compilerflavor="$LCOMPILER-$LCOMPILERVER"
          ;;
        *)
          compilerflavor=unknown-unknown ;;
      esac
      ;;

    ppc64le)
      if [[ -n $OPAL_PREFIX ]]; then
        mpiflavor="${OPAL_PREFIX//*\//}"
      elif [[ -n $LMOD_MPI_NAME ]]; then
        mpiver="${LMOD_MPI_VERSION//-*/}"
        mpiflavor="$LMOD_MPI_NAME-$mpiver"
      else
        mpiflavor=unknown
      fi
      case $CC in
        *gcc*)
          LCOMPILER=gnu
          LCOMPILERVER=$("$CC" --version | head -n 1 | sed -e 's/.*\([0-9][.][0-9][.][0-9]\).*/\1/')
          compilerflavor="$LCOMPILER-$LCOMPILERVER"
          ;;
        *xlc*)
          LCOMPILER=ibm
          LCOMPILERVER=$("$CC" --version | head -n 1 | sed -e 's/.*V\([0-9][0-9][.][0-9][.][0-9]\) .*/\1/')
          compilerflavor="$LCOMPILER-$LCOMPILERVER"
          ;;
        *) compilerflavor=unknown-unknown ;;
      esac
      ;;

    *)
      # CCS-NET machines or generic Linux?
      if [[ $MPI_NAME ]]; then
        mpiflavor="$MPI_NAME-$MPI_VERSION"
      elif [[ $LMOD_MPI_NAME ]]; then
        mpiflavor="$LMOD_MPI_NAME"
        if [[ $LMOD_MPI_VERSION ]]; then
          mpiflavor+="-${LMOD_MPI_VERSION//-*/}"
        fi
      else
        mpiflavor="unknown"
      fi
      if [[ $LCOMPILER ]]; then
        compilerflavor="$LCOMPILER-$LCOMPILERVER"
      else
        compilerflavor="unknown"
      fi
      ;;
  esac
  echo "$platform-$mpiflavor-$compilerflavor"
}

#--------------------------------------------------------------------------------------------------#
# returns a path to a directory
function selectscratchdir ()
{
  local scratchdirs
  scratchdirs=$(df --output=pcent,target 2>&1 | grep -c unrecognized)
  if [[ $scratchdirs == 0 ]]; then
    scratchdirs=$(df --output=pcent,target 2> /dev/null | grep scratch | grep -v netscratch | sort -g)
    if [[ -z "$scratchdirs" ]]; then
      scratchdirs=$(df --output=pcent,target 2> /dev/null | grep workspace | sort -g)
    fi;
  else
    scratchdirs=$(df -a 2> /dev/null | grep net/scratch | awk '{ print $4 " "$5 }' | sort -g)
    if [[ -z "$scratchdirs" ]]; then
      scratchdirs=$(df -a 2> /dev/null | grep lustre/scratch | awk '{ print $4 " "$5 }' | sort -g)
    fi;
  fi;
  for item in $scratchdirs; do
    if [[ $(echo "$item" | grep -c %) == 1 ]]; then
      continue;
    fi;
    if [[ -d $item/users ]]; then
      item2="${item}/users";
      mkdir -p "$item2/$USER" &>/dev/null;
      if [[ -w $item2/$USER ]]; then
        echo "$item2";
        return;
      fi;
    fi;
    mkdir -p "$item/$USER" &>/dev/null;
    if [[ -w $item/$USER ]]; then
      echo "$item";
      return;
    fi;
    mkdir -p "$item/yellow/$USER" &>/dev/null;
    if [[ -w $item/yellow/$USER ]]; then
      echo "$item/yellow";
      return;
    fi;
  done;
  item="/netscratch/$USER;"
  mkdir -p "$item" &>/dev/null;
  if [[ -w $item ]]; then
    echo "$item";
    return;
  fi
}

#--------------------------------------------------------------------------------------------------#
function lookupppn()
{
  # https://hpc.lanl.gov/index.php?q=summary_table
  local target
  local ppn
  target=$(uname -n | sed -e 's/[.].*//')
  ppn=1
  case "${target}" in
    t[rt]-fe* | t[rt]-login*)
      if [[ "$CRAY_CPU_TARGET" == "haswell" ]]; then
        ppn=32
      elif [[ "$CRAY_CPU_TARGET" == "knl" ]]; then
        ppn=68
      else
        echo "ERROR: Expected CRAY_CPU_TARGET to be set in the environment."
        exit 1
      fi
      ;;
    fi* | ic* | sn* | cy* )         ppn=36 ;;
    rzansel* | rzmanta* | sierra* ) ppn=40 ;;
    *) ppn=$(grep -c processor /proc/cpuinfo) ;;
  esac
  echo "$ppn"
}

function npes_build
{
  local np
  np=1
  if [[ ${PBS_NP} ]]; then
    np="${PBS_NP}"
  elif [[ ${SLURM_NPROCS} ]]; then
    np="${SLURM_NPROCS}"
  elif [[  ${SLURM_CPUS_ON_NODE} ]]; then
    np="${SLURM_CPUS_ON_NODE}"
  elif [[ ${SLURM_TASKS_PER_NODE} ]]; then
    np="${SLURM_TSKS_PER_NODE}"
  elif [[ ${LSB_DJOB_NUMPROC} ]]; then
    np="${LSB_DJOB_NUMPROC}"
  elif [[ -f /proc/cpuinfo ]]; then
    # lscpu=`lscpu | grep "CPU(s):" | head -n 1 | awk '{ print $2 }'`
    np=$(grep -c processor /proc/cpuinfo)
  fi
  echo "$np"
}

function npes_test
{
  # assume no parallel testing capability.
  local np
  np=1

  # allow specialization per machine (if needed)
  local target
  target=$(uname -n | sed -e 's/[.].*//')
  case "${target}" in

    # current LSF only allows one executable to run under 'jsrun' at a time.
    # rzansel* | rzmanta* | sierra*) ppn=1 ;;

    *)

      # use lscpu if it is available.
      if ! [[ $(which lscpu 2>/dev/null) == 0 ]]; then
        local cps
        local ns
        # number of cores per socket
        cps=$(lscpu | grep "^Core(s)" | awk '{ print $4 }')
        # number of sockets
        ns=$(lscpu | grep "^Socket(s):" | awk '{ print $2 }')
        np=$(( cps * ns))

      else

        if [[ ${PBS_NP} ]]; then
          np="${PBS_NP}"
        elif [[ ${SLURM_NPROCS} ]]; then
          np="${SLURM_NPROCS}"
        elif [[  ${SLURM_CPUS_ON_NODE} ]]; then
          np="${SLURM_CPUS_ON_NODE}"
        elif [[ ${SLURM_TASKS_PER_NODE} ]]; then
          np="${SLURM_TSKS_PER_NODE}"
        elif [[ $(uname -p) == "ppc" ]]; then
          # sinfo --long --partition=pdebug (show limits)
          np=64
        elif [[ -f /proc/cpuinfo ]]; then
          # lscpu=`lscpu | grep "CPU(s):" | head -n 1 | awk '{ print $2 }'`
          np=$(grep -c processor /proc/cpuinfo)
        fi

      fi
      ;;
  esac
  echo "$np"
}

#--------------------------------------------------------------------------------------------------#
# Configure, Build and Run the Tests
#--------------------------------------------------------------------------------------------------#
function install_versions
{
  local config_step
  local build_step
  local test_step
  config_step=0
  build_step=0
  test_step=0

  # shellcheck disable=2154
  if [[ -z ${steps} ]]; then
    echo "You must provide variable steps."
    echo "E.g.: steps=\"configure build test\""
    return
  else
    for s in $steps; do
      case $s in
        config) config_step=1 ;;
        build)  build_step=1  ;;
        test)   test_step=1   ;;
      esac
    done
  fi
  if [[ -z ${buildflavor} ]]; then
    buildflavor=$(flavor)
  fi
  # shellcheck disable=2154
  if [[ -z ${rttversion} ]]; then
    echo "You must provide variable rttversion."
    # echo "E.g.: VERSIONS=( \"debug\" \"opt\" )"
    return
  fi
  # shellcheck disable=2154
  if [[ -z "${options}" ]]; then
    echo "You must provide variable option."
    return
  fi
  # shellcheck disable=2154
  if [[ -z "${package}" ]]; then
    echo "You must provide variable package."
    echo "E.g.: package=\"draco\""
    return
  fi
  # shellcheck disable=2154
  if [[ -z "${install_prefix}" ]]; then
    echo "You must provide variable install_prefix."
    echo "E.g.: install_prefix=/usr/projects/draco/${pdir:-notset}/$buildflavor"
    return
  fi
  if ! [[ ${build_pe} ]]; then
    build_pe=$(npes_build)
  fi
  if ! [[ ${test_pe} ]]; then
    test_pe=$(npes_test)
  fi

  # Echo environment before we start:
  echo
  echo
  echo "# Environment"
  echo "# ------------"
  echo
  run "module list"
  run "printenv"
  echo "---"
  echo "Environment size = $(printenv | wc -c)"

  echo
  echo
  echo "# Begin release: $buildflavor/$rttversion"
  echo "# ------------"
  echo

  # Create install directory
  install_dir="$install_prefix/$rttversion"
  if ! [[ -d $install_dir ]]; then
    run "mkdir -p $install_dir" || (echo "Could not create $install_dir" && return 1)
  fi

  # try to locate the souce directory
  if ! [[ -d "${source_prefix:=notset}" ]]; then
    echo "You must set the variable 'source_prefix'"
    return
  fi
  if [[ -f "${source_prefix}/source/ChangeLog" ]]; then
    source_dir="$source_prefix/source"
  else
    local possible_source_dirs
    possible_source_dirs=$(/bin/ls -1 "$source_prefix/source")
    for dir in $possible_source_dirs; do
      if [[ -f "$source_prefix/source/$dir/ChangeLog" ]]; then
        source_dir="$source_prefix/source/$dir"
        break
      fi
    done
  fi
  if ! [[ -f "$source_dir/CMakeLists.txt" ]]; then
    echo "Could not find sources. Tried looking at $source_prefix/source/"
    exit 1
  fi
  # source_dir="$source_prefix/source/$package"
  build_dir="${build_prefix:=notset}/${rttversion:=notset}/${package:0:1}"

  # Purge any existing files before running cmake to configure the build directory.
  if [[ "$config_step" == 1 ]]; then
    if [[ -d "${build_dir}" ]]; then
      run "rm -rf ${build_dir}"
    fi
    run "mkdir -p $build_dir" || (echo "Could not create directory $build_dir." && return 2)
  fi

  run "cd $build_dir"
  if [[ "$config_step" == 1 ]]; then
    run "cmake -DCMAKE_INSTALL_PREFIX=$install_dir $options $CONFIG_EXTRA $source_dir" \
      || (echo "Could not configure in $build_dir from source at $source_dir" && return 3)
  fi
  if [[ "$build_step" == 1 ]]; then
    run "${MAKE_COMMAND:-make} -j $build_pe -l $build_pe install"  \
      || (echo "Could not build code/tests in $build_dir" && return 4)
  fi
  if [[ "$test_step" == 1 ]]; then
    # run all tests
    run "${CTEST_COMMAND:-ctest} -j $test_pe --output-on-failure" \
      || run "${CTEST_COMMAND:-ctest} -j $test_pe --output-on-failure --rerun-failed"
  fi
  wait
  if ! [[ "${build_permissions:-notset}" = "notset" ]]; then
    run "chmod -R $build_permissions $build_dir"
  fi
}

#--------------------------------------------------------------------------------------------------#
# If $jobids is set, wait for those jobs to finish before setting groups and permissions.
function publish_release()
{
  echo " "
  echo "Waiting batch jobs to finish ..."
  echo "   Running jobs = ${jobids:-no jobs}"

  case $(osName) in
    toss* | cle* ) SHOWQ=squeue ;;
    darwin| ppc64) SHOWQ=squeue ;;
    ppc64le )      SHOWQ=bjobs  ;;
  esac

  # wait for jobs to finish
  for jobid in $jobids; do
    while [[ $("${SHOWQ}" | grep -c "$jobid") -gt 0 ]]; do
      "${SHOWQ}" | grep "$jobid"
      sleep 5m
    done
    echo "   Job $jobid is complete."
  done

  echo " "
  echo "Updating file permissions ..."
  echo " "

  # shellcheck disable=2154
  if [[ -n $install_permissions ]]; then
    # Set access to top level install dir.
    if [[ -d $install_prefix ]]; then
      run "chgrp -R ${install_group:-ccsrad} $source_prefix"
      run "chmod -R $install_permissions $source_prefix"
      run "find $source_prefix -type d -exec chmod g+s {} +"
    fi
  fi
}

#--------------------------------------------------------------------------------------------------#
# Pause until the 'last modified' timestamp of file $1 to be $2 seconds old.
function allow_file_to_age
{
  if [[ ! $2 ]]; then
    echo "ERROR: This function requires two arguments: a filename and an age value (sec)."
    exit 1
  fi

  # If file does not exist, no need to wait.
  if [[ ! -f $1 ]]; then
    return
  fi

  local timediff
  local print_message
  local pr_last_check
  local timediff
  local timeleft

  # assume file was last modified 0 seconds ago.
  timediff=0

  # If no changes for $2 seconds, continue
  # else, wait until until file, $1, hasn't been touched for $2 seconds.
  print_message=1
  while [[ $timediff -lt $2 ]]; do
    eval "$(date +'now=%s')"
    pr_last_check=$(date +%s -r "$1")
    # shellcheck disable=2154
    timediff=$(("$now"-"$pr_last_check"))
    timeleft=$(("$2"-"$timediff"))
    if [[ $timeleft -gt 0 ]]; then
      if [[ $print_message == 1 ]]; then
        echo "The log file $1 was recently modified."
        echo "To avoid colliding with another running test we are waiting"
        print_message=0
      fi
      echo "... $timeleft seconds"
    fi
    sleep 30s
  done
}


#--------------------------------------------------------------------------------------------------#
# Canonicalize_filename:
# See http://stackoverflow.com/questions/1055671/how-can-i-get-the-behavior-of-gnus-readlink-f-on-a-mac
#--------------------------------------------------------------------------------------------------#
canonicalize_filename () {
  local target_file
  local physical_directory
  local result
  target_file=$1
  physical_directory=""
  result=""

  # Need to restore the working directory after work.
  pushd "$(pwd)" > /dev/null || \
    (echo "canonicalize_filename error 1: cannot access PWD" && return 1)

  cd "$(dirname "$target_file")" || \
    (echo "canonicalize_filename error 2: cannot access target_file PWD" && return 2)
  target_file=$(basename "$target_file")

  # Iterate down a (possible) chain of symlinks
  while [ -L "$target_file" ]; do
    target_file=$(readlink "$target_file")
    cd "$(dirname "$target_file")" || \
      (echo "canonicalize_filename error 3: cannot access target_file PWD" && return 3)
    target_file=$(basename "$target_file")
  done

  # Compute the canonicalized name by finding the physical path for the directory we're in and
  # appending the target file.
  physical_directory=$(pwd -P)
  result="$physical_directory"/"$target_file"

  # restore the working directory after work.
  popd > /dev/null || (echo "canonicalize_filename error 4: cannot popd" && return 4)

  echo "$result"
}

#--------------------------------------------------------------------------------------------------#
# Compare versions
#
# Example:
#
# EMACSVER=`emacs --version | head -n 1 | sed -e 's/.*Emacs //'`
# if `version_gt "24.0.0" $EMACSVER` ; then echo yes; fi
#--------------------------------------------------------------------------------------------------#
function version_gt()
{
  test "$(printf '%s\n' "$@" | sort -V | head -n 1)" != "$1";
}

#--------------------------------------------------------------------------------------------------#
# check whether the given file matches any of the set extensions
#
# Example:
#
# FILE_EXTS=".f90 .F90"
# FILE_ENDINGS="_f.h _f77.h _f90.h"
# for file in $modified_files; do
#   if ! matches_extension $file; then
#      continue;
#   fi
#   <do stuff>
# done
#--------------------------------------------------------------------------------------------------#
matches_extension() {
  local filename
  local extension
  local end
  local ext
  filename=$(basename "$1")
  extension=".${filename##*.}"
  for end in $FILE_ENDINGS; do [[ "$filename" == *"$end" ]] && return 1; done
  for ext in $FILE_EXTS; do [[ "$ext" == "$extension" ]] && return 0; done
  return 1
}

#--------------------------------------------------------------------------------------------------#
# relcreatesymlinks
#
# Creates symlinks to an installation to satisfy requests from clients.  For example:
#
#   cts1-openmpi-2.1.2-gcc-7.4.0    -> snow-openmpi-2.1.2-gcc-7.4.0
#   grizzly-openmpi-2.1.2-gcc-7.4.0 -> snow-openmpi-2.1.2-gcc-7.4.0
#
# This function is designed to be called from release-<machine>.sh.  The following variables should
# be set in the environment:
#
# - source_prefix :: e.g. /usr/projects/draco/draco-NN_NN_NN
# - evironments   :: a list of bash function names. These functiosn should establish build
#                    environments (see cts1-env.sh, etc.)
# - siblings      :: a space delimited list of machine names that are considered to be equivalent
#                    (e.g.: "fire ice cyclone")
#
#--------------------------------------------------------------------------------------------------#
function relcreatesymlinks() {

  OLDPWD=$(pwd)
  machfam=$(machineFamily)

  if ! [[ -d $source_prefix ]]; then
    echo "relcreatesymlinks:: Expected variable source_prefix to be set" && return 1
  fi
  if [[ ${siblings:=notset} == "notset" ]]; then
    echo "relcreatesymlinks:: Expected variable siblings to be set" && return 2
  fi
  if [[ ${environments:=notset} == "notset" ]]; then
    echo "relcreatesymlinks:: Expected variable environments to be set" && return 3
  fi

  run "cd $source_prefix"

  for env in $environments; do

    if ! [[ $(fn_exists "$env") ]]; then
      echo -n "relcreatesymlinks:: Attempted to load environment $env, but this function is "
      echo "not defined". && return 4
    fi
    # establish environment
    $env
    buildflavor=$(flavor)

    # create symlinks
    me="${buildflavor//-*/}"
    familymatch="${siblings// /|}"
    case $me in
      @($familymatch) )
        linkname="${buildflavor//${me}/${machfam}}"
        run "ln -s $buildflavor $linkname"
        for m in $siblings; do
          linkname="${buildflavor//${me}/${m}}"
          linkname="${linkname//yellow/}"
          if ! [[ -d $source_prefix/$linkname ]]; then
            run "ln -s $buildflavor $linkname"
          fi
        done
        ;;
      esac
    done

  run "cd $OLDPWD"
}

 #--------------------------------------------------------------------------------------------------#
 export canonicalize_filename
 export die
 export flavor
 export install_versions
 export job_launch_sanity_checks
 export machineName
 export matches_extension
 export npes_build
 export npes_test
 export osName
 export relcreatesymlinks
 export run
 export selectscratchdir
 export version_gt

 #--------------------------------------------------------------------------------------------------#
 # End common.sh
 #--------------------------------------------------------------------------------------------------#
