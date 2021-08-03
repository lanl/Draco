#!/bin/bash
#-*- Mode: bash -*-
#--------------------------------------------------------------------------------------------------#
# File  : environment/bashrc/.bashrc
# Date  : Tuesday, May 31, 2016, 14:48 pm
# Author: Kelly Thompson
# Note  : Copyright (C) 2016-2020, Triad National Security, LLC., All rights are reserved.
#
# Bash configuration file upon bash shell startup
#
# Instructions (customization):
#
# 1. Setup
#    - Copy logic from draco/environment/bashrc/sample.bashrc and
#      draco/environment/bashrc/sample.bash_profile.
# 2. Override settings using the code found in the sample.bashrc.
#--------------------------------------------------------------------------------------------------#

#uncomment to debug this script.
#export verbose=true

#--------------------------------------------------------------------------------------------------#
# ENVIRONMENTS for interactive sessions
#--------------------------------------------------------------------------------------------------#

# If this is an interactive shell then the environment variable $- should contain an "i":
case ${-} in
  *i*)
    export INTERACTIVE=true
    [[ "${verbose:=false}" == "true" ]] && echo "in draco/environment/bashrc/.bashrc"

    # Shell options
    shopt -s checkwinsize # autocorrect window size
    shopt -s cdspell      # autocorrect spelling errors on cd command line.
    shopt -s histappend   # append to the history file, don't overwrite it

    # More environment setup ----------------------------------------------------------------------#

    # don't put duplicate lines or lines starting with space in the history. See bash(1) for more
    # options
    HISTCONTROL=ignoreboth

    # for setting history length see HISTSIZE and HISTFILESIZE in bash(1)
    HISTSIZE=1000
    HISTFILESIZE=2000

    # Prevent creation of core files (ulimit -a to see all limits).
    # ulimit -c 0

    # colored GCC warnings and errors
    export GCC_COLORS='error=01;31:warning=01;35:note=01;36:caret=01;32:locus=01:quote=01'
    export TERM=xterm-256color

    # Attempt to find DRACO
    if ! [[ $DRACO_SRC_DIR ]]; then
      _BINDIR=$(dirname "${BASH_ARGV[0]}")
      DRACO_SRC_DIR=$(cd "$_BINDIR/../.." || exit; pwd)
      export DRACO_SRC_DIR
      export DRACO_ENV_DIR="${DRACO_SRC_DIR}/environment"
    fi
    [[ "${DRACO_ENV_DIR:-notset}" == "notset" ]] && \
      echo "DRACO_ENV_DIR not set, some draco env will be disabled."

    #----------------------------------------------------------------------------------------------#
    # Common aliases
    #----------------------------------------------------------------------------------------------#

    # help shellcheck linter.
    # shellcheck source=/dev/null
    source "${DRACO_ENV_DIR}/bashrc/bash_aliases.sh"

    # If this is an xterm set the title to user@host:dir
    case "$TERM" in
      xterm*|rxvt*) echo -ne "\033]0;${nodename:-xterm}\007" ;;
      *) ;;
    esac
    ;; # end case 'interactive'

  #------------------------------------------------------------------------------------------------#
  # ENVIRONMENTS for non interactive sessions
  #------------------------------------------------------------------------------------------------#

  *) # Not an interactive shell (e.g. A PBS shell?)
    export INTERACTIVE=false
    ;;
esac

#--------------------------------------------------------------------------------------------------#
# ENVIRONMENTS - bash functions (all interactive sessions)
#--------------------------------------------------------------------------------------------------#

# Bash functions are not inherited by subshells.
if [[ ${INTERACTIVE} ]]; then

  # Common bash functions and alias definitions
  # shellcheck source=/dev/null
  [[ -f "${DRACO_ENV_DIR}/bashrc/bash_functions.sh" ]] && \
    source "${DRACO_ENV_DIR}/bashrc/bash_functions.sh"
  # shellcheck source=/dev/null
  [[ -f "${DRACO_ENV_DIR}/../tools/common.sh" ]] && source "${DRACO_ENV_DIR}/../tools/common.sh"

  # aliases and bash functions for working with slurm
  if !  [[ $(which squeue 2>&1 | grep -c "no squeue") == 1 ]] &&
    [[ $(which squeue | grep -c squeue) -gt 0 ]]; then
    # shellcheck source=/dev/null
    [[ -f "${DRACO_ENV_DIR}/bashrc/bashrc_slurm" ]] && source "${DRACO_ENV_DIR}/bashrc/bashrc_slurm"
  fi
fi

#--------------------------------------------------------------------------------------------------#
# ENVIRONMENTS - once per login
#--------------------------------------------------------------------------------------------------#

if [[ ${INTERACTIVE} == true ]]; then

  # Append PATHS (not linux specific, not ccs2 specific).
  add_to_path "${DRACO_ENV_DIR}/bin"
  add_to_path "${DRACO_SRC_DIR}/tools"

  # Tell wget to use LANL's www proxy (see
  # trac.lanl.gov/cgi-bin/ctn/trac.cgi/wiki/SelfHelpCenter/ProxyUsage)
  # export http_proxy=http://wpad.lanl.gov/wpad.dat
  current_domain=$(awk '/^domain/ {print $2}' /etc/resolv.conf)
  #  found=`nslookup proxyout.lanl.gov | grep -c Name`
  #  if test ${found} == 1; then
  if [[ ${current_domain} == "lanl.gov" ]]; then
    export http_proxy=http://proxyout.lanl.gov:8080
    export https_proxy=$http_proxy
    export HTTP_PROXY=$http_proxy
    export HTTPS_PROXY=$http_proxy
    export no_proxy=".lanl.gov,127.0.0.1,localhost"
    export NO_PROXY=$no_proxy
  fi

  # cd paths - disable here, let user choose in ~/.bashrc
  CDPATH=

  # Hooks for clang-format as git commit hook:
  # Possible values: ON, TRUE, OFF, FALSE, DIFF (the default value is ON).
  export DRACO_AUTO_CLANG_FORMAT=ON

  # Parse the setup scripts, but don't actually load any modules.  This allows the developer to run
  # 'dracoenv' or 'rdde' later to load modules.
  if [[ -z "$DRACO_ENV_LOAD" ]]; then
    export DRACO_ENV_LOAD=ON
  fi

  # Silence warnings from GTK/Gnome
  export NO_AT_BRIDGE=1

  #------------------------------------------------------------------------------------------------#
  # ENVIRONMENTS - machine specific settings
  #------------------------------------------------------------------------------------------------#
  target=$(uname -n | sed -e s/[.].*//)

  case ${target} in

    # Darwin Heterogeneous Cluster (GPU, ARM, P9, etc.)
    # wiki: https://darwin.lanl.gov
    darwin-fe* | cn[0-9]*)
      # shellcheck source=/dev/null
      source "${DRACO_ENV_DIR}/bashrc/.bashrc_darwin_fe" ;;

    # Badger | Cyclone | Fire | Grizzly | Ice | Snow
    ba* | cy* | fi* | gr* | ic* | sn* )
      # shellcheck source=/dev/null
      source "${DRACO_ENV_DIR}/bashrc/.bashrc_cts1" ;;

    # wtrw and rfta
    red-wtrw* | rfta* | redcap* )
      # shellcheck source=/dev/null
      source "${DRACO_ENV_DIR}/bashrc/.bashrc_rfta" ;;

    # capulin, thunder, trinitite (tt-rfe) | trinity (tr-fe)
    cp-rfe* | cp-rlogin* | th-login* |tt-fey* | tt-rfe* | tt-login* | tr-fe* | tr-login* | nid* )
      # shellcheck source=/dev/null
      source "${DRACO_ENV_DIR}/bashrc/.bashrc_cray" ;;

    # LLNL ATS-2
    rzmanta* | rzansel* | sierra*)
      # shellcheck source=/dev/null
      source "${DRACO_ENV_DIR}/bashrc/.bashrc_ats2" ;;

    # LLNL ATS-4
    rznevada*)
      # shellcheck source=/dev/null
      source "${DRACO_ENV_DIR}/bashrc/.bashrc_ats4" ;;

    # LAP Virtual Machine
    vc*)
      # shellcheck source=/dev/null
      source "${DRACO_ENV_DIR}/bashrc/.bashrc_vm" ;;

    ccscs[1-9]* | ccsnet[1-9]*)
      # shellcheck source=/dev/null
      source "${DRACO_ENV_DIR}/bashrc/.bashrc_ccsnet" ;;

    # Assume personal workstation
    *) export NoModules=1 ;;

  esac

  # shellcheck source=/dev/null
  [[ -f "${DRACO_ENV_DIR}/bashrc/bash_functions2.sh" ]] && \
    source "${DRACO_ENV_DIR}/bashrc/bash_functions2.sh"

fi

# provide some bash functions (dracoenv, rmdracoenv) for non-interactive
# sessions.

# shellcheck source=/dev/null
[[ -f "${DRACO_ENV_DIR}/bashrc/bash_functions2.sh" ]] && \
  source "${DRACO_ENV_DIR}/bashrc/bash_functions2.sh"

if [[ "${verbose:=false}" == "true" ]]; then
  echo "in draco/environment/bashrc/.bashrc ... done";
fi

#--------------------------------------------------------------------------------------------------#
# end of .bashrc
#--------------------------------------------------------------------------------------------------#
