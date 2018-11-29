#!/bin/bash
##-*- Mode: bash -*-
##---------------------------------------------------------------------------##
## File  : environment/bashrc/bash_functions2.sh
## Date  : Thursday, Dec 21, 2017, 12:14 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2017, Los Alamos National Security, LLC.
##         All rights are reserved.
##
##  Bash configuration file upon bash shell startup
##---------------------------------------------------------------------------##

# Generic functions for loading/unloaded the default set of modules.
fn_exists=`type dracoenv 2>/dev/null | head -n 1 | grep -c 'is a function'`
if test $fn_exists = 0; then
  # only define if they do not already exist...
  function dracoenv ()
  {
    for m in $dracomodules; do
	module load $m
    done
    unset MPI_ROOT
  }
  function rmdracoenv ()
  {
    # unload in reverse order.
    mods=( ${dracomodules} )
    for ((i=${#mods[@]}-1; i>=0; i--)); do
      loaded=`echo $LOADEDMODULES | grep -c ${mods[$i]}`
      if test $loaded = 1; then
        module unload ${mods[$i]}
      fi
    done
  }
fi

##---------------------------------------------------------------------------##
## end of environment/bashrc/bash_functions2.sh
##---------------------------------------------------------------------------##
