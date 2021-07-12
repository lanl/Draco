#!/bin/bash
## -*- Mode: bash -*-
##-----------------------------------------------------------------------------------------------##
## File  : environment/bashrc/bash_functions.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2020, Triad National Security, LLC., All rights are reserved.
##-----------------------------------------------------------------------------------------------##
##
## Summary: Misc bash functions useful during development of code.
##
## Functions
## ---------
##
## whichall <string> - show path of all commands in path that match <string>
##
## cleanemacs        - recursively remove ~ files, .flc files and .rel files.
##
## findsymbol <sym>  - search all libraries (.so and .a files) in the current directory for symbol
##                     <sym>.
##
## npwd              - function used to set the prompt under bash.
##
## xfstatus          - report status of transfer.lanl.gov
##
## rm_from_path      - remove a directory from $PATH
##
## add_to_path       - add a directory to $PATH
##
## proxy             - (un)set http_proxy variables
##
## fn_exists         - test if a bash function is defined
##
## run               - echo then evaluate a bash command
##
## rdde              - reload the default draco environment
##
## qrm               - quick remove (for lustre filesystems).
##-----------------------------------------------------------------------------------------------##

# pragma once
if [[ $(type dbf_pragma_once 2>&1 | grep -c "bash: type: ") != 0 ]]; then
  function dbf_pragma_once () {
    echo "dbfpo"
  }

if [[ "${verbose:=false}" ==  "true" ]]; then
  echo "in draco/environment/bashrc/bash_functions.sh";
fi

##-----------------------------------------------------------------------------------------------##
## Find all matches in PATH (not just the first one)
##-----------------------------------------------------------------------------------------------##

function whichall ()
{
  for dir in ${PATH//:/ }; do
    if [ -x "$dir/$1" ]; then echo "$dir/$1"; fi;
  done
}
export -f whichall

##-----------------------------------------------------------------------------------------------##
## Recursively delete all ~ files.
##-----------------------------------------------------------------------------------------------##

function cleanemacs
{
  echo "Cleaning up XEmacs temporary and backup files.";
  find . -name '*~' -exec echo rm -rf {} \;
  find . -name '*~' -exec rm -rf {} \;
  find . -name '.*~' -exec echo rm -rf {} \;
  find . -name '.*~' -exec rm -rf {} \;
  # find . -name '*.flc' -exec echo rm -rf {} \;
  # find . -name '*.flc' -exec rm -rf {} \;
  # find . -name '*.rel' -exec echo rm -rf {} \;
  # find . -name '*.rel' -exec rm -rf {} \;
  echo "done."
}
export -f cleanemacs

##-----------------------------------------------------------------------------------------------##
## Used for formatting PROMPT.
## $HOME            -> ~
## ...scratch...    -> #
## .../projects/... -> @
##-----------------------------------------------------------------------------------------------##
parse_git_branch() {
  git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/ (\1)/'
}
export parse_git_branch

function npwd()
{
  # Optional arguments:
  #   $1 - number of chars to print.
  #   $2 - scratch location

  # How many characters of the $PWD should be kept
  local pwdmaxlen=40
  if [[ $1 ]]; then pwdmaxlen=$1; fi

  # local regHome=$(echo ${HOME} | sed -e 's/.*\///')

  local scratchdirs="/scratch:/lustre/ttscratch1:/lustre/scratch[34]/yellow"
  if [[ $2 ]]; then scratchdirs=$2; fi

  # Indicator that there has been directory truncation:
  local trunc_symbol="..."
  # substitute ~ for $HOME to shorten the full path
  # shellcheck disable=SC2001
  newPWD=$(echo "${PWD}" | sed -e "s%$HOME%~%")
  for dir in ${scratchdirs//:/ }; do
    newPWD="${newPWD//${dir}\/${USER}/#}"
  done

  local devdirs=/usr/projects/jayenne/devs
  newPWD="${newPWD//${devdirs}\/${USER}/@}"

  if [ ${#newPWD} -gt "$pwdmaxlen" ]; then
    local pwdoffset=$(( ${#newPWD} - pwdmaxlen ))
    newPWD="${trunc_symbol}${newPWD:$pwdoffset:$pwdmaxlen}"
  fi
  echo "${newPWD}"
}
export -f npwd

# Alternate
function npwd_alt()
{
  # Home directory
  local regHome=${HOME}
  # Home directory might be prepended with NFS mount location.
  local netHome=/net/machname.lanl.gov
  #   How many characters of the $PWD should be kept
  local pwdmaxlen=40
  #   Indicator that there has been directory truncation:
  local trunc_symbol="..."
  # substitute ~ for $HOME to shorten the full path
  newPWD="${PWD//${netHome}/}"
  newPWD="${newPWD//${regHome}/~}"
  if [[ ${#newPWD} -gt pwdmaxlen ]] ; then
    local pwdoffset=$(( ${#newPWD} - pwdmaxlen ))
    newPWD="${trunc_symbol}${newPWD:$pwdoffset:$pwdmaxlen}"
  fi
  echo "${newPWD}"
}
export -f npwd_alt

# Another alternate
# PROMPT_COMMAND='DIR=`pwd|sed -e "s!$HOME!~!"`; if [ ${#DIR} -gt 30 ]; then CurDir=${DIR:0:12}...${DIR:${#DIR}-15}; else CurDir=$DIR; fi'
# PS1="[\$CurDir] \$ "

##-----------------------------------------------------------------------------------------------##
## Usage:
##    findsymbol <symbol>
##
## Searches all .a and .so files in local directory for symbol
## <symbol>.  If found, the script provides the name of the library
## that contains the symbol.
##-----------------------------------------------------------------------------------------------##

function findsymbol()
{
  local nm_opt a_libs so_libs libs symbol gres
  nm_opt='-a'
  a_libs=$(/usr/bin/ls -1 ./*.a)
  if [[ -z "$a_libs" ]]; then a_libs=""; fi
  so_libs=$(/usr/bin/ls -1 ./*.so)
  if [[ -z "$so_libs" ]]; then so_libs=""; fi
  libs="$a_libs $so_libs"
  echo -e "\nSearching..."
  symbol=" T "
  for lib in $libs; do
    gres=$(nm "$nm_opt" "$lib" | grep "$1" | grep "$symbol")
    if ! [[ "$gres" == "" ]]; then
      echo-e "\nFound \"$symbol\" in $lib:\n     $gres"
    fi
  done
  echo " "
}
export -f findsymbol

##-----------------------------------------------------------------------------------------------##
## Transfer 2.0 (Mercury replacement)
## Ref: http://transfer.lanl.gov
##
## Examples:
##   xfpush foo.txt
##   xfstatus
##   xfpull foo.txt
##-----------------------------------------------------------------------------------------------##

function xfstatus()
{
  ssh red@transfer.lanl.gov myfiles
}
export -f xfstatus

##-----------------------------------------------------------------------------------------------##
## If string is found in PATH, remove it.
##-----------------------------------------------------------------------------------------------##
function rm_from_path ()
{
  export PATH="${PATH//${1}/}"
}
export -f rm_from_path

##-----------------------------------------------------------------------------------------------##
## If path is a directory add it to PATH (if not already in PATH)
##
## Use:
##   add_to_path <path> TEXINPUTS|BSTINPUTS|BIBINPUTS|PATH
##-----------------------------------------------------------------------------------------------##
function add_to_path ()
{
  case $2 in
    TEXINPUTS)
      if [ -d "$1" ] && [[ ":${TEXINPUTS}:" != *":$1:"* ]]; then
        TEXINPUTS="${TEXINPUTS:+${TEXINPUTS}:}$1"; fi ;;
    BSTINPUTS)
      if [ -d "$1" ] && [[ ":${BSTINPUTS}:" != *":$1:"* ]]; then
        BSTINPUTS="${BSTINPUTS:+${BSTINPUTS}:}$1"; fi ;;
    BIBINPUTS)
      if [ -d "$1" ] && [[ ":${BIBINPUTS}:" != *":$1:"* ]]; then
        BIBINPUTS="${BIBINPUTS:+${BIBINPUTS}:}$1"; fi ;;
    *)
      if [ -d "$1" ] && [[ ":${PATH}:" != *":$1:"* ]]; then
        PATH="${PATH:+${PATH}:}$1"; fi ;;
  esac
}
export -f add_to_path

##-----------------------------------------------------------------------------------------------##
## Toggle LANL proxies on/off
## https://wiki.archlinux.org/index.php/proxy_settings
##-----------------------------------------------------------------------------------------------##
function proxy()
{
  if [[ ! ${http_proxy} ]]; then
    # proxies not set, set them
    export http_proxy=http://proxyout.lanl.gov:8080
    export https_proxy=$http_proxy
    export HTTP_PROXY=$http_proxy
    export HTTPS_PROXY=$http_proxy
    # export http_no_proxy="*.lanl.gov"
    export no_proxy="localhost,127.0.0.1,.lanl.gov"
    export NO_PROXY=$no_proxy
  else
    # proxies are set, kill them
    unset http_proxy
    unset https_proxy
    unset HTTP_PROXY
    unset HTTPS_PROXY
    #unset http_no_proxy
    unset no_proxy
    unset NO_PROXY
  fi
}
export -f proxy

##-----------------------------------------------------------------------------------------------##
## Test to determine if named bash function exists in the current environment.
##-----------------------------------------------------------------------------------------------##
function fn_exists()
{
  type "$1" 2>/dev/null | grep -q 'is a function'
  res=$?
  echo $res
  return $res
}
export -f fn_exists

##-----------------------------------------------------------------------------------------------##
## Echo commands before execution (used in scripts)
##-----------------------------------------------------------------------------------------------##
function run () {
  echo "$1"
  if [[ "${dry_run:-notset}" != "notset" ]]; then eval "$1"; fi
}
export -f run

##-----------------------------------------------------------------------------------------------##
## Reset the draco developer environment
##-----------------------------------------------------------------------------------------------##
function rdde ()
{
  unset DRACO_BASHRC_DONE
  # shellcheck source=environment/bashrc/.bashrc
  source "${DRACO_ENV_DIR}/bashrc/.bashrc"
}
export -f rdde

#--------------------------------------------------------------------------------------------------#
# Quick remove: instead of 'rm -rf', mv the directory to .../trash/tmpname
#
# Use: qrm dir1 dir2 ...
#--------------------------------------------------------------------------------------------------#
function qrm ()
{
  # must provide at least one directory
  if [[ ! $1 ]]; then
    printf "You must provide a single argument to this function that is the name of the "
    printf "directory to delete.\n"
    return
  fi

  for dir in "$@"; do

    # fully qualified directory name
    fqd=$(cd "$dir" && pwd)

    if [[ ! -d "$dir" ]]; then
      echo "$dir doesn't look like a directory. aborting."
      return
    fi

    if [[ $(echo "$fqd" | grep -c scratch) == 0 ]]; then
      echo "This command should only be used for scratch directories."
      return
    fi

    # Identify the scratch system
    trashdir="${fqd//${USER}*/${USER}/trash}"

    # ensure trash folder exists.
    if ! [[ -d "$trashdir" ]]; then mkdir -p "$trashdir"; fi
    if ! [[ -d "$trashdir" ]]; then
      echo "FATAL ERROR: Unable access trashdir = $trashdir"
      return
    fi

    # We rename/move the old directory to a random name
    TMPDIR=$(mktemp -d "$trashdir/XXXXXXXXXX") || { echo "Failed to create temporary directory"; \
      return; }

    # do it
    mv "$dir" "$TMPDIR/."

  done
}
export -f qrm

if [[ "${verbose:=false}" ==  "true" ]]; then
  echo "in draco/environment/bashrc/bash_functions.sh ... done";
fi

fi # end pragma once

#--------------------------------------------------------------------------------------------------#
# End environment/bashrc/bash_functions
#--------------------------------------------------------------------------------------------------#
