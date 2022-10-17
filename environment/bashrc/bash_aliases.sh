#!/bin/bash
##-*- Mode: bash -*-
#--------------------------------------------------------------------------------------------------#
## File  : environment/bashrc/bashrc_aliases.sh
## Date  : Tuesday, Sep 01, 2020, 19:01 pm
## Author: Kelly Thompson <kgt@lanl.gov>
## Note  : Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved.
##
## bashrc_aliases is sourced by interactive shells from .bash_profile and .bashrc
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
# Draco Dev Env Customizations
#--------------------------------------------------------------------------------------------------#

# Generic Settings
alias ll='\ls -Flh'
alias lt='\ls -Flth'
alias ls='\ls -F'
alias la='\ls -A'
alias l.='\ls -hd .*'
alias lt.='ls -Flth .*'

# alias a2ps='a2ps --sides=duplex --medium=letter'
alias btar='tar --use-compress-program /usr/bin/bzip2'
alias cpuinfo='cat /proc/cpuinfo'
alias df='df -h'
alias dirs='dirs -v'
alias du='du -h --max-depth=1 --exclude=.snapshot'
alias less='/usr/bin/less -r'
alias mdstat='cat /proc/mdstat'
alias meminfo='cat /proc/meminfo'
alias mroe='more'
nodename=$(uname -n | sed -e 's/[.].*//g')
alias resettermtitle='echo -ne "\033]0;${nodename}\007"'

# Module related:
alias moduel='module'
alias ma='module avail'
alias mica='module --ignore_cache avail'
alias mls='module list'
alias mld='module load'
alias mul='module unload'
alias msh='module show'

# machines
# alias sierra='ssh -t redcap ssh sierra.llnl.gov'
alias rzansel='ssh -t ihpc-gate.lanl.gov ssh rzansel.llnl.gov'
alias rznevada='ssh -t ihpc-gate.lanl.gov ssh rznevada.llnl.gov'
alias rzvernal='ssh -t ihpc-gate.lanl.gov ssh rzvernal.llnl.gov'

#--------------------------------------------------------------------------------------------------#
# Color Prompt
#--------------------------------------------------------------------------------------------------#

# Provide special ls commands if this is a color-xterm or compatible terminal.

# 1. Does the current terminal support color?
if [ -x /usr/bin/tput ] && tput setaf 1 >&/dev/null; then
  # We have color support; assume it's compliant with Ecma-48 (ISO/IEC-6429). (Lack of such support
  # is extremely rare, and such a case would tend to support setf rather than setaf.)
  color_prompt=yes
fi

# 2. Override color_prompt for special values of $TERM
case "$TERM" in
  xterm-color|*-256color) color_prompt=yes;;
  emacs|dumb)
    color_prompt=no
    unset LS_COLORS
    ;;
esac

if [[ "${color_prompt:-no}" == "yes" ]]; then

  # Use custom colors if provided.
  if [[ -r ~/.dircolors ]]; then
    eval "$(dircolors -b ~/.dircolors)"
  else
    eval "$(dircolors -b)"
  fi

  # append --color option to some aliased commands

  alias ll='\ls -Flh --color'
  alias lt='\ls -Flth --color'
  alias ls='\ls -F --color'
  alias la='\ls -A --color'
  alias l.='\ls -hd --color .*'
  alias lt.='ls -Flth --color .*'

  alias grep='grep --color=auto'
  alias fgrep='fgrep --color=auto'
  alias egrep='egrep --color=auto'

  # colored GCC warnings and errors
  # export GCC_COLORS='error=01;31:warning=01;35:note=01;36:caret=01;32:locus=01:quote=01'

  # Colorized prompt (might need some extra debian_chroot stuff -- see wls example).
  if [ -z "${debian_chroot:-}" ] && [ -r /etc/debian_chroot ]; then
    debian_chroot=$(cat /etc/debian_chroot)
  fi

  if [ "$color_prompt" = yes ]; then
    PS1='${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ '
  else
    PS1='${debian_chroot:+($debian_chroot)}\u@\h:\w\$ '
  fi

fi
unset color_prompt

#--------------------------------------------------------------------------------------------------#
# User Customizations
#--------------------------------------------------------------------------------------------------#

# Aliases -----------------------------------------------------------------------------------------#

# alias cmakedebug='unset MPI_ROOT;C_FLAGS=-Werror CXX_FLAGS=-Werror cmake -Wdeprecated -Wdev'
# alias cmakedebugfast='cmake -DBUILD_TESTING=OFF'
# alias cmakerel='cmake -DCMAKE_BUILD_TYPE=RELEASE'
# alias ehco='echo'
# alias em='emacsclient -c --alternate-editor=emacs'
# alias emacs='/bin/emacs -fn 6x13 &> /dev/null'
# alias functions='declare -f | grep \^[a-zA-Z]'
# alias gitk='gitk --all'
# alias moduel='module'
# alias reportinvalidsymlink='find . -xtype l'
# alias rzansel='ssh -t ihpc-gate1.lanl.gov ssh rzansel.llnl.gov'
# alias vi='emacs -nw'

# Kerberos
# alias kerb='kinit -f -l 8h -kt ~/.ssh/cron.keytab $USER@lanl.gov'

# add_to_path /scratch/vendors/bin PATH

#--------------------------------------------------------------------------------------------------#
# End environment/bashrc/bashrc_aliases.sh
#--------------------------------------------------------------------------------------------------#
