#!/bin/csh
##---------------------------------------------------------------------------##
## File  : environment/cshrc/.cshrc
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2017, Triad National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# Use: In ~/.cshrc add the following code:
#
# setenv DRACO_ENV_DIR ~/draco/environment
# source $DRACO_ENV_DIR/cshrc/.cshrc
#

if( ! $?DRACO_ENV_DIR )then
  if( -d /usr/projects/draco/environment )then
    setenv DRACO_ENV_DIR /usr/projects/draco/environment
  endif
  if( -d ~/draco/environment )then
    setenv DRACO_ENV_DIR ~/draco/environment
  endif
endif
if( ! $?DRACO_ENV_DIR )then
  echo "ERROR: DRACO_ENV_DIR has not been set."
  exit 1
endif
if( ! -d "${DRACO_ENV_DIR}" )then
  echo "ERROR: DRACO_ENV_DIR is not a valid directory (${DRACO_ENV_DIR})"
  exit 1
endif

# Setup PATH
setenv PATH $DRACO_ENV_DIR/bin:$PATH
if( -d ~/bin )then
    setenv PATH ~/bin:$PATH
endif

# Hooks for clang-format as git commit hook:
# Possible values: ON, TRUE, OFF, FALSE, DIFF (the default value is ON).
setenv DRACO_AUTO_CLANG_FORMAT ON

# Set term title
set myhost=`echo $HOST | sed -e 's/[.].*//g'`

# Aliases
alias btar 'tar --use-compress-program /usr/bin/bzip2'
alias cd.. 'cd ..'
alias cpuinfo 'cat /proc/cpuinfo'
alias df 'df -h'
alias dirs 'dirs -v'
alias dmesg 'dmesg -s 65536'
alias du 'du -s -h'
alias em 'emacs $* -g 90x55'
alias free 'free -m'
alias hosts 'cat /etc/hosts'
alias hpss 'echo Try using psi instead of hpss'
alias l. 'ls --color -aFl'
alias ldmesg 'dmesg -s 65536 | less'
alias less '/usr/bin/less -r'
alias ll 'ls --color -Fl'
alias ls 'ls --color -F'
alias lt 'ls --color -Flt'
alias lt. 'ls --color -aFlt'
alias mdstat 'cat /proc/mdstat'
alias meminfo 'cat /proc/meminfo'
alias mroe 'more'
alias print 'lp'
alias resettermtitle 'echo -ne "\033]0;${myhost}\007"'
alias rtt 'echo -ne "\033]0;${myhost}\007"'
alias sdk 'export DISPLAY 128.165.87.170:0.0;echo The DISPLAY value is now: $DISPLAY'
alias sdl 'DISPLAY 127.0.0.1:0.0;echo The DISPLAY value is now: $DISPLAY'
alias vi 'vim'
alias watchioblocks 'ps -eo stat,pid,user,command | egrep "^STAT|^D|^R"'
alias which 'alias | /usr/bin/which --tty-only --read-alias --show-dot --show-tilde'
alias wmdstat 'watch -n 2 "cat /proc/mdstat"'
alias xload 'xload -label ${myhost}'
