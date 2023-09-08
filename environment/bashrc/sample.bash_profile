# .bash_profile
#
# ~/.bash_profile is sourced only for login shells.
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
# CCS-2 standard setup
#
# - If this is an interactive shell then the environment variable $- should
#   contain an "i":
#--------------------------------------------------------------------------------------------------#
case ${-} in
  *i*) export INTERACTIVE=true  ;;
  *)   export INTERACTIVE=false ;;
esac

#--------------------------------------------------------------------------------------------------#
# Draco Developer Environment
#--------------------------------------------------------------------------------------------------#
export DRACO_ENV_DIR=${HOME}/draco/environment
if [[ -f ${DRACO_ENV_DIR}/bashrc/.bashrc ]]; then
  source ${DRACO_ENV_DIR}/bashrc/.bashrc
fi

#--------------------------------------------------------------------------------------------------#
# User Customizations
#--------------------------------------------------------------------------------------------------#
# module unload cmake ('module' may not be available for remote ssh connections)
# if [[ $INTERACTIVE == true ]]; then

  # User Customization goes here

  # export USERNAME=<moniker>
  # export NAME="First Last"
  # export EDITOR="emacs -nw"
  # export LPDEST=gumibears
  # export NO_AT_BRIDGE=1          # Silence warnings from GTK/Gnome
  # export LIBGL_ALWAYS_INDIRECT=1 # Fixes issues with FastX and OpenGL

  # aliases and functions
  # source ~/.bash_aliases
  # source ~/.bash_functions

  # Prompt ----------------------------------------------------------------------#
  # - see http://www.tldp.org/HOWTO/Bash-Prompt-HOWTO/

  # if test "$TERM" = emacs || \
  #   test "$TERM" = dumb  || \
  #   test -z "`declare -f npwd | grep npwd`"; then
  #   export PS1="\h:\w [\!] % "
  #   export LS_COLORS=''
  # else
  #   found=`declare -f npwd | wc -l`
  #   found2=`declare -f parse_git_branch | wc -l`
  #   if test ${found} != 0 && test ${found2} != 0; then
  #     export PS1="\[\033[34m\]\h:\[\033[32m\]\$(npwd)\[\033[35m\]\$(parse_git_branch)\[\033[00m\] [\!] % "
  #   fi
  # fi

  # Modules ------------------------------------------------------------ #

  # target="`uname -n | sed -e s/[.].*//`"
  # case ${target} in
  #   ro-rfe* | ro-login*)
  #     module swap intel intel/15.0.3
  #     ;;
  #   *)
  #     module swap cmake cmake/3.4.0
  #     module list
  #     ;;
  # esac

  # LaTeX
  # add_to_path $mydir TEXINPUTS
  # add_to_path $mydir BSTINPUTS
  # add_to_path $mydir BIBINPUTS

# fi

#--------------------------------------------------------------------------------------------------#
# end ~/.bash_profile
#--------------------------------------------------------------------------------------------------#
