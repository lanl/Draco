# .bashrc
#
# ~/.bashrc is sourced only for non-login shells.
#-------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------#
# CAI-2 standard setup
#
# - If this is an interactive shell then the environment variable $- should
#   contain an "i":
#--------------------------------------------------------------------------------------------------#
case ${-} in
*i*)
   export INTERACTIVE=true
   # If this is an interactive shell and DRACO_ENV_DIR isn't set. Assume that we
   # need to source the .bash_profile.
   if test -z "${DRACO_ENV_DIR}" && test -f ${HOME}/.bash_profile; then
     source $HOME/.bash_profile
   fi
   ;;
*) # Not an interactive shell
  export INTERACTIVE=false
  export DRACO_ENV_DIR=$HOME/draco/environment
  ;;
esac

#--------------------------------------------------------------------------------------------------#
# Draco developer environment
#--------------------------------------------------------------------------------------------------#
if [[ -f ${DRACO_ENV_DIR}/bashrc/.bashrc ]]; then
  # Don't autoload the modules.  Allow me to run 'dracoenv' or 'rmdracoenv'
  # later (or 'rdde').
  export DRACO_ENV_LOAD=OFF
  source ${DRACO_ENV_DIR}/bashrc/.bashrc
fi

#--------------------------------------------------------------------------------------------------#
# User customizations
#--------------------------------------------------------------------------------------------------#
#if test "$INTERACTIVE" = true; then

  # aliases, bash functions
  # source ~/sample.bash_aliases
  # source ~/sample.bash_functions

  # Set terminal title
  # echo -ne "\033]0;${nodename}\007"

  # other personalized settings

  # case `uname -n | sed -e s/[.].*//` in
  # ccscs[123456789])
  #    export CDPATH=.:/home/$USER:/scratch:/scratch/$USER:/scratch/vendors/Modules
  #    module load global exuberant-ctags
  #    ;;
  # esac

# fi

#--------------------------------------------------------------------------------------------------#
# End of ~/.bashrc
#--------------------------------------------------------------------------------------------------#
