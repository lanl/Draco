#!/bin/bash

# Install git pre-commit hooks. The root directory of the target local git repository is expected as
# a parameter.

# Example:
# .../draco/environment/git/install-hooks.sh $HOME/ccs-codes/jayenne

# This file is part of a set of unofficial pre-commit hooks available at github.
# Link:    https://github.com/githubbrowser/Pre-commit-hooks
# Contact: David Martin, david.martin.mailbox@googlemail.com

###########################################################
# CONFIGURATION:
# select which pre-commit hooks are going to be installed
HOOKS="pre-commit pre-commit-clang-format pre-commit-autopep8 pre-commit-flake8 pre-commit-fprettify"
HOOKS="$HOOKS pre-commit-cmake-format pre-commit-cmake-lint"
TOOLS="common.sh"
###########################################################
# There should be no need to change anything below this line.

# exit on error
set -e

SCRIPTPATH="$(dirname -- "$0")"

# shellcheck source=environment/git/canonicalize_filename.sh
source "$SCRIPTPATH/canonicalize_filename.sh"

TOOLSPATH="$(canonicalize_filename "$SCRIPTPATH/../../tools")"

# If argument is provided assume it is the top level of a git repository.
if [[ $1 ]]; then
  dotgitdir="$(canonicalize_filename "$1/.git")"
else
  dotgitdir="$(canonicalize_filename "$SCRIPTPATH/../../.git")"
fi
export dotgitdir
if ! [[ -d "$dotgitdir" ]] ; then
  echo "Error: $dotgitdir does not exist."
  echo "Are you sure $dotgitdir is the root directory of your local git repository?"
  echo "Example: ./install-hooks $HOME/draco"
  exit 1
fi
if ! [[ -d "$dotgitdir/hooks" ]]; then
  run "mkdir -p $dotgitdir/hooks"
fi

# copy hooks to the directory specified as parameter
copy_hooks() {
  for hook in $HOOKS; do
    echo "Copying $hook to $dotgitdir/hooks/."
    cp -i -- "$SCRIPTPATH/$hook" "$dotgitdir/hooks/." || true
  done
  echo "Copying common.sh to $dotgitdir/hooks/."
  cp -i -- "$TOOLSPATH/common.sh" "$dotgitdir/hooks/." || true

  for hook in $HOOKS $TOOLS; do
    chmod +x "$dotgitdir/hooks/$hook"
  done
}

echo ""
echo "Git pre-commit hook installation."
echo ""

# create hooks subfolder if it does not yet exist
if [ -d "$dotgitdir/hooks" ] ; then
  mkdir -p -- "$dotgitdir/hooks"
fi

echo "Copying canonicalize_filename.sh to $dotgitdir/hooks/."
cp -i -- "$SCRIPTPATH/canonicalize_filename.sh" "$dotgitdir/hooks" || true
copy_hooks
echo -ne "\nFinished installation.\n"

# ------------------------------------------------------------------------------------------------ #
# End install-hooks.sh
# ------------------------------------------------------------------------------------------------ #
