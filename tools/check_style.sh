#!/bin/bash
# -*- Mode: sh -*-

# ------------------------------------------------------------------------------------------------ #
# File  : tools/check_style.sh
# Date  : Tuesday, May 31, 2016, 14:48 pm
# Author: Kelly Thompson <kgt@lanl.gov>
# Note  : Copyright (C) 2016-2021, Triad National Security, LLC., All rights are reserved.
#
# Runs various lint programs in the current directory and list locally modified files that are not
# compliant with the current coding standard (see .clang_format in the top level source directory.)
# - clang-format for C and C++ code.
# - emacs for Fortran90
# - cmake-format and cmake-tidy for CMake scripts.
# ------------------------------------------------------------------------------------------------ #

#--------------------------------------------------------------------------------------------------#
# Environment
#--------------------------------------------------------------------------------------------------#

# Enable job control
set -m

# protect temp files
umask 0077

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  rscriptdir="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink
  # file was located
  [[ $SOURCE != /* ]] && SOURCE="$rscriptdir/$SOURCE"
done
rscriptdir="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"

if [[ -f "${rscriptdir}/common.sh" ]]; then
  # shellcheck source=tools/common.sh
  source "${rscriptdir}/common.sh"
else
  echo " "
  echo "FATAL ERROR: Unable to locate Draco's bash functions: "
  echo "   looking for .../tools/common.sh"
  echo "   searched rscriptdir = $rscriptdir"
  exit 1
fi
export rscriptdir

#--------------------------------------------------------------------------------------------------#
# Support functions
#--------------------------------------------------------------------------------------------------#
print_use()
{
    echo " "
    echo "Usage: ${0##*/} -f -t"
    echo " "
    echo "All arguments are optional."
    echo "  -f Show diff and fix files (when possible)."
    echo -n "  -t Run as a pre-commit check, print list of non-conformant files and return with"
    echo "     exit code = 1 (implies -d)."
    echo " "
}

#--------------------------------------------------------------------------------------------------#
# Sanity Checks
#--------------------------------------------------------------------------------------------------#

# clang-format must be in the PATH
if [[ ${CLANG_FORMAT_VER} ]]; then
  cfver="-${CLANG_FORMAT_VER}"
else
  cfver=""
fi
# Assume applications have version postfix.
gcf=$(which "git-clang-format${cfver}")
cf=$(which "clang-format${cfver}")
# if not found, try to find applications w/o version postfix.
if ! [[ -f "${gcf}" ]]; then
  gcf=$(which git-clang-format)
fi
if ! [[ -f "${cf}" ]]; then
  gcf=$(which clang-format)
fi
# if still not found, abort.
if [[ ! "${gcf}" ]]; then
   echo "ERROR: git-clang-format${cfver} was not found in your PATH."
   echo "pwd="
   pwd
   echo "which git-clang-format${cfver}"
   echo "$gcf"
   exit 1
#else
#  echo "Using $gcf --binary $cf"
fi
if [[ ! "${cf}" ]]; then
   echo "ERROR: clang-format${cfver} was not found in your PATH."
   echo "pwd="
   pwd
   echo "which clang-format${cfver}"
   echo "$cf"
   echo "which git"
   which git
   exit 1
fi

#--------------------------------------------------------------------------------------------------#
# Default values
#--------------------------------------------------------------------------------------------------#
fix_mode=0     # show diffs AND modify code
foundissues=0  # 0 == ok

#--------------------------------------------------------------------------------------------------#
# Command options
#--------------------------------------------------------------------------------------------------#

while getopts ":fht" opt; do
case $opt in
f) fix_mode=1 ;; # also modify code (as possible)
h) print_use; exit 0 ;;
t) fix_mode=0 ;;
\?) echo "" ;echo "invalid option: -$OPTARG"; print_use; exit 1 ;;
:)  echo "" ;echo "option -$OPTARG requires an argument."; print_use; exit 1 ;;
esac
done

#--------------------------------------------------------------------------------------------------#
# Test C++ code with git-clang-format
#--------------------------------------------------------------------------------------------------#

echo -ne "\n--------------------------------------------------------------------------------\n"
echo -ne "Checking modified C/C++ code for style conformance...\n\n"

patchfile_c=$(mktemp /tmp/gcf.patch.XXXXXXXX)

# don't actually modify the files (originally we compared to branch 'develop', but let's try
# ORIG_HEAD or maybe use CI variables like TRAVIS_BRANCH or CI_MERGE_REQUEST_TARGET_BRANCH_NAME).
#run "git branch -a"
#echo "TRAVIS_BRANCH = $TRAVIS_BRANCH"
#echo "CI_MERGE_REQUEST_TARGET_BRANCH_NAME = $CI_MERGE_REQUEST_TARGET_BRANCH_NAME"
target_branch=develop
if [[ -n ${TRAVIS_BRANCH} ]]; then
  target_branch=${TRAVIS_BRANCH}
elif [[ -n ${CI_MERGE_REQUEST_TARGET_BRANCH_NAME} ]]; then
  target_branch=${CI_MERGE_REQUEST_TARGET_BRANCH_NAME}
fi
echo "Looking at code changes compared to target branch = $target_branch"
cmd="${gcf} --binary ${cf} -f --diff --extensions hh,cc,cu $target_branch"
run "${cmd}" &> "${patchfile_c}"

# if the patch file has the string "no modified files to format", the check passes.
if [[ $(grep -c "no modified files" "${patchfile_c}") != 0 ]] || \
   [[ $(grep -c "clang-format did not modify any files" "${patchfile_c}") != 0 ]]; then
  echo "PASS: Changes to C++ sources conform to this project's style requirements."
else
  foundissues=1
  echo "FAIL: some C++ files do not conform to this project's style requirements:"
  # Modify files, if requested
  if [[ ${fix_mode} == 1 ]]; then
    echo -e "      The following patch has been applied to your file.\n"
    run "git apply $patchfile_c"
    cat "${patchfile_c}"
  else
    echo -ne "      run clang-formt -i < file> to fix files, or"
    echo -e  " run ${0##*/} with option -f to automatically apply this patch.\n"
    cat "${patchfile_c}"
  fi
fi
rm -f "${patchfile_c}"

# ------------------------------------------------------------------------------------------------ #
# Liast of modified files
# - Used for cmake-format and Emacs F90 processing.
# ------------------------------------------------------------------------------------------------ #

# staged files
modifiedfiles=$(git diff --name-only --cached)
# unstaged files
modifiedfiles="${modifiedfiles} $(git diff --name-only)"
# all files
modifiedfiles=$(echo "${modifiedfiles}" | sort -u)

#--------------------------------------------------------------------------------------------------#
# Test CMake script for formatting issues
#--------------------------------------------------------------------------------------------------#

CMF=$(which cmake-format)
if [[ $CMF ]]; then
  CMFVER=$("$CMF" --version)
  if [[ $(version_gt "0.5.0" "${CMFVER}") ]]; then
    echo "WARNING: Your version of cmake-format is too old. Expecting v 0.5+. Pre-commit-hook"
    echo "         partially disabled (cmake-format, cmake-lint)"
    unset CMF
  fi
  if ! [[ -f "${rscriptdir}/../.cmake-format.py" ]]; then unset CMF; fi
fi

DIFFVER=$(diff --version | head -n 1 | sed -e 's/.*[ ]\([0-9.]\)/\1/')
[[ $(version_gt "3.4" "${DIFFVER}") ]] && DIFFCOLOR="--color"

if [[ -x $CMF ]]; then

  echo -ne "\n--------------------------------------------------------------------------------\n"
  echo -ne "Checking modified CMake code for style conformance...\n\n"

  patchfile_cmf=$(mktemp /tmp/cmf.patch.XXXXXXXX)

  # file types to parse.
  FILE_EXTS=".cmake"
  FILE_ENDINGS_INCLUDE="CMakeLists.txt"
  # FILE_ENDINGS_EXLCUDE=".cmake.in"

  # Loop over all modified cmake files.  Create one patch containing all changes to these files
  for file in $modifiedfiles; do

    # ignore file if we do check for file extensions and the file does not match any of the
    # extensions specified in $FILE_EXTS
    if ! matches_extension "$file"; then continue; fi

    file_nameonly=$(basename "${file}")
    tmpfile1="/tmp/cmf-${file_nameonly}"
    cp -f "${file}" "${tmpfile1}"
    $CMF -c "${rscriptdir}/../.cmake-format.py" -i "${tmpfile1}" &> /dev/null
    # color output is possible if diff -version >= 3.4 with option `--color`
    diff ${DIFFCOLOR} -u "${file}" "${tmpfile1}" | \
      sed -e "1s|--- |--- a/|" -e "2s|+++ ${tmpfile1}|+++ b/${file}|" >> "$patchfile_cmf"
    rm "${tmpfile1}"

  done

  unset FILE_EXTS
  unset FILE_ENDINGS_INCLUDE
  unset FILE_ENDINGS_EXCLUDE

  # If the patch file is size 0, then no changes are needed.
  if [[ -s "$patchfile_cmf" ]]; then
    foundissues=1
    echo "FAIL: some CMake files do not conform to this project's style requirements:"
    # Modify files, if requested
    if [[ ${fix_mode} == 1 ]]; then
      echo -e "      The following patch has been applied to your file.\n"
      run "git apply $patchfile_cmf"
      cat "${patchfile_cmf}"
    else
      echo -e "      run cmake-format -i <file> to fix files, or"
      echo -e "      run ${0##*/} with option -f to automatically apply this patch.\n"
      cat "${patchfile_cmf}"
    fi
  else
    echo "PASS: Changes to CMake sources conform to this project's style requirements."
  fi
  rm -f "${patchfile_cmf}"

fi

#--------------------------------------------------------------------------------------------------#
# Test CMake script for lint issues
#--------------------------------------------------------------------------------------------------#

CML=$(which cmake-lint)
if [[ $CML ]]; then
  CMLVER=$("$CML" --version)
  if [[ $(version_gt "0.5.0" "${CMLVER}") ]]; then
    echo "WARNING: Your version of cmake-lint is too old. Expecting v 0.5+. Pre-commit-hook"
    echo "         partially disabled (cmake-format, cmake-lint)"
    unset CML
  fi
  if ! [[ -f "${rscriptdir}/../.cmake-format.py" ]]; then unset CML; fi
fi

if [[ -x $CML ]]; then

  echo -ne "\n--------------------------------------------------------------------------------\n"
  echo -ne "Checking modified CMake code for lint conformance...\n\n"

  # file types to parse.
  FILE_EXTS=".cmake"
  FILE_ENDINGS_INCLUDE="CMakeLists.txt"
  # FILE_ENDINGS_EXLCUDE=".cmake.in"

  cml_issues=0
  # Loop over all modified cmake files.  Create one patch containing all changes to these files
  for file in $modifiedfiles; do

    # ignore file if we do check for file extensions and the file does not match any of the
    # extensions specified in $FILE_EXTS
    if ! matches_extension "$file"; then continue; fi

    printf "==> cmake-lint %s\n" "$file"
    $CML --suppress-decoration "${file}" && echo -ne "==> cmake-lint ${file} ... OK\n" || cml_issues=1

  done

  unset FILE_EXTS
  unset FILE_ENDINGS_INCLUDE
  unset FILE_ENDINGS_EXCLUDE

  # If the patch file is size 0, then no changes are needed.
  if [[ ${cml_issues} -gt 0 ]]; then
    foundissues=1
    echo -e "FAIL: some CMake files do not conform to this project's style requirements:"
    echo -e "      You must fix these issues manuall. Run cmake-lint <file> to view issues again."
  else
    echo  "PASS: Changes to CMake sources conform to this project's style requirements."
  fi

fi

#--------------------------------------------------------------------------------------------------#
# Check mode (Test F90 code indentation with emacs and bash)
#--------------------------------------------------------------------------------------------------#

# Defaults ----------------------------------------
EMACS=$(which emacs)
if [[ $EMACS ]]; then
  EMACSVER=$("$EMACS" --version | head -n 1 | sed -e 's/.*Emacs //')
  if [[ $(version_gt "24.0.0" "${EMACSVER}") ]]; then
    echo "WARNING: Your version of emacs is too old. Expecting v 24.0+. Pre-commit-hook partially"
    echo "         disabled (f90 indentation)"
    unset EMACS
  fi
fi

if [[ -x "$EMACS" ]]; then

  echo -ne "\n--------------------------------------------------------------------------------\n"
  echo -e "Checking modified F90 code for style conformance (indentation).\n"

  patchfile_f90=$(mktemp /tmp/emf90.patch.XXXXXXXX)

  # file types to parse.
  FILE_EXTS=".f90 .F90 .f .F"
  # FILE_ENDINGS_INCLUDE="CMakeLists.txt"
  # FILE_ENDINGS_EXLCUDE=".cmake.in"

  # Loop over all modified F90 files.  Create one patch containing all changes to these files
  for file in $modifiedfiles; do

    # ignore file if we do check for file extensions and the file does not match any of the
    # extensions specified in $FILE_EXTS
    if ! matches_extension "$file"; then continue; fi

    file_nameonly=$(basename "${file}")
    tmpfile1="/tmp/f90-format-${file_nameonly}"
    cp -f "${file}" "${tmpfile1}"
    "$EMACS" -batch "${tmpfile1}" -l "${rscriptdir}/../environment/git/f90-format.el" \
      -f emacs-format-f90-sources &> /dev/null
    # color output is possible if diff -version >= 3.4 with option `--color`
    diff ${DIFFCOLOR} -u "${file}" "${tmpfile1}" | \
      sed -e "1s|--- |--- a/|" -e "2s|+++ ${tmpfile1}|+++ b/${file}|" >> "$patchfile_f90"
    rm "${tmpfile1}"

  done

  # If the patch file is size 0, then no changes are needed.
  if [[ -s "$patchfile_f90" ]]; then
    foundissues=1
    echo -n "FAIL: some F90 files do not conform to this project's style requirements:"
    # Modify files, if requested
    if [[ "${fix_mode}" == 1 ]]; then
      echo -e "      The following patch has been applied to your file.\n"
      run "git apply $patchfile_f90"
      cat "$patchfile_f90"
    else
      echo -ne "      run ${0##*/} with option -f to automatically apply this patch.\n"
      cat "$patchfile_f90"
    fi
  else
    echo -n "PASS: Changes to F90 sources conform to this project's style requirements."
  fi
  rm -f "${patchfile_f90}"

fi

#--------------------------------------------------------------------------------------------------#
# Check mode (Test F90 code line length with bash)
#--------------------------------------------------------------------------------------------------#

echo -ne "\n--------------------------------------------------------------------------------"
echo -e "\nChecking modified F90 code for style conformance (line length)..\n"

tmpfile2=$(mktemp /tmp/f90-format-line-len.XXXXXXXX)

# Loop over all modified F90 files.  Create one patch containing all changes to these files
for file in $modifiedfiles; do

  # ignore file if we do check for file extensions and the file does not match any of the extensions
  # specified in $FILE_EXTS
  if ! matches_extension "$file"; then continue; fi

  header_printed=0
  lineno=0

  # shellcheck disable=SC2162
  while read line; do
    (( lineno++ ))
    # Exceptions:
    # - Long URLs
    exception=$(echo "${line}" | grep -i -c http)
    if [[ ${#line} -gt 100 && ${exception} == 0 ]]; then
      if [[ ${header_printed} == 0 ]]; then
        {
          echo -e "\nFile: ${file} [code line too long]\n";
          echo "  line   length content";
          echo -n "  ------ ------ -------------------------------------------------------------";
          echo "-------------------"
        } >> "${tmpfile2}"
        header_printed=1
      fi
      printf "  %-6s %-6s %s\n" "${lineno}" "${#line}" "${line}" >> "${tmpfile2}"
    fi
    # reset exception flag
    exception=0
  done < "${file}"
done

unset FILE_EXTS
unset FILE_ENDINGS_INCLUDE
unset FILE_ENDINGS_EXCLUDE

# If there are issues, report them
if [[ $(wc -l < "${tmpfile2}") -gt 0 ]]; then
    foundissues=1
    echo -ne "FAIL: some F90 files do not conform to this project's style requirements:\n"
    cat "$tmpfile2"
    echo -ne "\nPlease reformat lines listed above to fit into 80 columns and attempt running\n"
    echo -ne "${0##*/} again. These issues cannot be fixed with the -f option."
fi

#--------------------------------------------------------------------------------------------------#
# Done
#--------------------------------------------------------------------------------------------------#

# Return code: 0==ok,1==bad
exit $foundissues

#--------------------------------------------------------------------------------------------------#
# End check_style.sh
#--------------------------------------------------------------------------------------------------#
