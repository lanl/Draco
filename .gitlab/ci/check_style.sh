#!/bin/bash
# -*- Mode: sh -*-

# ------------------------------------------------------------------------------------------------ #
# File  : .gitlab/ci/check_style.sh
# Date  : Thursday, May 20, 2021, 16:31 pm
# Author: Kelly Thompson <kgt@lanl.gov>
# Note  : Copyright (C) 2020-2023, Triad National Security, LLC., All rights are reserved.
#
# Runs various lint programs in the current directory and list locally modified files that are not
# compliant with the current coding standard (see .clang_format in the top level source directory.)
# - clang-format for C and C++ code. (.clang-format)
# - fprettify for Fortran90 (.fprettify.rc)
# - cmake-format and cmake-tidy for CMake scripts (.cmake-format.py)
# ------------------------------------------------------------------------------------------------ #

#--------------------------------------------------------------------------------------------------#
# Environment
#--------------------------------------------------------------------------------------------------#

# Enable job control
set -m

# protect temp files
umask 0077

# load some common bash functions
rscriptdir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" || exit )
if ! [[ -d $rscriptdir ]]; then
  rscriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
fi
if [[ -f "$rscriptdir/common.sh" ]]; then
  # shellcheck source=.gitlab/ci/common.sh
  source "$rscriptdir/common.sh"
else
  echo " "
  echo "FATAL ERROR: Unable to locate Draco's bash functions: "
  echo "   looking for .../tools/common.sh, searched rscriptdir = $rscriptdir"
  exit 1
fi
export rscriptdir

#--------------------------------------------------------------------------------------------------#
# Support functions
#--------------------------------------------------------------------------------------------------#
print_use()
{
    echo -e "\nUsage: ${0##*/} -f -t\n"
    echo "All arguments are optional."
    echo "  -f Show diff and fix files (when possible)."
    echo -n "  -t Run as a pre-commit check, print list of non-conformant files and return with"
    echo -e "     exit code = 1 (implies -d).\n"
}


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

# common ancestor of feature branch (MR) and develop
echo -ne "\n--------------------------------------------------------------------------------\n"
currentbranch=$(git rev-parse --abbrev-ref HEAD)
comancest=$(git merge-base origin/develop "${currentbranch}")
echo "==> current branch       : ${currentbranch}"
echo -ne "==> common ancestor hash : ${comancest}\n"

# The copyright block check needs the full history
run "git fetch --depth=1000000"

# ------------------------------------------------------------------------------------------------ #
# List of modified files
# ------------------------------------------------------------------------------------------------ #

modifiedfiles=$(git diff "${comancest}" --name-only)
# all files
modifiedfiles=$(echo "${modifiedfiles}" | sort -u)
echo "Files modified in this PR:"
for file in $modifiedfiles; do
  echo "  - $file"
done

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
# Test C++ code with git-clang-format
#--------------------------------------------------------------------------------------------------#

echo -ne "\n--------------------------------------------------------------------------------\n"
echo -ne "Checking modified C/C++ code for style conformance...\n\n"

mkdir -p "/tmp/$USER" || die "Could not create /tmp/$USER"
patchfile_c=$(mktemp "/tmp/$USER/gcf.patch.XXXXXXXX")

# don't actually modify the files (originally we compared to branch 'develop', but let's try
# ORIG_HEAD or maybe use CI variables like TRAVIS_BRANCH or CI_MERGE_REQUEST_TARGET_BRANCH_NAME).
target_branch=$comancest
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
    echo -e  "      run ${0##*/} with option -f to automatically apply this patch.\n"
    cat "${patchfile_c}"
  fi
fi
rm -f "${patchfile_c}"

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
  if ! [[ -f "${rscriptdir}/../../.cmake-format.py" ]]; then
    unset CMF;
    echo "WARNING: Missing .cmake-format.py configuration file, skipping cmake-format checks."
  fi
fi

# DIFFVER=$(diff --version | head -n 1 | sed -e 's/.*[ ]\([0-9.]\)/\1/')
# [[ $(version_gt "3.4" "${DIFFVER}") ]] && DIFFCOLOR="--color"

if [[ -x "${CMF}" ]]; then

  echo -ne "\n--------------------------------------------------------------------------------\n"
  echo -ne "Checking modified CMake code for style conformance...\n\n"

  patchfile_cmf=$(mktemp "/tmp/$USER/cmf.patch.XXXXXXXX")

  # file types to parse.
  FILE_EXTS=".cmake"
  FILE_ENDINGS_INCLUDE="CMakeLists.txt"
  # FILE_ENDINGS_EXCLUDE=".cmake.in"

  # Loop over all modified cmake files.  Create one patch containing all changes to these files
  for file in $modifiedfiles; do

    # ignore file if we do check for file extensions and the file does not match any of the
    # extensions specified in $FILE_EXTS
    if ! matches_extension "$file"; then continue; fi

    # ignore file if it no longer exists
    if ! [[ -f "${file}" ]]; then continue; fi

    file_nameonly=$(basename "${file}")
    tmpfile1="/tmp/$USER/cmf-${file_nameonly}"
    echo "==> cmake-format ${file}"
    cp -f "${file}" "${tmpfile1}"
    $CMF -c "${rscriptdir}/../../.cmake-format.py" -i "${tmpfile1}" &> /dev/null
    # color output is possible if diff -version >= 3.4 with option `--color`

    # shellcheck disable=SC2086
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

else
  echo -ne "\n==> Skipping cmake-format checks\n\n"
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
  if ! [[ -f "${rscriptdir}/../../.cmake-format.py" ]]; then unset CML; fi
fi

if [[ -x $CML ]]; then

  echo -ne "\n--------------------------------------------------------------------------------\n"
  echo -ne "Checking modified CMake code for lint conformance...\n\n"

  # file types to parse.
  FILE_EXTS=".cmake"
  FILE_ENDINGS_INCLUDE="CMakeLists.txt"
  FILE_ENDINGS_EXCLUDE="draco-nightly.cmake"
  export FILE_EXTS FILE_ENDINGS_INCLUDE FILE_ENDINGS_EXCLUDE

  cml_issues=0
  # Loop over all modified cmake files.  Create one patch containing all changes to these files
  for file in $modifiedfiles; do

    # ignore file if we do check for file extensions and the file does not match any of the
    # extensions specified in $FILE_EXTS
    if ! matches_extension "$file"; then continue; fi

    # ignore file if it no longer exists
    if ! [[ -f "${file}" ]]; then continue; fi

    printf "==> cmake-lint %s\n" "$file"
    $CML --suppress-decoration --max-branches 20 "${file}" \
      && echo -ne "==> cmake-lint ${file} ... OK\n" || cml_issues=1

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

else
  echo -ne "\n==> Skipping cmake-lint checks\n\n"
fi

#--------------------------------------------------------------------------------------------------#
# Check mode (Test F90 code indentation with fprettify)
#--------------------------------------------------------------------------------------------------#

# Defaults ----------------------------------------
FPY=$(which fprettify)
if [[ $FPY ]]; then
  FPYVER=$("$FPY" --version | tail -n 1 | sed -e 's/.*fprettify //')
  if [[ $(version_gt "0.3.5" "${FPYVER}") ]]; then
    echo "WARNING: Your version of fprettify is too old. Expecting v 0.3.6+. Pre-commit-hook"
    echo "         for f90 indentation will be disabled."
    unset FPY
  fi
fi

if [[ -x "$FPY" ]]; then

  echo -ne "\n--------------------------------------------------------------------------------\n"
  echo -e "Checking modified F90 code for style conformance.\n"

  patchfile_f90=$(mktemp "/tmp/$USER/fpy.patch.XXXXXXXX")
  lintfile_f90=$(mktemp "/tmp/$USER/fpy.lint.XXXXXXXX")

  # file types to parse.
  FILE_EXTS=".f90 .F90 .f .F"
  # FILE_ENDINGS_INCLUDE="CMakeLists.txt"
  # FILE_ENDINGS_EXCLUDE=".cmake.in"
  export FILE_EXTS

  # Loop over all modified F90 files.  Create one patch containing all changes to these files
  for file in $modifiedfiles; do

    # ignore file if we do check for file extensions and the file does not match any of the
    # extensions specified in $FILE_EXTS
    if ! matches_extension "$file"; then continue; fi

    # ignore file if it no longer exists
    if ! [[ -f "${file}" ]]; then continue; fi

    file_nameonly=$(basename "${file}")
    tmpfile1="/tmp/$USER/f90-format-${file_nameonly}"

    # 1. These issues can be fixed automatically.
    "$FPY" -sS "${file}" &> "${tmpfile1}"
    diff -u "${file}" "${tmpfile1}" | \
      sed -e "1s|--- |--- a/|" -e "2s|+++ ${tmpfile1}|+++ b/${file}|" >> "$patchfile_f90"
    rm "${tmpfile1}"

    # 2. These issues can not be fixed automatically.
    cp "${file}" "${tmpfile1}"
    "$FPY" -S "${tmpfile1}"
    "$FPY" "${tmpfile1}" >> "${lintfile_f90}" 2>&1

  done

  # If the patch file is size 0, then no changes are needed.
  if [[ -s "$patchfile_f90" ]] || [[ -s "$lintfile_f90" ]]; then
    foundissues=1
    echo -ne "FAIL: some F90 files do not conform to this project's style requirements:\n"
    # Modify files, if requested
    if [[ -s "$patchfile_f90" ]]; then
      if [[ "${fix_mode}" == 1 ]]; then
        run "git apply $patchfile_f90"
        echo -ne "\n      Changes have been made to your F90 files to meet style guidelines."
        echo -ne "\n      Please check the updated files and add them to your commit.\n"
      else
        echo -ne "      run ${0##*/} with option -f to automatically apply this patch.\n"
        cat "$patchfile_f90"
      fi
    fi
    if [[ -s "$lintfile_f90" ]]; then

      echo -e "      The following F90 style errors must be fixed.\n"
      cat "${lintfile_f90}"
      echo -ne "\n      Please reformat lines listed above and attempt running\n"
      echo -ne "      ${0##*/} again. These issues cannot be fixed with the -f option."
    fi
  else
    echo -n "PASS: Changes to F90 sources conform to this project's style requirements."
  fi
  rm -f "${patchfile_f90}"
  rm -f "${lintfile_f90}"

fi

# ------------------------------------------------------------------------------------------------ #
# Check copyright block
# ------------------------------------------------------------------------------------------------ #

echo -ne "\n--------------------------------------------------------------------------------\n"
echo -e "Checking modified code for copyright block conformance.\n"

patchfile_cb=$(mktemp "/tmp/$USER/copyright_block.patch.XXXXXXXX")

# file types to parse.
FILE_EXTS=".c .cc .cmake .h .hh .in .f90 .F90 .f .F .py .txt"
#FILE_ENDINGS_INCLUDE="_f.h _f77.h _f90.h"
FILE_ENDINGS_EXCLUDE="ChangeLog Release.cc Release.hh check_style.sh imc_info.cc .in kat_main.h"
FILE_ENDINGS_EXCLUDE+=" rngNxW.h util.h terminal.h util_cpu.h fpe_trap.hh"
export FILE_EXTS FILE_ENDINGS_EXCLUDE

# Loop over all modified files.  Create one patch containing all changes to these files
for file in $modifiedfiles; do

  # ignore file if we do check for file extensions and the file does not match any of the
  # extensions specified in $FILE_EXTS
  if ! matches_extension "$file"; then continue; fi

  # If this PR deletes a file, skip it
  if ! [[ -f "${file}" ]]; then continue; fi
  # Gold files, python test inputs, and python run files don't have a copyright line
  if [[ "${file}" =~ "golds/" || "${file}" =~ "mcgrid/test/" || "${file}" =~ "ncpd/test/" \
        || "${file}" =~ "run/" ]]; then
    continue;
  fi

  printf "==> Checking copyright block for %s\n" "$file"

  file_nameonly=$(basename "${file}")
  tmpfile1="/tmp/$USER/copyright-${file_nameonly}"

  # Copy the file and attempt update it.
  cp "${file}" "${tmpfile1}"

  today=$(date +%Y)

  # This data was found in the header comments.  It might be a single year or a range.
  crl=$(grep Copyright "${tmpfile1}")
  # shellcheck disable=SC2001
  create_date=$(echo "${crl}" | sed -e 's/.* \([0-9][0-9]*\).*/\1/')

  # These dates are reported by git
  git_last_mod_date=$(git log -1 "${file}" | grep Date | \
                              sed -e 's/.* \([0-9][0-9][0-9][0-9]\).*/\1/')
  git_create_date=$(git log "${file}" | grep Date | tail -n 1 | \
                            sed -e 's/.* \([0-9][0-9][0-9][0-9]\).*/\1/')

  # Sanity Checks
  if [[ "${create_date}" =~ "Copyright" ]]; then
    echo "Failed to parse copyright line (err 1)"
    exit 1
  fi
  # [[ "${mod_date}" =~ "Copyright" ]] && echo "Failed to parse copyright line" && exit 1
  if [[ "${git_last_mod_date}" =~ "Copyright" ]]; then
    echo "Failed to parse copyright line (err 2)"
    exit 1
  fi
  if [[ "${git_create_date}" =~ "Copyright" ]]; then
    echo "Failed to parse copyright line (err 3)"
    exit 1
  fi
  if [[ "${create_date}" -gt "${today}" ]] || [[ "${create_date}" -lt "1990" ]]; then
    echo "Existing copyright date range is corrupt. Please fix $file manually."
    echo "Create date > today or Create date < 1990."
    exit 1
  fi
  if [[ "${git_create_date}" -gt "${today}" ]] || [[ "${git_create_date}" -lt "1990" ]]; then
    echo "Existing copyright date range is corrupt. Please fix $file manually."
    echo "Git create date > today or Git create date < 1990."
    exit 1
  fi

  # We converted from CVS to svn in 2010. This is the oldest create date that git will report.  In
  # this case older data is lost, so just use whatever is in the file as the create date.
  if [[ "${git_create_date}" -lt "2011" ]] && [[ "${create_date}" -lt "${git_create_date}" ]]; then
    git_create_date="${create_date}"
  fi

  # Expected Copyright line:
  ecrl="Copyright (C) "
  if [[ "${git_create_date}" != "${today}" ]]; then
    ecrl+="${git_create_date}-"
  fi
  ecrl+="${today} Triad National Security, LLC., All rights reserved."

  # If existing copyright spans two lines, reduce it to one line.
  twolines=$(grep -A 1 Copyright "${file}" | tail -n 1 | grep -c reserved)
  twolines_closes_cpp_comment=$(grep -A 1 Copyright "${file}" | tail -n 1 | grep -c '[*]/')
  if [[ $twolines -gt 0 ]]; then
    if [[ $twolines_closes_cpp_comment -gt 0 ]]; then
      sed -i 's%^.*All [Rr]ights [Rr]eserved[.]*.*[*]/$% */%' "${tmpfile1}"
    else
      sed -i '/All rights reserved/d' "${tmpfile1}"
    fi
  fi

  # Do we have terminating comement character on the 'copyright' line.  If so, keep it.
  ecm=""
  if [[ $(echo "${crl}" | grep -c "\\\*/") -gt 0 ]]; then
    ecm=" */"
  fi

  # Replace copyright with new one
  sed -i "s%Copyright.*%${ecrl}${ecm}%" "${tmpfile1}"
  diff -u "${file}" "${tmpfile1}" | \
    sed -e "1s|--- |--- a/|" -e "2s|+++ ${tmpfile1}|+++ b/${file}|" >> "$patchfile_cb"
  rm "${tmpfile1}"

  unset today
  unset crl
  unset create_date
  unset git_last_mod_date
  unset git_create_date
  unset ecrl
  unset twolines
  unset twolines_closes_cpp_comment
  unset ecm

done

# If the patch file is size 0, then no changes are needed.
if [[ -s "$patchfile_cb" ]]; then
  foundissues=1
  # Modify files, if requested
  if [[ "${fix_mode}" == 1 ]]; then
    run "git apply $patchfile_cb"
    echo -ne "\n      Changes have been made to your files to meet Copyright block guidelines."
    echo -ne "\n      Please check the updated files and add them to your commit.\n"
  else
    echo -ne "FAIL: some files do not conform to this project's Copyright block requirements:\n"
    echo -ne "      run ${0##*/} with option -f to automatically apply this patch.\n"
    cat "$patchfile_cb"
  fi
else
  echo -n "PASS: Changes to sources conform to this project's Copyright block requirements."
fi
rm -f "${patchfile_cb}"

#--------------------------------------------------------------------------------------------------#
# Done
#--------------------------------------------------------------------------------------------------#

echo -ne "\n--------------------------------------------------------------------------------\n"
echo "Done"
echo -ne "\n--------------------------------------------------------------------------------\n"

# Return code: 0==ok,1==bad
exit $foundissues

#--------------------------------------------------------------------------------------------------#
# End check_style.sh
#--------------------------------------------------------------------------------------------------#
