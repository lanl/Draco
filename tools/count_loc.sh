#! /bin/sh
###################################################################
#  NAME
#	count_loc.sh
#
#  SYNOPSIS
#       count_loc.sh [dir]
#
#  DESCRIPTION
#       A script to estimate lines of code.
#       The directory to be scanned is given in [dir]
#       If [dir] is absent, then the current directory is examined.
#       Attempts are made to find all files of certain types and
#       to strip out comments and blank lines, however this
#       should not be considered to be exhaustive.
#
#       Author: Mark G. Gray, Randy M. Roberts, John M. McGhee
#               Los Alamos National Laboratory
#       Date:   Wen Sep 29 10:44:35 MST 1999
#
#
###################################################################
# $Id$
###################################################################

#Count lines of code in a C++ file by counting semicolons
cpp_loc()
{
  grep ";" $* | wc -l
}

#Count lines of code in a script file, stripping any blank lines 
#and/or comments
script_loc()
{
  awk '/^[ \t]*$/ || $1 ~ /#/ {next} {print}' $* | wc -l
}

#Count lines of code in a html file, stripping any blank lines
#and the first line of any comments
html_loc()
{
  awk '/^[ \t]*$/ || $1 ~ /<!--/ {next} {print}' $* | wc -l
}

#Count lines of code in a tex file, stripping any blank lines
#and/or comments
latex_loc()
{
  awk '/^[ \t]*$/ || $1 ~ /%/ {next} {print}' $* | wc -l
}

#Count lines of code in a elisp file, stripping any blank lines
#and/or comments
elisp_loc()
{
  awk '/^[ \t]*$/ || $1 ~ /;;/ {next} {print}' $* | wc -l
}

#Define a fucntion which  lists the names of any 
#shell script files
list_shell_scripts()
{
for file in $*
do 
    if file $file | grep 'shell script' >/dev/null
    then
      echo $file
    fi
done
}

#Define a fucntion which  lists the names of any 
#perl script files
list_perl_scripts()
{
for file in $*
do 
    if file $file | grep 'perl script' >/dev/null
    then
      echo $file
    fi
done
}

#Get the command line arguments (if any)
if [ $# = 0 ]
then
  TMP="."
elif [ $# = 1 ]
then
  TMP="$1"
else
  TMP="$1"
  echo "Warning -- multiple command line arguments found,"
  echo "           using first argument and ignoring others"
  echo
fi

#Check for existence of the directory to be scanned
if [ ! -d $TMP ];
then
 echo $TMP " is not a directory."
 exit 1
fi

#Get the full pathname of the directory
CODE_DIR=`(cd $TMP; pwd)`

# Print the header
echo
echo "Counting lines of source code -- " 
echo "Directory: " $CODE_DIR 
echo "Date     : " `date`
echo

# Go to the directory to be scanned.
cd $CODE_DIR

#Define C++ files 
# (Note: the output from the find command using this filter
#        could be saved into a list to speed things up)

#------------------------- C++ Source Code ------------------------------

#Scan for C++ source code
echo "Total C++ source: "
find . -type f \( -name '*.cc' -o -name '*.hh' \) ! -name '#*' -exec cat {} \; |  cpp_loc
echo

#------------------------ C++ Test Code ----------------------------------

#Scan for C++ test code
echo "C++ source in test directories: "
for i in `find . -name test -type d -print`
do
  find ${i} -type f \( -name '*.cc' -o -name '*.hh' \) ! -name '#*' -exec cat {} \; 
done | cpp_loc
echo

#Scan for C++ Design-by-Contract specifications
echo "C++ contract specifications:"
for i in `find . -type f \( -name '*.cc' -o -name '*.hh' \) ! -name '#*' -print` 
do
  awk '$1 ~ /Assert|Require|Ensure|Check|Insist/ {print}' ${i}
done | cpp_loc
echo

#------------------------ Documentation ----------------------------------

#Scan for C++ comments (finds both C and C++ style comments)
echo "C++ comments:"
for i in `find . -type f \( -name '*.cc' -o -name '*.hh' \) ! -name '#*' -print`
do 
  awk '$1~/\/\// {print} /\/\*/, /\*\// {print}' ${i}
done | wc -l
echo

#Scan for LaTeX source code
echo 'LaTeX Documentation:'
find . -name "*.tex" ! -name '#*' -type f -exec cat {} \; |  latex_loc
echo

#Scan for html source
echo 'html source: '
find . -type f \( -name '*.html' -o -name '*.htm' -o -name '*.shml'  \) ! -name '#*' -exec cat {} \; |  html_loc
echo

#------------------------ Scripts ----------------------------------

#Scan for executable script source code. 
#(such as sh, bash, csh, etc.)
echo Executable shell script source:
files=`find . -type f -perm -100  ! -name 'configure' ! -name '*~' ! -name '#*' -print`
script_loc `list_shell_scripts $files` /dev/null
echo

#Scan for Perl source code. 
echo Executable Perl source:
files=`find . -type f -perm -100  ! -name '*~' ! -name '#*' -print`
script_loc `list_perl_scripts $files` /dev/null
echo

#Scan for Python source code
echo "Python source:"
find . -name '*.py' ! -name '#*' -type f -exec cat {} \; | script_loc
echo

#Scan for Expect source code
echo "Expect source:"
find . -name '*.exp' ! -name '#*' -type f -exec cat {} \; | script_loc
echo

#------------------------ Templates ----------------------------------

#Scan for Elisp source code
echo "Elisp source:"
find . -name '*.el' ! -name '#*' -type f  -exec cat {} \; | elisp_loc
echo


#################### end of count_loc.sh ######################################
