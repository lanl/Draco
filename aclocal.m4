# generated automatically by aclocal 1.7.3 -*- Autoconf -*-

# Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002
# Free Software Foundation, Inc.
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

dnl-------------------------------------------------------------------------dnl
dnl ac_conf.m4
dnl
dnl Service macros used in configure.ac's throughout Draco.
dnl
dnl Thomas M. Evans
dnl 1999/02/04 01:56:19
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl AC_DRACO_PREREQ
dnl
dnl Checks the configure version
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DRACO_PREREQ], [dnl

   # we need at least autoconf 2.53 to work correctly
   AC_PREREQ(2.53)

])

dnl-------------------------------------------------------------------------dnl
dnl AC_NEEDS_LIBS
dnl
dnl add DRACO-dependent libraries necessary for a package
dnl usage: configure.ac
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_NEEDS_LIBS, [dnl
   if test ${has_libdir:=no} != "yes" ; then
       DRACO_LIBS="${DRACO_LIBS} -L\${libdir}"
       has_libdir="yes"
   fi

   for lib in $1
   do
       # temporary string to keep line from getting too long
       draco_depends="\${libdir}/lib\${LIB_PREFIX}${lib}\${libsuffix}"
       DRACO_DEPENDS="${DRACO_DEPENDS} ${draco_depends}"
       DRACO_LIBS="${DRACO_LIBS} -l\${LIB_PREFIX}${lib}"
   done

   # Keep a list of component dependencies free of other tags or paths.
   DEPENDENT_COMPONENTS="$1"

])

dnl-------------------------------------------------------------------------dnl
dnl AC_NEEDS_LIBS_TEST
dnl
dnl add DRACO-dependent libraries necessary for a package test
dnl usage: configure.ac
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_NEEDS_LIBS_TEST, [dnl
   DRACO_TEST_LIBS="${DRACO_TEST_LIBS} -L\${libdir}"
   for lib in $1
   do
       # temporary string to keep line from getting too long
       draco_test_depends="\${libdir}/lib\${LIB_PREFIX}${lib}\${libsuffix}"
       DRACO_TEST_DEPENDS="${DRACO_TEST_DEPENDS} ${draco_test_depends}"
       DRACO_TEST_LIBS="${DRACO_TEST_LIBS} -l\${LIB_PREFIX}${lib}"
   done
])

dnl-------------------------------------------------------------------------dnl
dnl AC_RUNTESTS
dnl
dnl add DRACO-package tests (default to use DejaGnu)
dnl usage: in configure.ac:
dnl AC_RUNTESTS(testexec1 testexec2 ... , {nprocs1 nprocs2 ... | scalar})
dnl where serial means run as serial test only.
dnl If compiling with scalar c4 then nprocs are ignored.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_RUNTESTS, [dnl
	test_alltarget="$test_alltarget $1"
        
	test_nprocs="$2"

	if test -z "${test_nprocs}" ; then
	    AC_MSG_ERROR("No procs choosen for the tests!")
        fi
])

dnl-------------------------------------------------------------------------dnl
dnl AC_TESTEXE
dnl
dnl determines what type of executable the tests are, for example, you 
dnl can set the executable to some scripting extension, like python.
dnl the default is an executable binary
dnl options are PYTHON
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_TESTEXE, [dnl
   test_exe="$1"
])

dnl-------------------------------------------------------------------------dnl
dnl AC_INSTALL_EXECUTABLE
dnl
dnl where executables will be installed
dnl usage: configure.ac
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_INSTALL_EXECUTABLE, [ dnl
   install_executable="\${bindir}/\${package}"
   installdirs="${installdirs} \${bindir}"
   alltarget="${alltarget} bin/\${package}"
])

dnl-------------------------------------------------------------------------dnl
dnl AC_INSTALL_LIB
dnl
dnl where libraries will be installed
dnl usage: configure.ac
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_INSTALL_LIB, [ dnl
   install_lib='${libdir}/lib${LIB_PREFIX}${package}${libsuffix}'
   installdirs="${installdirs} \${libdir}"
   alltarget="${alltarget} lib\${LIB_PREFIX}\${package}\${libsuffix}"

   # test will need to link this library
   PKG_DEPENDS='../lib${LIB_PREFIX}${package}${libsuffix}'
   PKG_LIBS='-L.. -l${LIB_PREFIX}${package}'
])

dnl-------------------------------------------------------------------------dnl
dnl AC_INSTALL_HEADERS
dnl
dnl where headers will be installed 
dnl usage: configure.ac
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_INSTALL_HEADERS, [ dnl
   install_headers="\${installheaders}"
   installdirs="${installdirs} \${includedir} \${includedir}/\${package}"
])

dnl-------------------------------------------------------------------------dnl
dnl AC_CHECK_TOOLS
dnl
dnl Find tools used by the build system (latex, bibtex, python, etc)
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DRACO_CHECK_TOOLS], [dnl

   dnl
   dnl TOOL CHECKS
   dnl
   
   dnl check for and assign the path to python
   AC_PATH_PROG(PYTHON_PATH, python, null)
   if test "${PYTHON_PATH}" = null ; then
       AC_MSG_ERROR("No valid Python found!")
   fi
   
   dnl check for and assign the path to perl
   AC_PATH_PROG(PERL_PATH, perl, null)
   if test "${PERL_PATH}" = null ; then
       AC_MSG_WARN("No valid Perl found!")
   fi

   dnl check for CVS
   AC_PATH_PROG(CVS_PATH, cvs, null)
   if test "${CVS_PATH}" = null ; then
       AC_MSG_WARN("No valid CVS found!")
   fi

   dnl check for and assign the path to ghostview
   AC_CHECK_PROGS(GHOSTVIEW, ghostview gv, null)
   if test "${GHOSTVIEW}" = null ; then
       AC_MSG_WARN("No valid ghostview found!")
   fi

   dnl check for and assign the path to latex
   AC_CHECK_PROGS(LATEX, latex, null)
   if test "${LATEX}" = null ; then
       AC_MSG_WARN("No valid latex found!")
   fi
   AC_SUBST(LATEXFLAGS)

   dnl check for and assign the path to bibtex
   AC_CHECK_PROGS(BIBTEX, bibtex, null)
   if test "${BIBTEX}" = null ; then
       AC_MSG_WARN("No valid bibtex found!")
   fi
   AC_SUBST(BIBTEXFLAGS)

   dnl check for and assign the path to xdvi
   AC_CHECK_PROGS(XDVI, xdvi, null)
   if test "${XDVI}" = null ; then
       AC_MSG_WARN("No valid xdvi found!")
   fi
   AC_SUBST(XDVIFLAGS)

   dnl check for and assign the path to dvips
   AC_CHECK_PROGS(DVIPS, dvips, null)
   if test "${DVIPS}" = null ; then
       AC_MSG_WARN("No valid dvips found!")
   fi
   AC_SUBST(DVIPSFLAGS)

   dnl check for and assign the path for printing (lp)
   AC_CHECK_PROGS(LP, lp lpr, null)
   if test "${LP}" = null ; then
       AC_MSG_WARN("No valid lp or lpr found!")
   fi
   AC_SUBST(LPFLAGS)

   dnl check for and assign the path for doxygen
   AC_PATH_PROG(DOXYGEN_PATH, doxygen, null)
   if test "${DOXYGEN_PATH}" = null ; then
       AC_MSG_WARN("No valid Doxygen found!")
   fi
   AC_SUBST(DOXYGEN_PATH)

])

dnl-------------------------------------------------------------------------dnl
dnl AC_ASCI_WHITE_TEST_WORK_AROUND_PREPEND
dnl
dnl changes compiler from newmpxlC to newxlC so that tests can be run
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_ASCI_WHITE_TEST_WORK_AROUND_PREPEND], [dnl

   # change compiler
   if test "${CXX}" = newmpxlC; then
       white_compiler='newmpxlC'
       CXX='newxlC'
       AC_MSG_WARN("Changing to ${CXX} compiler for configure tests.")
   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_ASCI_WHITE_TEST_WORK_AROUND_APPEND
dnl
dnl changes compiler back to newmpxlC
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_ASCI_WHITE_TEST_WORK_AROUND_APPEND], [dnl

   # change compiler back
   if test "${white_compiler}" = newmpxlC; then
       CXX='newmpxlC'
       AC_MSG_WARN("Changing back to ${CXX} compiler.")
   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_HEAD_MAKEFILE
dnl 
dnl Builds default makefile in the head directory
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_HEAD_MAKEFILE], [dnl

   AC_FIND_TOP_SRC($srcdir, package_top_srcdir)
   AC_DBS_VAR_SUBSTITUTIONS
   AC_CONFIG_FILES([Makefile:config/Makefile.head.in])

])

dnl-------------------------------------------------------------------------dnl
dnl AC_SRC_MAKEFILE
dnl 
dnl Builds default makefile in the src directory
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_SRC_MAKEFILE], [dnl

   AC_FIND_TOP_SRC($srcdir, package_top_srcdir)
   AC_DBS_VAR_SUBSTITUTIONS
   AC_CONFIG_FILES([Makefile:../config/Makefile.src.in])

])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_conf.m4
dnl-------------------------------------------------------------------------dnl


dnl-------------------------------------------------------------------------dnl
dnl ac_local.m4
dnl
dnl Macros used internally within the Draco build system.
dnl
dnl Thomas M. Evans
dnl 1999/02/04 01:56:22
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl AC_WITH_DIR
dnl
dnl Define --with-xxx[=DIR] with defaults to an environment variable.
dnl       Usage: AC_WITH_DIR(flag, CPPtoken, DefaultValue, HelpStr)
dnl                for environment variables enter \${ENVIRONVAR} for
dnl                DefaultValue
dnl usage: in aclocal.m4
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_WITH_DIR, [dnl

 dnl
 dnl  The following M4 macros will be expanded into the body of AC_ARG_WITH
 dnl
 dnl AC_PACKAGE is the flag with all dashes turned to underscores
 dnl AC_WITH_PACKAGE will be substituted to the autoconf shell variable
 dnl    with_xxx
 dnl AC_CMDLINE is the shell command to strip double and trailing slashes
 dnl    from directory names.

 define([AC_PACKAGE], [translit($1, [-], [_])])dnl
 define([AC_WITH_PACKAGE], [with_]AC_PACKAGE)dnl
 define([AC_CMDLINE],dnl
[echo "$]AC_WITH_PACKAGE[" | sed 's%//*%/%g' | sed 's%/$%%'])dnl

 AC_ARG_WITH($1,
   [  --with-$1[=DIR]    $4 ($3 by default)],
   if test $AC_WITH_PACKAGE != "no" ; then
      if test $AC_WITH_PACKAGE = "yes" ; then
         # following eval needed to remove possible '\' from $3
         eval AC_WITH_PACKAGE=$3
      fi

      # this command removes double slashes and any trailing slash

      AC_WITH_PACKAGE=`eval AC_CMDLINE`
      if test "$AC_WITH_PACKAGE:-null}" = "null" ; then
         { echo "configure: error: --with-$1 directory is unset" 1>&2; \
           exit 1; }
      fi
      if test ! -d $AC_WITH_PACKAGE ; then
         { echo "configure: error: $AC_WITH_PACKAGE: invalid directory" 1>&2; \
           exit 1; }
      fi

      # this sets up the shell variable, with the name of the CPPtoken,
      # and that we later will do an AC_SUBST on.
      $2="${AC_WITH_PACKAGE}/"

      # this defines the CPP macro with the directory and single slash appended.
      AC_DEFINE_UNQUOTED($2, ${AC_WITH_PACKAGE}/)dnl

      # print a message to the users (that can be turned off with --silent)

      echo "$2 has been set to $$2" 1>&6

   fi)

   AC_SUBST($2)dnl

])
	
dnl-------------------------------------------------------------------------dnl
dnl AC_VENDORLIB_SETUP(1,2)
dnl
dnl set up for VENDOR_LIBS or VENDOR_TEST_LIBS
dnl usage: in aclocal.m4
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_VENDORLIB_SETUP, [dnl

   # $1 is the vendor_<> tag (equals pkg or test)
   # $2 are the directories added 

   if test "${$1}" = pkg ; then
       VENDOR_LIBS="${VENDOR_LIBS} $2"
   elif test "${$1}" = test ; then
       VENDOR_TEST_LIBS="${VENDOR_TEST_LIBS} $2"
   fi
])

dnl-------------------------------------------------------------------------dnl
dnl AC_FIND_TOP_SRC(1,2)
dnl 
dnl Find the top source directory of the package by searching upward
dnl from the argument directory. The top source directory is defined
dnl as the one with a 'config' sub-directory.
dnl
dnl Note: This function will eventually quit if the searched for
dnl directory is not above the argument. It does so when $temp_dir
dnl ceases to be a valid directory, which only seems to happen after a
dnl LOT of ..'s are added to it.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_FIND_TOP_SRC, [dnl
   
   # $1 is the component's source directory
   # $2 is the variable to store the package's main source directory in.

   temp_dir=$1
   AC_MSG_CHECKING([package top source directory])
   while test -d $temp_dir -a ! -d $temp_dir/config ; do   
       temp_dir="${temp_dir}/.."
   done
   if test -d $temp_dir; then
       $2=`cd $temp_dir; pwd;`
       AC_MSG_RESULT([$$2])
   else
       AC_MSG_ERROR('Could not find package top source directory')
   fi
])

dnl-------------------------------------------------------------------------dnl
dnl DO VARIABLE SUBSTITUTIONS ON AC_OUTPUT
dnl
dnl These are all the variable substitutions used within the draco
dnl build system
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_VAR_SUBSTITUTIONS], [dnl

   # these variables are declared "precious", meaning that they are
   # automatically substituted, put in the configure --help, and
   # cached 
   AC_ARG_VAR(CC)dnl
   AC_ARG_VAR(CFLAGS)dnl

   AC_ARG_VAR(CXX)dnl
   AC_ARG_VAR(CXXFLAGS)dnl

   AC_ARG_VAR(LD)dnl
   AC_ARG_VAR(LDFLAGS)dnl

   AC_ARG_VAR(AR)dnl
   AC_ARG_VAR(ARFLAGS)dnl

   AC_ARG_VAR(CPPFLAGS)dnl

   # dependency rules
   AC_SUBST(DEPENDENCY_RULES)

   # other compiler substitutions
   AC_SUBST(STRICTFLAG)dnl
   AC_SUBST(PARALLEL_FLAG)dnl
   AC_SUBST(RPATH)dnl
   AC_SUBST(LIB_PREFIX)dnl

   # install program
   AC_SUBST(INSTALL)dnl
   AC_SUBST(INSTALL_DATA)dnl

   # files to install
   : ${installfiles:='${install_executable} ${install_lib} ${install_headers}'}
   AC_SUBST(installfiles)dnl
   AC_SUBST(install_executable)dnl
   AC_SUBST(install_lib)dnl
   AC_SUBST(install_headers)dnl
   AC_SUBST(installdirs)dnl

   # package libraries
   AC_SUBST(alltarget)dnl
   AC_SUBST(libsuffix)dnl
   AC_SUBST(dirstoclean)dnl
   AC_SUBST(package)dnl
   AC_SUBST(DRACO_DEPENDS)dnl
   AC_SUBST(DRACO_LIBS)dnl
   AC_SUBST(VENDOR_DEPENDS)dnl
   AC_SUBST(VENDOR_INC)dnl
   AC_SUBST(VENDOR_LIBS)dnl
   AC_SUBST(ARLIBS)dnl

   # package testing libraries
   AC_SUBST(PKG_DEPENDS)dnl
   AC_SUBST(PKG_LIBS)dnl
   AC_SUBST(DRACO_TEST_DEPENDS)dnl
   AC_SUBST(DRACO_TEST_LIBS)dnl
   AC_SUBST(VENDOR_TEST_DEPENDS)dnl
   AC_SUBST(VENDOR_TEST_LIBS)dnl
   AC_SUBST(ARTESTLIBS)dnl
   AC_SUBST(test_alltarget)dnl
   AC_SUBST(test_flags)dnl
   AC_SUBST(test_scalar)dnl
   AC_SUBST(test_nprocs)dnl
   AC_SUBST(test_output_files)dnl

   # libraries
   AC_ARG_VAR(LIBS)dnl

   # configure options
   AC_SUBST(configure_command)dnl

   # directories in source tree
   AC_SUBST(package_top_srcdir)
   
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_local.m4
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl ac_vendors.m4
dnl
dnl Macros for each vendor that is used supported by the Draco build
dnl system.
dnl
dnl Thomas M. Evans
dnl 1999/02/04 01:56:22
dnl-------------------------------------------------------------------------dnl


dnl-------------------------------------------------------------------------dnl
dnl AC_MPI_SETUP
dnl
dnl MPI implementation (off by default)
dnl MPI is an optional vendor
dnl
dnl we wait to set the basic MPI libraries (if it is on) until
dnl after checking the C4 status; these functions are performed
dnl in ac_dracoenv.m4, section SYSTEM-SPECIFIC SETUP; we do this
dnl here because each platform has different mpi options for
dnl vendors and mpich
dnl
dnl note that we used to do this in a function called AC_COMM_SET;
dnl however, there are too many platform-dependent variables 
dnl to continue doing this; so we do all these operations in the
dnl platform specific section of ac_dracoenv.m4
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_MPI_SETUP], [dnl

   dnl define --with-mpi
   AC_ARG_WITH(mpi,
      [  --with-mpi=[vendor,mpich,lampi] 
	                  determine MPI implementation (vendor on SGI,SUN; mpich on LINUX)])

   dnl define --with-mpi-inc and --with-mpi-lib
   AC_WITH_DIR(mpi-inc, MPI_INC, \${MPI_INC_DIR},
	       [tell where MPI includes are])
   AC_WITH_DIR(mpi-lib, MPI_LIB, \${MPI_LIB_DIR},
	       [tell where MPI libs are])

   # determine if this package is needed for testing or for the
   # package
   vendor_mpi=$1 

   # set default value for with_mpi which is no
   if test "${with_mpi:=no}" = yes ; then 
       with_mpi='vendor'
   fi

   # if the user sets MPI_INC and MPI_LIB directories then turn on  
   # with_mpi and set it to vendor if with_mpi=no to begin with
   if test "${with_mpi}" = no ; then
       if test -n "${MPI_INC}" ; then
	   with_mpi='vendor'
       elif test -n "${MPI_LIB}" ; then
	   with_mpi='vendor'
       fi
   fi
   
   # if c4=mpi and with-mpi=no explicitly then 
   # define them (mpi gets set to vendor by default)
   if test "$with_c4" = mpi ; then
       if test "$with_mpi" = no ; then
	   with_mpi='vendor'
       fi
   fi

]) 


AC_DEFUN([AC_MPI_FINALIZE], [dnl

   # only add stuff if mpi is not no and the vendor is defined
   if test "${with_mpi}" != no && test -n "${vendor_mpi}"; then

       # include path
       if test -n "${MPI_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${MPI_INC}"
       fi
   
       # libraries
       if test -n "${MPI_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_mpi, -L${MPI_LIB} ${mpi_libs})
       elif test -z "${MPI_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_mpi, ${mpi_libs})
       fi

       # add MPI directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${MPI_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${MPI_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_SPRNG_SETUP
dnl
dnl SPRNG LIBRARY SETUP (on by default -lfg)
dnl SPRNG is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_SPRNG_SETUP], [dnl

   dnl define --with-sprng
   AC_ARG_WITH(sprng,
      [  --with-sprng[=lib]      determine the rng lib (lfg is default)])
	
   dnl define --with-sprng-inc and --with-sprng-lib
   AC_WITH_DIR(sprng-inc, SPRNG_INC, \${SPRNG_INC_DIR},
	       [tell where SPRNG includes are])
   AC_WITH_DIR(sprng-lib, SPRNG_LIB, \${SPRNG_LIB_DIR},
	       [tell where SPRNG libraries are])

   # determine if this package is needed for testing or for the 
   # package
   vendor_sprng=$1

   # choices are with_sprng = lfg, lcg, yes, or no

   # default (sprng is set to lfg by default)
   if test "${with_sprng:=lfg}" = yes ; then
       with_sprng='lfg'
   fi

])


AC_DEFUN([AC_SPRNG_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_sprng}"; then

       # include path
       if test -n "${SPRNG_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${SPRNG_INC}"
       fi
   
       # libraries
       if test -n "${SPRNG_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_sprng, -L${SPRNG_LIB} -l${with_sprng})
       elif test -z "${SPRNG_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_sprng, -l${with_sprng})
       fi

       # add sprng directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${SPRNG_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${SPRNG_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_AZTEC_SETUP
dnl
dnl AZTEC SETUP (on by default)
dnl AZTEC is a required vendor
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_AZTEC_SETUP], [dnl

   dnl define --with-aztec
   AC_ARG_WITH(aztec,
      [  --with-aztec=[lib]      determine the aztec lib (aztec is the default)])
 
   dnl define --with-aztec-inc
   AC_WITH_DIR(aztec-inc, AZTEC_INC, \${AZTEC_INC_DIR},
	       [tell where AZTEC includes are])

   dnl define --with-aztec-lib
   AC_WITH_DIR(aztec-lib, AZTEC_LIB, \${AZTEC_LIB_DIR},
	       [tell where AZTEC libraries are])

   # set default value of aztec includes and libs
   if test "${with_aztec:=aztec}" = yes ; then
       with_aztec='aztec'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_aztec=$1

])


AC_DEFUN([AC_AZTEC_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_aztec}" ; then

       # include path
       if test -n "${AZTEC_INC}"; then 
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${AZTEC_INC}"
       fi

       # library path
       if test -n "${AZTEC_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_aztec, -L${AZTEC_LIB} -l${with_aztec})
       elif test -z "${AZTEC_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_aztec, -l${with_aztec})
       fi

       # add AZTEC directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${AZTEC_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${AZTEC_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_GSL_SETUP
dnl
dnl GSL SETUP (on by default)
dnl GSL is a required vendor
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_GSL_SETUP, [dnl

   dnl define --with-gsl
   AC_ARG_WITH(gsl,
      [  --with-gsl=[gsl] 
                       determine GSL lib (gsl is default)])
 
   dnl define --with-gsl-inc
   AC_WITH_DIR(gsl-inc, GSL_INC, \${GSL_INC_DIR},
	       [tell where GSL includes are])

   dnl define --with-gsl-lib
   AC_WITH_DIR(gsl-lib, GSL_LIB, \${GSL_LIB_DIR},
	       [tell where GSL libraries are])

   # set default value of gsl includes and libs
   if test "${with_gsl:=gsl}" = yes ; then
       with_gsl='gsl'
   fi

   # if atlas is available use it's version of cblas, 
   # otherwise use the version provided by GSL
   if test "${with_lapack}" = atlas; then
       gsl_libs='-lgsl'
   else
       gsl_libs='-lgsl -lgslcblas'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_gsl=$1
])


AC_DEFUN([AC_GSL_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_gsl}"; then

       # include path
       if test -n "${GSL_INC}"; then 
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${GSL_INC}"
       fi

       # library path
       if test -n "${GSL_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gsl, -L${GSL_LIB} ${gsl_libs})
       elif test -z "${GSL_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gsl, ${gsl_libs})
       fi

       # add GSL directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${GSL_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${GSL_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_TRILINOS_SETUP
dnl
dnl TRILINOS SETUP (on by default)
dnl TRILINOS is a required vendor
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_TRILINOS_SETUP], [dnl

   dnl define --with-trilinos
   AC_ARG_WITH(trilinos,
      [  --with-trilinos=[lib]    determine the trilinos implementation (aztecoo is default)])
 
   dnl define --with-trilinos-inc
   AC_WITH_DIR(trilinos-inc, TRILINOS_INC, \${TRILINOS_INC_DIR},
	       [tell where TRILINOS includes are])

   dnl define --with-trilinos-lib
   AC_WITH_DIR(trilinos-lib, TRILINOS_LIB, \${TRILINOS_LIB_DIR},
	       [tell where TRILINOS libraries are])

   # set default value of trilinos includes and libs
   if test "${with_trilinos:=aztecoo}" = yes ; then
       with_trilinos='aztecoo'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_trilinos=$1

])


AC_DEFUN([AC_TRILINOS_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_trilinos}" ; then

       # include path
       if test -n "${TRILINOS_INC}"; then 
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${TRILINOS_INC}"
       fi

       # library path
       if test -n "${TRILINOS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_trilinos, -L${TRILINOS_LIB} -l${with_trilinos} -lepetra -ltriutils)
       elif test -z "${TRILINOS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_trilinos, -l${with_trilinos} -lepetra -ltriutils)
       fi

       # add TRILINOS directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${TRILINOS_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${TRILINOS_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_SCALAPACK_SETUP
dnl
dnl SCALAPACK SETUP 
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_SCALAPACK_SETUP], [dnl

   dnl define --with-scalapack
   AC_ARG_WITH(scalapack,
      [  --with-scalapack=[scalapack] ])
 
   dnl define --with-scalapack-lib
   AC_WITH_DIR(scalapack-lib, SCALAPACK_LIB, \${SCALAPACK_LIB_DIR},
	       [tell where SCALAPACK libraries are])

   # set default value of scalapack includes and libs
   if test "${with_scalapack:=scalapack}" = yes ; then
       with_scalapack='scalapack'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_scalapack=$1

])


AC_DEFUN([AC_SCALAPACK_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_scalapack}" ; then

       # no includes for scalapack

       # library path
       if test -n "${SCALAPACK_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_scalapack, -L${SCALAPACK_LIB} -lscalapack)
       elif test -z "${SCALAPACK_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_scalapack, -lscalapack)
       fi

       # add SCALAPACK directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${SCALAPACK_LIB}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_BLACS_SETUP
dnl
dnl BLACS SETUP 
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_BLACS_SETUP], [dnl

   dnl define --with-blacs
   AC_ARG_WITH(blacs,
      [  --with-blacs=[blacs] ])
 
   dnl define --with-blacs-lib
   AC_WITH_DIR(blacs-lib, BLACS_LIB, \${BLACS_LIB_DIR},
	       [tell where BLACS libraries are])

   # set default value of blacs includes and libs
   if test "${with_blacs:=blacs}" = yes ; then
       with_blacs='blacs'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_blacs=$1

])


AC_DEFUN([AC_BLACS_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_blacs}" ; then

       # no includes for blacs

       # library path
       if test -n "${BLACS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_blacs, -L${BLACS_LIB} -lblacsF77init -lblacsCinit -lblacs -lblacsCinit -lblacs)
       elif test -z "${BLACS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_blacs, -lblacsF77init -lblacsCinit -lblacs -lblacsCinit -lblacs)
       fi

       # add BLACS directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${BLACS_LIB}"

   fi

])
dnl-------------------------------------------------------------------------dnl
dnl AC_HYPRE_SETUP
dnl
dnl HYPRE SETUP 
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_HYPRE_SETUP], [dnl

   dnl define --with-hypre
   AC_ARG_WITH(hypre,
      [  --with-hypre=[hypre] ])
 
   dnl define --with-hypre-inc
   AC_WITH_DIR(hypre-inc, HYPRE_INC, \${HYPRE_INC_DIR},
	       [tell where HYPRE includes are])

   dnl define --with-hypre-lib
   AC_WITH_DIR(hypre-lib, HYPRE_LIB, \${HYPRE_LIB_DIR},
	       [tell where HYPRE libraries are])

   # set default value of hypre includes and libs
   if test "${with_hypre:=hypre}" = yes ; then
       with_hypre='hypre'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_hypre=$1

])


AC_DEFUN([AC_HYPRE_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_hypre}" ; then

       # include path
       if test -n "${HYPRE_INC}"; then 
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${HYPRE_INC}"
       fi

       # library path
       if test -n "${HYPRE_LIB}" ; then

	   AC_VENDORLIB_SETUP(vendor_hypre, -L${HYPRE_LIB} -lHYPRE_parcsr_ls -lHYPRE_DistributedMatrixPilutSolver -lHYPRE_ParaSails -lHYPRE_Euclid -lHYPRE_MatrixMatrix -lHYPRE_DistributedMatrix -lHYPRE_IJ_mv -lHYPRE_parcsr_mv -lHYPRE_seq_mv -lHYPRE_krylov -lHYPRE_utilities)

       elif test -z "${HYPRE_LIB}" ; then

	   AC_VENDORLIB_SETUP(vendor_hypre, -lHYPRE_parcsr_ls -lHYPRE_DistributedMatrixPilutSolver -lHYPRE_ParaSails -lHYPRE_Euclid -lHYPRE_MatrixMatrix -lHYPRE_DistributedMatrix -lHYPRE_IJ_mv -lHYPRE_parcsr_mv -lHYPRE_seq_mv -lHYPRE_krylov -lHYPRE_utilities)

       fi

       # add HYPRE directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${HYPRE_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${HYPRE_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_METIS_SETUP
dnl
dnl METIS SETUP (on by default)
dnl METIS is a required vendor
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_METIS_SETUP], [dnl

   dnl define --with-metis
   AC_ARG_WITH(metis,
      [  --with-metis=[lib]    the metis implementation])
 
   dnl define --with-metis-inc
   AC_WITH_DIR(metis-inc, METIS_INC, \${METIS_INC_DIR},
	       [tell where METIS includes are])

   dnl define --with-metis-lib
   AC_WITH_DIR(metis-lib, METIS_LIB, \${METIS_LIB_DIR},
	       [tell where METIS libraries are])

   # set default value of metis includes and libs
   if test "${with_metis:=metis}" = yes ; then
       with_metis='metis'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_metis=$1

])


AC_DEFUN([AC_METIS_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_metis}" ; then

       # include path
       if test -n "${METIS_INC}"; then 
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${METIS_INC}"
       fi

       # library path
       if test -n "${METIS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_metis, -L${METIS_LIB} -l${with_metis})
       elif test -z "${METIS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_metis, -l${with_metis})
       fi

       # add METIS directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${METIS_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${METIS_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_PCG_SETUP
dnl
dnl PCG LIBRARY SETUP (on by default)
dnl PCG is a required vendor
dnl
dnl note that we add some system-specific libraries for this
dnl vendor in AC_DRACO_ENV; also, the user must set up LAPACK for
dnl this to work
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_PCG_SETUP], [dnl

   dnl define --with-pcg
   AC_ARG_WITH(pcg,        
      [  --with-pcg[=lib]        determine the pcg lib name (pcg is default)])

   dnl define --with-pcg-lib
   AC_WITH_DIR(pcg-lib, PCG_LIB, \${PCG_LIB_DIR},
	       [tell where PCG libraries are])

   # determine if this package is needed for testing or for the 
   # package
   vendor_pcg=$1

   # pcg is set to libpcg by default
   if test "${with_pcg:=pcg}" = yes ; then
       with_pcg='pcg'
   fi

])


AC_DEFUN([AC_PCG_FINALIZE], [dnl

   # set up the libraries
   if test -n "${vendor_pcg}"; then

       # library path
       if test -z "${PCG_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_pcg, -l${with_pcg})
       elif test -n "${PCG_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_pcg, -L${PCG_LIB} -l${with_pcg})
       fi

       # add PCG directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${PCG_LIB}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_GANDOLF_SETUP
dnl
dnl GANDOLF LIBRARY SETUP (on by default)
dnl GANDOLF is a required vendor
dnl
dnl SGI needs "-lfortran" on the load line when including libgandolf.a.
dnl This library is added to ${LIBS} in AC_DRACO_ENV.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_GANDOLF_SETUP], [dnl

   dnl define --with-gandolf
   AC_ARG_WITH(gandolf,        
      [  --with-gandolf[=lib]    determine the gandolf lib name (gandolf is default)])

   dnl define --with-gandolf-lib
   AC_WITH_DIR(gandolf-lib, GANDOLF_LIB, \${GANDOLF_LIB_DIR},
	       [tell where GANDOLF libraries are])

   # determine if this package is needed for testing or for the 
   # package
   vendor_gandolf=$1

   # gandolf is set to libgandolf by default
   if test "${with_gandolf:=gandolf}" = yes ; then
       with_gandolf='gandolf'
   fi

])


AC_DEFUN([AC_GANDOLF_FINALIZE], [dnl

   # set up the libraries
   if test -n "${vendor_gandolf}"; then

       # set up library paths
       if test -z "${GANDOLF_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gandolf, -l${with_gandolf})
       elif test -n "${GANDOLF_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gandolf, -L${GANDOLF_LIB} -l${with_gandolf})
       fi

       # add GANDOLF directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${GANDOLF_LIB}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_EOSPAC5_SETUP
dnl
dnl EOSPAC5 LIBRARY SETUP (on by default)
dnl EOSPAC5 is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_EOSPAC5_SETUP], [dnl

   dnl define --with-eospac
   AC_ARG_WITH(eospac,        
      [  --with-eospac[=lib]     determine the eospac lib name (eospac is default)])

   dnl define --with-eospac-lib
   AC_WITH_DIR(eospac-lib, EOSPAC5_LIB, \${EOSPAC5_LIB_DIR},
	       [tell where EOSPAC5 libraries are])

   # determine if this package is needed for testing or for the 
   # package (valid values are pkg or test)
   vendor_eospac=$1

   # eospac is set to libeospac by default
   if test "${with_eospac:=eospac}" = yes ; then
       with_eospac='eospac'
   fi

])


AC_DEFUN([AC_EOSPAC5_FINALIZE], [dnl

   # set up the libraries
   if test -n "${vendor_eospac}"; then

       # set up library paths
       if test -z "${EOSPAC5_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_eospac, -l${with_eospac})
       elif test -n "${EOSPAC5_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_eospac, -L${EOSPAC5_LIB} -l${with_eospac})
       fi

       # add EOSPAC5 directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${EOSPAC5_LIB}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_LAPACK_SETUP
dnl
dnl LAPACK SETUP (on by default)
dnl LAPACK is a required vendor
dnl
dnl NOTE: this also sets up the BLAS
dnl
dnl note that we add system specific libraries to this list in
dnl ac_platforms.m4
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_LAPACK_SETUP], [dnl

   dnl define --with-lapack
   AC_ARG_WITH(lapack,
      [  --with-lapack=[vendor,atlas]
                          determine LAPACK implementation (vendor default)])

   dnl define --with-lapack-lib
   AC_WITH_DIR(lapack-lib, LAPACK_LIB, \${LAPACK_LIB_DIR}, 
	       [tell where LAPACK libs are])

   # determine if this package is needed for testing or for the 
   # package
   vendor_lapack=$1

   # lapack is set to vendor by default
   if test "${with_lapack:=vendor}" = yes ; then
       with_lapack='vendor'
   fi

   # define the atlas libraries (these are system independent)
   if test "${with_lapack}" = atlas; then
       lapack_libs='-llapack -lf77blas -lcblas -latlas'
   fi
])


AC_DEFUN([AC_LAPACK_FINALIZE], [dnl

   # set up lapack libraries
   if test -n "${vendor_lapack}"; then

       # set libraries
       if test -z "${LAPACK_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_lapack, ${lapack_libs})
       elif test -n "${LAPACK_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_lapack, -L${LAPACK_LIB} ${lapack_libs})
       fi

       # add LAPACK directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${LAPACK_LIB}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_GRACE_SETUP
dnl
dnl GRACE SETUP (on by default)
dnl GRACE is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_GRACE_SETUP], [dnl

   dnl define --with-grace
   AC_ARG_WITH(grace,
      [  --with-grace=[lib]      determine the grace lib (grace_np is the default])
 
   dnl define --with-grace-inc
   AC_WITH_DIR(grace-inc, GRACE_INC, \${GRACE_INC_DIR},
	       [tell where GRACE includes are])

   dnl define --with-grace-lib
   AC_WITH_DIR(grace-lib, GRACE_LIB, \${GRACE_LIB_DIR},
	       [tell where GRACE libraries are])

   # set default value of grace includes and libs
   if test "${with_grace:=grace_np}" = yes ; then
       with_grace='grace_np'
   fi

   # define GRACE header file
   GRACE_H="<${with_grace}.h>"
   AC_DEFINE_UNQUOTED(GRACE_H, ${GRACE_H})dnl

   # determine if this package is needed for testing or for the 
   # package
   vendor_grace=$1

])


AC_DEFUN([AC_GRACE_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_grace}" ; then

       # include path
       if test -n "${GRACE_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${GRACE_INC}"
       fi

       # library path
       if test -n "${GRACE_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_grace, -L${GRACE_LIB} -l${with_grace})
       elif test -z "${GRACE_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_grace, -l${with_grace})
       fi

       # add GRACE directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${GRACE_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${GRACE_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_SPICA_SETUP
dnl
dnl SPICA LIBRARY SETUP (on by default -lSpicaCSG)
dnl SPICA is an optional vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_SPICA_SETUP], [dnl

   dnl define --with-spica
   AC_ARG_WITH(spica,
      [  --with-spica[=yes]                 spica is on by default])
	
   dnl define --with-spica-inc and --with-spica-lib
   AC_WITH_DIR(spica-inc, SPICA_INC, \${SPICA_INC_DIR},
	       [tell where SPICA includes are])
   AC_WITH_DIR(spica-lib, SPICA_LIB, \${SPICA_LIB_DIR},
	       [tell where SPICA libraries are])

   # determine if this package is needed for testing or for the 
   # package
   vendor_spica=$1

   # define variable if spica is on
   if test "${with_spica:=yes}" != no; then
       AC_DEFINE([USE_SPICA])
   fi
])


AC_DEFUN([AC_SPICA_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_spica}"; then

       # include path
       if test -n "${SPICA_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${SPICA_INC}"
       fi
   
       # libraries
       if test -n "${SPICA_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_spica, -L${SPICA_LIB} -lSpicaCSG)
       elif test -z "${SPICA_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_spica, -lSpicaCSG)
       fi

       # add spica directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${SPICA_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${SPICA_INC}"

   fi
])

dnl-------------------------------------------------------------------------dnl
dnl AC_XERCES_SETUP
dnl
dnl XERCES LIBRARY SETUP
dnl xerces is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_XERCES_SETUP], [dnl

   dnl define --with-xerces
   AC_ARG_WITH(xerces,
      [  --with-xerces[=lib]      determine the XERCES xml lib (xerces-c is default)])
	
   dnl define --with-xerces-inc and --with-xerces-lib
   AC_WITH_DIR(xerces-inc, XERCES_INC, \${XERCES_INC_DIR},
	       [tell where XERCES includes are])
   AC_WITH_DIR(xerces-lib, XERCES_LIB, \${XERCES_LIB_DIR},
	       [tell where XERCES libraries are])

   # determine if this package is needed for testing or for the 
   # package
   vendor_xerces=$1

   # default (xerces is set to xerces-c by default)
   if test "${with_xerces:=xerces-c}" = yes ; then
       with_xerces='xerces-c'
   fi
])


AC_DEFUN([AC_XERCES_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_xerces}"; then

       # include path
       if test -n "${XERCES_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${XERCES_INC}"
       fi
   
       # libraries
       if test -n "${XERCES_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_xerces, -L${XERCES_LIB} -l${with_xerces})
       elif test -z "${XERCES_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_xerces, -l${with_xerces})
       fi

       # add xerces directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${XERCES_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${XERCES_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_HDF5_SETUP
dnl
dnl HDF5 SETUP (on by default; 'mpi' if mpi in use, else 'serial')
dnl HDF5 is an optional vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_HDF5_SETUP], [dnl

   dnl define --with-hdf5
   AC_ARG_WITH(hdf5,
      [  --with-hdf5=[serial,mpi]      determine hdf5 implementation (default:  'mpi' if mpi in use, else 'serial')])
 
   dnl define --with-hdf5-inc
   AC_WITH_DIR(hdf5-inc, HDF5_INC, \${HDF5_INC_DIR},
	       [tell where HDF5 includes are])

   dnl define --with-hdf5-lib
   AC_WITH_DIR(hdf5-lib, HDF5_LIB, \${HDF5_LIB_DIR},
	       [tell where HDF5 libraries are])

   # default (mpi if mpi is in use, else serial)
   if test "${with_hdf5:=no}" = yes ; then
       if test "${with_mpi}" != no ; then
	   with_hdf5='mpi'
       else
	   with_hdf5='serial'
       fi
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_hdf5=$1

   # define variable if hdf5 is on
   if test "${with_hdf5:=yes}" != no; then
       AC_DEFINE([USE_HDF5])
   fi

])


AC_DEFUN([AC_HDF5_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_hdf5}" ; then

       # include path
       if test -n "${HDF5_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${HDF5_INC}"
       fi

       # library path
       if test -n "${HDF5_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_hdf5, -L${HDF5_LIB} -lhdf5 -lz)
       elif test -z "${HDF5_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_hdf5, -lhdf5 -lz)
       fi

       # add HDF5 directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${HDF5_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${HDF5_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_UDM_SETUP
dnl
dnl UDM SETUP (on by default; 'mpi' if mpi in use, else 'serial')
dnl UDM is an optional vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_UDM_SETUP], [dnl

   dnl define --with-udm
   AC_ARG_WITH(udm,
      [  --with-udm=[serial,mpi]      determine udm implementation (default:  'mpi' if mpi in use, else 'serial')])
 
   dnl define --with-udm-inc
   AC_WITH_DIR(udm-inc, UDM_INC, \${UDM_INC_DIR},
	       [tell where UDM includes are])

   dnl define --with-udm-lib
   AC_WITH_DIR(udm-lib, UDM_LIB, \${UDM_LIB_DIR},
	       [tell where UDM libraries are])

   # default (mpi if mpi is in use, else serial)
   if test "${with_udm:=no}" = yes ; then
       if test "${with_mpi}" != no ; then
	   with_udm='mpi'
       else
	   with_udm='serial'
       fi
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_udm=$1

   # define variable if udm is on
   if test "${with_udm:=no}" != no; then
       AC_DEFINE([USE_UDM])
   fi

])


AC_DEFUN([AC_UDM_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_udm}" ; then

       # include path
       if test -n "${UDM_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${UDM_INC}"
           # set extra #define if using udm in parallel
           if test "${with_udm}" = mpi ; then
               AC_DEFINE(UDM_HAVE_PARALLEL)
           fi
       fi

       # library path
       if test -n "${UDM_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_udm, -L${UDM_LIB} -ludm)
       elif test -z "${UDM_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_udm, -ludm)
       fi

       # add UDM directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${UDM_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${UDM_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_DLOPEN_SETUP
dnl
dnl This is an optional vendor.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DLOPEN_SETUP], [dnl

   dnl define --enable-dlopen
   AC_ARG_ENABLE(dlopen,
      [  --enable-dlopen          Enable dlopen (default: on if --enable-shared, off otherwise)])

   # determine if this package is needed for testing or for the
   # package.
   vendor_dlopen=$1 

   # set default value for enable_dlopen, which is the value of enable_shared.
   if test "${enable_shared}" = yes ; then
       if test "${enable_dlopen:=yes}" != no ; then 
	   enable_dlopen=yes
       fi
   else
       if test "${enable_dlopen:=no}" != no ; then 
	   enable_dlopen=yes
       fi
   fi

   # turn off dlopen if not using shared libraries.
   if test "${enable_shared}" != yes ; then
       if test "${enable_dlopen}" = yes ; then
	   AC_MSG_WARN("Must specify --enable-shared when using --enable-dlopen.")
           AC_MSG_WARN("   dlopen disabled.")
       fi
       enable_dlopen=no
   fi

   if test "${enable_dlopen}" = yes ; then
       AC_DEFINE(USE_DLOPEN)
   fi
]) 


AC_DEFUN([AC_DLOPEN_FINALIZE], [dnl
   # Libraries are platform-specific; done in ac_platforms.
])

dnl-------------------------------------------------------------------------dnl
dnl AC_VENDOR_FINALIZE
dnl
dnl Run at the end of the environment setup to add defines required by
dnl the vendors.  We do this to allow platform specific mods to the 
dnl vendor defines BEFORE they are added to CCPFLAGS, etc. 
dnl
dnl This macro needs to be updated when new vendors are added.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_VENDOR_FINALIZE], [dnl

   # call finalize functions for each vendor, the order is important
   # each vendor setup is appended to the previous; thus, the calling
   # level goes from high to low
   AC_TRILINOS_FINALIZE
   AC_GSL_FINALIZE

   AC_AZTEC_FINALIZE
   AC_PCG_FINALIZE
   AC_HYPRE_FINALIZE
   AC_SCALAPACK_FINALIZE
   AC_BLACS_FINALIZE
   AC_LAPACK_FINALIZE
   AC_EOSPAC5_FINALIZE
   AC_GANDOLF_FINALIZE
   AC_SPRNG_FINALIZE
   AC_GRACE_FINALIZE
   AC_METIS_FINALIZE
   AC_SPICA_FINALIZE
   AC_XERCES_FINALIZE

   AC_UDM_FINALIZE
   AC_HDF5_FINALIZE

   AC_MPI_FINALIZE
   AC_DLOPEN_FINALIZE

   # print out vendor include paths
   AC_MSG_CHECKING("vendor include paths")
   if test -n "${VENDOR_INC_DIRS}"; then
       AC_MSG_RESULT("${VENDOR_INC_DIRS}")
   else
       AC_MSG_RESULT("no vendor include dirs defined")
   fi

   # print out vendor lib paths
   AC_MSG_CHECKING("vendor lib paths")
   if test -n "${VENDOR_LIB_DIRS}"; then
       AC_MSG_RESULT("${VENDOR_LIB_DIRS}")
   else
       AC_MSG_RESULT("no vendor lib dirs defined")
   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_ALL_VENDORS_SETUP
dnl
dnl DRACO INCLUSIVE VENDOR MACRO
dnl-------------------------------------------------------------------------dnl
dnl allows one to include all vendor macros by calling this macro.
dnl designed for draco/configure.in and draco/src/configure.in

AC_DEFUN(AC_ALL_VENDORS_SETUP, [dnl

   dnl include all macros for easy use in top-level configure.in's
   AC_MPI_SETUP(pkg)
   AC_SPRNG_SETUP(pkg)
   AC_PCG_SETUP(pkg)
   AC_AZTEC_SETUP(pkg)
   AC_GSL_SETUP(pkg)
   AC_TRILINOS_SETUP(pkg)
   AC_METIS_SETUP(pkg)
   AC_LAPACK_SETUP(pkg)
   AC_GANDOLF_SETUP(pkg)
   AC_EOSPAC5_SETUP(pkg)
   AC_GRACE_SETUP(pkg)
   AC_SPICA_SETUP(pkg)
   AC_XERCES_SETUP(pkg)
   AC_HDF5_SETUP(pkg)
   AC_UDM_SETUP(pkg)
   AC_DLOPEN_SETUP(pkg)
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_vendors.m4
dnl-------------------------------------------------------------------------dnl


dnl-------------------------------------------------------------------------dnl
dnl ac_dracoarg.m4
dnl
dnl Declarations of Draco configure options (with some default
dnl settings). 
dnl
dnl Thomas M. Evans
dnl 1999/02/04 01:56:20
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl AC_DRACO_ARGS
dnl
dnl Declaration of Draco non-vendor configure options. This macro can 
dnl be called to fill out configure help screens
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_ARGS, [dnl

   dnl
   dnl Library prefix
   dnl
     
   AC_ARG_WITH(lib-prefix,
      [  --with-lib-prefix[=library prefix]
                          give prefix to libraries (default rtt_)])

   # default for lib_prefix is rtt_
   LIB_PREFIX="${with_lib_prefix:=rtt_}"
   if test "${LIB_PREFIX}" = no ; then
       LIB_PREFIX=''
   fi

   dnl
   dnl c4 toggle (scalar by default)
   dnl

   dnl define --with-c4
   AC_ARG_WITH(c4, 
      [  --with-c4[=scalar,mpi,shmem]   
		          turn on c4 (default scalar) ])

   # give with-c4 implied argument
   if test "${with_c4:=scalar}" = yes ; then
       with_c4='scalar'
   fi

   dnl
   dnl DBC toggle
   dnl

   dnl defines --with-dbc
   AC_ARG_WITH(dbc,
      [  --with-dbc[=level]      set Design-by-Contract])
	
   if test "${with_dbc}" = yes ; then
       with_dbc='7'
   elif test "${with_dbc}" = no ; then
       with_dbc='0'
   fi
	
   dnl
   dnl SHARED versus ARCHIVE libraries
   dnl

   dnl defines --enable-shared
   AC_ARG_ENABLE(shared,
      [  --enable-shared         turn on shared libraries (.a default)])

   dnl
   dnl CHOOSE A C++ COMPILER
   dnl

   dnl defines --with-cxx
   AC_ARG_WITH(cxx,
      [  --with-cxx[=gcc,sgi,kcc,compaq,guide]                                    
                          choose a c++ compiler (defaults are machine dependent)])

   dnl the default is gcc
   if test "${with_cxx}" = yes ; then
       with_cxx='gcc'
   fi

   dnl
   dnl STATIC VERSUS DYNAMIC LINKING
   dnl

   dnl defines --enable-static-ld
   AC_ARG_ENABLE(static-ld,
      [  --enable-static-ld      use (.a) libraries if possible])

   dnl
   dnl ANSI STRICT COMPLIANCE
   dnl

   dnl defines --enable-strict-ansi
   AC_ARG_ENABLE(strict-ansi,
      [  --disable-strict-ansi   turn off strict ansi compliance])

   dnl
   dnl ONE_PER INSTANTIATION FLAG
   dnl

   dnl defines --enable-one-per
   AC_ARG_ENABLE(one-per,
      [  --disable-one-per       turn off --one_per flag])

   dnl
   dnl COMPILER OPTIMZATION LEVEL
   dnl

   dnl defines --with-opt
   AC_ARG_WITH(opt,
      [  --with-opt[=0,1,2,3]    set optimization level (0 by default)])

   if test "${with_opt}" = yes ; then
       with_opt='0'
   fi

   dnl defines --enable-debug
   AC_ARG_ENABLE(debug,
      [  --enable-debug          turn on debug (-g) option])

   dnl
   dnl POSIX SOURCE
   dnl

   dnl defines --with-posix
   AC_ARG_WITH(posix,
      [  --with-posix[=num]      give posix source (system-dependent defaults)])

   dnl
   dnl ADD TO CPPFLAGS
   dnl
   
   dnl defines --with-cppflags
   AC_ARG_WITH(cppflags,
      [  --with-cppflags[=flags] add flags to \$CPPFLAGS])

   dnl
   dnl ADD TO CXXFLAGS
   dnl
   
   dnl defines --with-cxxflags
   AC_ARG_WITH(cxxflags,
      [  --with-cxxflags[=flags] add flags to \$CXXFLAGS])

   dnl
   dnl ADD TO CFLAGS
   dnl
   
   dnl defines --with-cflags
   AC_ARG_WITH(cflags,
      [  --with-cflags[=flags]   add flags to \$CFLAGS])

   dnl
   dnl ADD TO F90FLAGS
   dnl
   
   dnl defines --with-f90flags
   AC_ARG_WITH(f90flags,
      [  --with-f90flags[=flags] add flags to \$F90FLAGS])

   dnl
   dnl ADD TO ARFLAGS
   dnl
   
   dnl defines --with-arflags
   AC_ARG_WITH(arflags,
      [  --with-arflags[=flags]  add flags to \$ARFLAGS])

   dnl
   dnl ADD TO LDFLAGS
   dnl
   
   dnl defines --with-ldflags
   AC_ARG_WITH(ldflags,
      [  --with-ldflags[=flags]  add flags to \$LDFLAGS])

   dnl 
   dnl ADD TO LIBRARIES
   dnl

   dnl defines --with-libs
   AC_ARG_WITH(libs,
      [  --with-libs=[libs]      add libs to \$LIBS])

   dnl
   dnl CHOSE BIT COMPILATION ON SGI'S
   dnl

   dnl defines --enable-32-bit
   AC_ARG_ENABLE(32-bit,
      [  --enable-32-bit         do 32-bit compilation (compiler dependent)])

   dnl defines --enable-64-bit
   AC_ARG_ENABLE(64-bit,
      [  --enable-64-bit         do 64-bit compilation (compiler dependent)])

   dnl
   dnl CHOSE MIPS INSTRUCTION SET ON SGI'S
   dnl

   dnl defines --with-mips
   AC_ARG_WITH(mips,
      [  --with-mips[=1,2,3,4]   set mips, mips4 by default (SGI ONLY)])

   if test "${with_mips}" = yes ; then
       with_mips='4'
   fi

   dnl 
   dnl STLport
   dnl

   dnl specify location of stlport installation.
   AC_ARG_WITH(stlport,
      [  --with-stlport        replace default STL with stlPort (off by default)])

   dnl Doxygen options

   AC_ARG_ENABLE(latex-doc,
      [  --enable-latex-doc      build latex docs with doxygen (off by default)],
      [AC_SUBST(latex_yes_no,'YES')],
      [AC_SUBST(latex_yes_no,'NO')])

   AC_ARG_WITH(doc-output,
      [  --with-doc-output=path  build documentation in path (prefix/documentation by default)],
      [AC_SUBST(doxygen_output_top,${with_doc_output})],
      [doxygen_output_top='DEFAULT'])

   dnl end of AC_DRACO_ARGS
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_dracoarg.m4
dnl-------------------------------------------------------------------------dnl


dnl ========================================================================
dnl 
dnl 	Author:	Mark G. Gray
dnl 		Los Alamos National Laboratory
dnl 	Date:	Wed Apr 19 16:39:19 MDT 2000
dnl 
dnl 	Copyright (c) 2000 U. S. Department of Energy. All rights reserved.
dnl 
dnl ========================================================================

dnl NAME

dnl	AC_WITH_F90, AC_F90_ENV

dnl SYNOPSIS/USAGE

dnl     AC_WITH_F90
dnl     AC_F90_ENV

dnl DESCRIPTION

dnl     AC_WITH_F90 sets the variable with_f90 to yes if it is not already 
dnl     set.

dnl     AC_F90_ENV set environment variables F90, F90FLAGS, F90EXT, 
dnl     F90FREE, F90FIXED, and MODFLAG for the compiler requested by 
dnl     with_f90.  If no specific compiler is requested, guess a compiler 
dnl     based on the target
dnl
========================================================================

dnl ### Ensure with_f90 set
AC_DEFUN(AC_WITH_F90, [dnl
   : ${with_f90:=yes}
    
   dnl turn off C++ compiler
   with_cxx='no'

   dnl defines --with-f90
   AC_ARG_WITH(f90,
       [  --with-f90[=XL,Fujitsu,Lahey,Portland,WorkShop,Cray,MIPS,Compaq,HP,Intel,NAG,Absoft]
                          choose an F90 compiler])
])

dnl
dnl CHOOSE A F90 COMPILER
dnl

AC_DEFUN(AC_F90_ENV, [dnl
   AC_REQUIRE([AC_CANONICAL_HOST])

   case "${with_f90:=yes}" in
   XL)
       AC_COMPILER_XL_F90
   ;;
   Fujitsu)
       AC_COMPILER_FUJITSU_F90
   ;;
   Lahey)
       AC_COMPILER_LAHEY_F90
   ;;
   Portland)
       AC_COMPILER_PORTLAND_F90
   ;;
   WorkShop)
       AC_COMPILER_WORKSHOP_F90
   ;;
   Cray)
      AC_COMPILER_CRAY_F90
   ;;
   MIPS)
       AC_COMPILER_MIPS_F90
   ;;
   Compaq)
       AC_COMPILER_COMPAQ_F90
   ;;
   HP)
       AC_COMPILER_HP_F90
   ;;
   Intel)
       AC_COMPILER_INTEL_F90
   ;;
   NAG)
       AC_COMPILER_NAG_F90
   ;;
   Absoft)
       AC_COMPILER_ABSOFT_F90
   ;;
   yes)				# guess compiler from target platform
       case "${host}" in   
       rs6000-ibm-aix*)
           AC_COMPILER_XL_F90
       ;;
       powerpc-ibm-aix*)
           AC_COMPILER_XL_F90
       ;;
       sparc-sun-solaris2.*)
           AC_COMPILER_WORKSHOP_F90
       ;;
       i?86-pc-linux*)
           AC_COMPILER_LAHEY_F90
       ;;
       ymp-cray-unicos*)
          AC_COMPILER_CRAY_F90
       ;;
       mips-sgi-irix*)
          AC_COMPILER_MIPS_F90
       ;;
       i??86-pc-cygwin*)
          AC_COMPILER_COMPAQ_F90
       ;;
       alpha*)
          AC_COMPILER_COMPAQ_F90
       ;;
       *hp-hpux*)
          AC_COMPILER_HP_F90
       ;;
       *)
          AC_MSG_ERROR([Cannot guess F90 compiler, set --with-f90])
       ;;
       esac
   ;;
   no)
   ;;
   *)
       AC_MSG_ERROR([Unrecognized F90 compiler, use --help])
   ;;
   esac

   AC_SUBST(F90FREE)
   AC_SUBST(F90FIXED)
   AC_SUBST(F90FLAGS)
   AC_SUBST(MODFLAG)
])

dnl-------------------------------------------------------------------------dnl
dnl IBM XLF95 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_XL_F90, [dnl

   # Check for working XL F90 compiler

  if test "${with_upslib:=no}" != "no"
  then
     AC_CHECK_PROG(F90, mpxlf95, mpxlf95, none)
     if test "${F90}" != mpxlf95
     then
         AC_MSG_ERROR([not found])
     fi
  else
     AC_CHECK_PROG(F90, xlf95, xlf95, none)
     if test "${F90}" != xlf95
     then
         AC_MSG_ERROR([not found])
     fi
  fi
  
   # FREE, FIXED AND MODULE FLAGS

   F90FREE='-qfree=f90'
   F90FIXED='-qfixed'
   MODFLAG='-I'

   # LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=

   # COMPILATION FLAGS

   if test "$F90FLAGS" = ""
   then
     # F90FLAGS="-qsuffix=f=f90 -qmaxmem=-1 -qextchk -qarch=pwr2 -bmaxstack:0x70000000 -bmaxdata:0x70000000 -qalias=noaryovrlp -qhalt=s ${F90FREE}"
       F90FLAGS="-qsuffix=f=f90 -qmaxmem=-1 -qextchk -qarch=auto -bmaxstack:0x70000000 -bmaxdata:0x70000000 -qalias=noaryovrlp -qnosave -qlanglvl=95pure -qzerosize ${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	   trapflags="-qinitauto=FF"
	   trapflags="${trapflags} -qflttrap=overflow:underflow:zerodivide:invalid:enable"
	   trapflags="${trapflags} -qsigtrap"
	   F90FLAGS="-g -d -C ${trapflags} -bloadmap:loadmap.dat ${F90FLAGS}"
       else
	 # F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
	   F90FLAGS="-O3 ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_XL_F90
])

dnl-------------------------------------------------------------------------dnl
dnl FUJITSU F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_FUJITSU_F90, [dnl

   # Check for working Fujitsu F90 compiler

   AC_CHECK_PROG(F90, f90, f90, none)
   if test "${F90}" = f90 && ${F90} -V 2>&1 | grep "Fujitsu"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG

   F90FREE='-Free'
   F90FIXED='-Fixed'
   MODFLAG='-I'

   # LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-static-flib'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="-X9 -Am ${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g -Haesu ${F90FLAGS}"
       else
	    F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_FUJITSU_F90
])

dnl-------------------------------------------------------------------------dnl
dnl LAHEY F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_LAHEY_F90, [dnl

   # Check for working Lahey F90 compiler

   AC_CHECK_PROG(F90, lf95, lf95, none)
   if test "${F90}" = lf95 && ${F90} --version 2>&1 | grep "Lahey"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG

   F90FREE='--nfix'
   F90FIXED='--fix'
   MODFLAG='-I'

   # LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-static-flib'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
     # F90FLAGS="--f95 ${F90FREE}"
       F90FLAGS="--staticlink --f95 --in --info --swm 2004,2006,2008,8202,8203,8204,8205,8206,8209,8220 ${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	  # F90FLAGS="-g --chk --trace ${F90FLAGS}"
	    F90FLAGS="-g --ap --chk --pca --private --trap --wo ${F90FLAGS}"
       else
	  # F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
	    F90FLAGS="-O --ntrace ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_LAHEY_F90
])

dnl-------------------------------------------------------------------------dnl
dnl PORTLAND F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_PORTLAND_F90, [dnl

   # Check for working Portland Group F90 compiler

   AC_CHECK_PROG(F90, pgf90, pgf90, none)
   if test "${F90}" = pgf90 && ${F90} --V 2>&1 | grep "Portland"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG

   F90FREE='-Mfreeform'
   F90FIXED='-Mnofreeform'
   MODFLAG='-module'

   # LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC=

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g -Mbounds -Mchkptr ${F90FLAGS}"
       else
	    F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_PORTLAND_F90
])

dnl-------------------------------------------------------------------------dnl
dnl COMPAQ F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_COMPAQ_F90, [dnl

   # Check for working compaq F90 compiler

   AC_CHECK_PROG(F90, f95, f95, none)
   if test "${F90}" = f95 && ${F90} -version 2>&1 | grep "Fortran"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG

   F90FREE=''
   F90FIXED=''
   MODFLAG='-I'

   # LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-non_shared'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
     # F90FLAGS="${F90FREE} -assume byterecl"
       F90FLAGS="${F90FREE} -assume byterecl -automatic -std -warn argument_checking"

       if test "${enable_debug:=no}" = yes
       then
	  # F90FLAGS="-g ${F90FLAGS}"
	    F90FLAGS="-g -check bounds -fpe2 ${F90FLAGS}"
       else
	  # F90FLAGS="-O ${F90FLAGS}"
	    F90FLAGS="-O5 -arch host -assume noaccuracy_sensitive -math_library accurate -tune host ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_COMPAQ_F90
])

dnl-------------------------------------------------------------------------dnl
dnl SUN WORKSHOP F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_WORKSHOP_F90, [dnl

   # Check for working WorkShop F90 compiler

   AC_CHECK_PROG(F90, f90, f90, none)
   if test "${F90}" = f90 && ${F90} -V 2>&1 | grep "WorkShop"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # Set F90FREE, F90FIXED, and MODFLAG

   F90FREE='-free'
   F90FIXED='-fixed'
   MODFLAG='-M'

   # Set LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-Bstatic'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g"
       else
	    F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_WORKSHOP_F90
])

dnl-------------------------------------------------------------------------dnl
dnl CRAY_F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_CRAY_F90, [dnl

   # Check for working Cray F90 compiler

   AC_CHECK_PROG(F90, f90, f90, none)
   if test "${F90}" = f90
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # FREE, FIXED AND MODULE FLAGS

   F90FREE='-f free'
   F90FIXED='-f fixed'
   MODFLAG='-p'

   # LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	   F90FLAGS="-g ${F90FLAGS}"
       else
	   F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_CRAY_F90
])

dnl-------------------------------------------------------------------------dnl
dnl IRIX MIPS F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_MIPS_F90, [dnl

   # Look for working MIPS compiler

   AC_CHECK_PROG(F90, f90, f90, none)
   if test "${F90}" = f90 && ${F90} -version 2>&1 | grep "MIPS"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # Set F90FREE, F90FIXED, and MODFLAG

   F90FREE='-freeform'
   F90FIXED='-col72'
   MODFLAG='-I'

   # LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
	#F90FLAGS="${F90FREE} -OPT:Olimit=0"
	F90FLAGS="${F90FREE} -mips4 -r10000 -DEBUG:fullwarn=ON:woff=878,938,1193,1438"

	if test "${enable_debug:=no}" = yes
	then
	  # F90FLAGS="-g ${F90FLAGS}"
	    F90FLAGS="-g -check_bounds -DEBUG:trap_uninitialized=ON ${F90FLAGS}"
	else
	  # F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
	    F90FLAGS="-O3 -OPT:IEEE_arithmetic=2:roundoff=2 ${F90FLAGS}"
	fi
   fi

   dnl end of AC_COMPILER_MIPS_F90
])

dnl-------------------------------------------------------------------------dnl
dnl HP F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_HP_F90, [dnl

   # CHECK FOR WORKING HP F90 COMPILER
   AC_CHECK_PROG(F90, f90, f90, none)
   if test "${F90}" = f90 && ${F90} +version 2>&1 | grep "HP"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG
   F90FREE='+source=free'
   F90FIXED='+source=fixed'
   MODFLAG='-I'

   # LINKER AND LIBRARY (AR)
   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='+noshared'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="${F90FREE} +U77"

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g -C ${F90FLAGS}"
       else
	    F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_HP_F90
])

dnl-------------------------------------------------------------------------dnl
dnl INTEL F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_INTEL_F90, [dnl

   # CHECK FOR WORKING INTEL F90 COMPILER
   AC_CHECK_PROG(F90, ifc, ifc, none)
   if test "${F90}" = ifc && ${F90} -V 2>&1 | grep "Intel"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG
   F90FREE='-FR'
   F90FIXED='-FI'
   MODFLAG='-I '
   MODSUFFIX='mod'

   # LINKER AND LIBRARY (AR)
   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-static'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="${F90FREE} -e95"

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g -C -implicitnone ${F90FLAGS}"
       else
	    F90FLAGS="-O3 -fno-alias -tpp7 -ipo -pad -align ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_INTEL_F90
])

dnl-------------------------------------------------------------------------dnl
dnl NAG F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_NAG_F90, [dnl

   # CHECK FOR WORKING NAG F90 COMPILER
   AC_CHECK_PROG(F90, f95, f95, none)
   if test "${F90}" = f95 && ${F90} -V 2>&1 | grep "NAGWare"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG
   F90FREE='-free'
   F90FIXED='-fixed'
   MODFLAG='-I '

   # LINKER AND LIBRARY (AR)
   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-unsharedf95'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="${F90FREE} -colour -info -target=native"

       if test "${enable_debug:=no}" = yes
       then
          # only use first line if memory error is suspected, too much output
          #   otherwise
	  # F90FLAGS="-g -C -mtrace=size -nan -u ${F90FLAGS}"
	    F90FLAGS="-g -C -nan -u ${F90FLAGS}"
       else
	    F90FLAGS="-O4 ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_NAG_F90
])

dnl-------------------------------------------------------------------------dnl
dnl ABSOFT F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_ABSOFT_F90, [dnl

   # CHECK FOR WORKING ABSOFT F90 COMPILER
   AC_CHECK_PROG(F90, f95, f95, none)
  
   # F90FREE, F90FIXED AND MODFLAG
   F90FREE=''
   F90FIXED=''
   MODFLAG='-p '

   # LINKER AND LIBRARY (AR)
   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC=''

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="-cpu:host -en"

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g -et -m0 -M399,1193,878 -Rb -Rc -Rs -Rp -trap=ALL ${F90FLAGS}"
       else
	    F90FLAGS="-O3 ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_ABSOFT_F90
])

dnl ========================================================================

