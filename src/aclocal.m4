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
dnl Service macros used in configure.in's throughout Draco.
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
dnl usage: configure.in
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
])

dnl-------------------------------------------------------------------------dnl
dnl AC_NEEDS_LIBS_TEST
dnl
dnl add DRACO-dependent libraries necessary for a package test
dnl usage: configure.in
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
dnl usage: in configure.in:
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
dnl usage: configure.in
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
dnl usage: configure.in
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
dnl usage: configure.in
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
dnl end of ac_conf.m4
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
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_MPI_SETUP, [dnl

   dnl define --with-mpi
   AC_ARG_WITH(mpi,
      [  --with-mpi=[vendor,mpich] 
	                  determine MPI implementation (vendor on SGI,SUN; mpich on LINUX)])

   dnl define --with-mpi-inc and --with-mpi-lib
   AC_WITH_DIR(mpi-inc, MPI_INC, \${MPI_INC_DIR},
	       [tell where MPI includes are])
   AC_WITH_DIR(mpi-lib, MPI_LIB, \${MPI_LIB_DIR},
	       [tell where MPI libs are])

   # define MPI include path
   if test -n "${MPI_INC}" ; then
       # remember that MPI_INC has the final slash
       MPI_H="\"${MPI_INC}mpi.h\""
   elif test -z "${MPI_INC}" ; then
       MPI_H="<mpi.h>"
   fi

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

   # add MPI directory to VENDOR_DIRS
   VENDOR_DIRS="${MPI_LIB} ${VENDOR_DIRS}"

   dnl we wait to set the basic MPI libraries (if it is on) until
   dnl after checking the C4 status; these functions are performed
   dnl in ac_dracoenv.m4, section SYSTEM-SPECIFIC SETUP; we do this
   dnl here because each platform has different mpi options for
   dnl vendors and mpich

   dnl note that we used to do this in a function called AC_COMM_SET;
   dnl however, there are too many platform-dependent variables 
   dnl to continue doing this; so we do all these operations in the
   dnl platform specific section of ac_dracoenv.m4
])

dnl-------------------------------------------------------------------------dnl
dnl AC_SPRNG_SETUP
dnl
dnl SPRNG LIBRARY SETUP (on by default -lfg)
dnl SPRNG is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_SPRNG_SETUP, [dnl

   dnl define --with-sprng
   AC_ARG_WITH(sprng,
      [  --with-sprng[=lib]      determine the rng lib (lfg is default)])
	
   dnl define --with-sprng-inc and --with-sprng-lib
   AC_WITH_DIR(sprng-inc, SPRNG_INC, \${SPRNG_INC_DIR},
	       [tell where SPRNG includes are])
   AC_WITH_DIR(sprng-lib, SPRNG_LIB, \${SPRNG_LIB_DIR},
	       [tell where SPRNG libraries are])

   # define SPRNG include path
   if test -n "${SPRNG_INC}" ; then
       # remember that SPRNG_INC has the final slash
       SPRNG_H="\"${SPRNG_INC}sprng.h\""
   elif test -z "${SPRNG_INC}" ; then
       SPRNG_H="<sprng.h>"
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_sprng=$1

   # choices are with_sprng = lfg, lcg, yes, or no

   # default (sprng is no and set to lfg by default)
   if test "${with_sprng:=lfg}" = yes ; then
       with_sprng='lfg'
   fi

   # set up the libraries
   if test "${with_sprng}" != no ; then
       if test -n "${SPRNG_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_sprng, -L${SPRNG_LIB} -l${with_sprng})
       elif test -z "${SPRNG_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_sprng, -l${with_sprng})
       fi
   fi

   # add sprng directory to VENDOR_DIRS
   VENDOR_DIRS="${SPRNG_LIB} ${VENDOR_DIRS}"
])

dnl-------------------------------------------------------------------------dnl
dnl AC_AZTEC_SETUP
dnl
dnl AZTEC SETUP (on by default)
dnl AZTEC is a required vendor
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_AZTEC_SETUP, [dnl

   dnl define --with-aztec
   AC_ARG_WITH(aztec,
      [  --with-aztec=[lib]      determine the aztec lib (aztec is the default])
 
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

   # define AZTEC include path
   if test -n "${AZTEC_INC}" ; then
       # remember that AZTEC_INC has the final slash
       AZTEC_H="\"${AZTEC_INC}az_aztec.h\""
       AZTEC_DEFS_H="\"${AZTEC_INC}az_aztec_defs.h\""
   elif test -z "${AZTEC_INC}" ; then
       AZTEC_H="<az_aztec.h>"
       AZTEC_DEFS_H="<az_aztec_defs.h>"
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_aztec=$1

   # set up the libraries
   if test "${with_aztec}" != no ; then
       if test -n "${AZTEC_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_aztec, -L${AZTEC_LIB} -l${with_aztec})
       elif test -z "${AZTEC_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_aztec, -l${with_aztec})
       fi
   fi

   # add AZTEC directory to VENDOR_DIRS
   VENDOR_DIRS="${AZTEC_LIB} ${VENDOR_DIRS}"

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
      [  --with-gsl=[lib]      determine the gsl lib (gsl is the default])
 
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

   # define GSL include path
   if test -n "${GSL_INC}" ; then
       # remember that GSL_INC has the final slash
       GSL_H="\"${GSL_INC}az_gsl.h\""
       GSL_DEFS_H="\"${GSL_INC}az_gsl_defs.h\""
   elif test -z "${GSL_INC}" ; then
       GSL_H="<az_gsl.h>"
       GSL_DEFS_H="<az_gsl_defs.h>"
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_gsl=$1

   # set up the libraries
   if test "${with_gsl}" != no ; then
       if test -n "${GSL_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gsl, -L${GSL_LIB} -l${with_gsl})
       elif test -z "${GSL_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gsl, -l${with_gsl})
       fi
   fi

   # add GSL directory to VENDOR_DIRS
   VENDOR_DIRS="${GSL_LIB} ${VENDOR_DIRS}"

])

dnl-------------------------------------------------------------------------dnl
dnl AC_GSLCBLAS_SETUP
dnl
dnl GSLCBLAS SETUP (on by default)
dnl GSLCBLAS is a required vendor
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_GSLCBLAS_SETUP, [dnl

   dnl define --with-gslcblas
   AC_ARG_WITH(gslcblas,
      [  --with-gslcblas=[lib]      determine the gslcblas lib (gslcblas is the default])
 
   dnl define --with-gslcblas-inc
   AC_WITH_DIR(gslcblas-inc, GSLCBLAS_INC, \${GSLCBLAS_INC_DIR},
	       [tell where GSLCBLAS includes are])

   dnl define --with-gslcblas-lib
   AC_WITH_DIR(gslcblas-lib, GSLCBLAS_LIB, \${GSLCBLAS_LIB_DIR},
	       [tell where GSLCBLAS libraries are])

   # set default value of gslcblas includes and libs
   if test "${with_gslcblas:=gslcblas}" = yes ; then
       with_gslcblas='gslcblas'
   fi

   # define GSLCBLAS include path
   if test -n "${GSLCBLAS_INC}" ; then
       # remember that GSLCBLAS_INC has the final slash
       GSLCBLAS_H="\"${GSLCBLAS_INC}az_gslcblas.h\""
       GSLCBLAS_DEFS_H="\"${GSLCBLAS_INC}az_gslcblas_defs.h\""
   elif test -z "${GSLCBLAS_INC}" ; then
       GSLCBLAS_H="<az_gslcblas.h>"
       GSLCBLAS_DEFS_H="<az_gslcblas_defs.h>"
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_gslcblas=$1

   # set up the libraries
   if test "${with_gslcblas}" != no ; then
       if test -n "${GSLCBLAS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gslcblas, -L${GSLCBLAS_LIB} -l${with_gslcblas})
       elif test -z "${GSLCBLAS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gslcblas, -l${with_gslcblas})
       fi
   fi

   # add GSLCBLAS directory to VENDOR_DIRS
   VENDOR_DIRS="${GSLCBLAS_LIB} ${VENDOR_DIRS}"

])

dnl-------------------------------------------------------------------------dnl
dnl AC_TRILINOS_SETUP
dnl
dnl TRILINOS SETUP (on by default)
dnl TRILINOS is a required vendor
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_TRILINOS_SETUP, [dnl

   dnl define --with-trilinos
   AC_ARG_WITH(trilinos,
      [  --with-trilinos=[lib]    determine the trilinos implementation (aztecoo is default])
 
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

   # define TRILINOS include path
   if test -n "${TRILINOS_INC}" ; then
       # remember that TRILINOS_INC has the final slash
       Trilinos_Util_H="\"${TRILINOS_INC}Trilinos_Util.h\""
       AztecOO_H="\"${TRILINOS_INC}AztecOO.h\""
       Epetra_MpiComm_H="\"${TRILINOS_INC}Epetra_MpiComm.h\"" 
       Epetra_Map_H="\"${TRILINOS_INC}Epetra_Map.h\""
       Epetra_Vector_H="\"${TRILINOS_INC}Epetra_Vector.h\""
       Epetra_CrsMatrix_H="\"${TRILINOS_INC}Epetra_CrsMatrix.h\""
       Epetra_LinearProblem_H="\"${TRILINOS_INC}Epetra_LinearProblem.h\""
       Epetra_IntVector_H="\"${TRILINOS_INC}Epetra_IntVector.h\""
       Epetra_Import_H="\"${TRILINOS_INC}Epetra_Import.h\""
       Epetra_Export_H="\"${TRILINOS_INC}Epetra_Export.h\""
       Epetra_CompObject_H="\"${TRILINOS_INC}Epetra_CompObject.h\""
       Epetra_Distributor_H="\"${TRILINOS_INC}Epetra_Distributor.h\""
       Epetra_DistObject_H="\"${TRILINOS_INC}Epetra_DistObject.h\""
       Epetra_MpiDistributor_H="\"${TRILINOS_INC}Epetra_MpiDistributor.h\""
       Epetra_BasicDirectory_H="\"${TRILINOS_INC}Epetra_BasicDirectory.h\""
       Epetra_Util_H="\"${TRILINOS_INC}Epetra_Util.h\""
       Epetra_Time_H="\"${TRILINOS_INC}Epetra_Time.h\""
dnl   elif test -z "${TRILINOS_INC}" ; then
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_trilinos=$1

   # set up the libraries
   if test "${with_trilinos}" != no ; then
       if test -n "${TRILINOS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_trilinos, -L${TRILINOS_LIB} -l${with_trilinos} -lepetra -ltriutils -ly12m)
       elif test -z "${TRILINOS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_trilinos, -l${with_trilinos} -lepetra -ltriutils -ly12m)
       fi
   fi

   # add TRILINOS directory to VENDOR_DIRS
   VENDOR_DIRS="${TRILINOS_LIB} ${VENDOR_DIRS}"

])

dnl-------------------------------------------------------------------------dnl
dnl AC_PCG_SETUP
dnl
dnl PCG LIBRARY SETUP (on by default)
dnl PCG is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_PCG_SETUP, [dnl

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

   # set up the libraries
   if test "${with_pcg}" != no ; then
       if test -z "${PCG_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_pcg, -l${with_pcg})
       elif test -n "${PCG_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_pcg, -L${PCG_LIB} -l${with_pcg})
       fi
   fi

   # add PCG directory to VENDOR_DIRS
   VENDOR_DIRS="${PCG_LIB} ${VENDOR_DIRS}"

   dnl note that we add some system-specific libraries for this
   dnl vendor in AC_DRACO_ENV; also, the user must set up LAPACK for
   dnl this to work
])

dnl-------------------------------------------------------------------------dnl
dnl AC_GANDOLF_SETUP
dnl
dnl GANDOLF LIBRARY SETUP (on by default)
dnl GANDOLF is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_GANDOLF_SETUP, [dnl

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

   # set up the libraries
   if test "${with_gandolf}" != no ; then
       if test -z "${GANDOLF_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gandolf, -l${with_gandolf})
       elif test -n "${GANDOLF_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gandolf, -L${GANDOLF_LIB} -l${with_gandolf})
       fi
   fi

   # add GANDOLF directory to VENDOR_DIRS
   VENDOR_DIRS="${GANDOLF_LIB} ${VENDOR_DIRS}"

   dnl SGI needs "-lfortran" on the load line when including libgandolf.a.
   dnl This library is added to ${LIBS} in AC_DRACO_ENV.
])

dnl-------------------------------------------------------------------------dnl
dnl AC_EOSPAC5_SETUP
dnl
dnl EOSPAC5 LIBRARY SETUP (on by default)
dnl EOSPAC5 is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_EOSPAC5_SETUP, [dnl

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

   # set up the libraries
   if test "${with_eospac}" != no ; then
       if test -z "${EOSPAC5_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_eospac, -l${with_eospac})
       elif test -n "${EOSPAC5_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_eospac, -L${EOSPAC5_LIB} -l${with_eospac})
       fi
   fi

   # add EOSPAC5 directory to VENDOR_DIRS
   VENDOR_DIRS="${EOSPAC5_LIB} ${VENDOR_DIRS}"

   dnl note that we add some system-specific libraries for this
   dnl vendor in AC_DRACO_ENV; also, the user must set up LAPACK for
   dnl this to work
])

dnl-------------------------------------------------------------------------dnl
dnl AC_LAPACK_SETUP
dnl
dnl LAPACK SETUP (on by default)
dnl LAPACK is a required vendor
dnl
dnl NOTE: this also sets up the BLAS
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_LAPACK_SETUP, [dnl

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

   # set up the libraries: if the option is vendor than we will take
   # care of things in dracoenv under each platform case; if the
   # option is atlas we will setup the basic library calls
   if test "${with_lapack}" = atlas ; then
       
       # if a library path has been defined use it otherwise assume
       # the libraries are in a default location
       if test -z "${LAPACK_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_lapack, -llapack -lf77blas -lcblas -latlas)
       elif test -n "${LAPACK_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_lapack, -L${LAPACK_LIB} -llapack -lf77blas -lcblas -latlas)
       fi

   fi

   # add LAPACK directory to VENDOR_DIRS
   VENDOR_DIRS="${LAPACK_LIB} ${VENDOR_DIRS}"

   dnl note that we add system specific libraries to this list in
   dnl dracoenv
])

dnl-------------------------------------------------------------------------dnl
dnl AC_GRACE_SETUP
dnl
dnl GRACE SETUP (on by default)
dnl GRACE is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_GRACE_SETUP, [dnl

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

   # define GRACE include path
   if test -n "${GRACE_INC}" ; then
       # remember that GRACE_INC has the final slash
       GRACE_H="\"${GRACE_INC}${with_grace}.h\""
   elif test -z "${GRACE_INC}" ; then
       GRACE_H="<${with_grace}.h>"
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_grace=$1

   # set up the libraries
   if test "${with_grace}" != no ; then
       if test -n "${GRACE_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_grace, -L${GRACE_LIB} -l${with_grace})
       elif test -z "${GRACE_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_grace, -l${with_grace})
       fi
   fi

   # add GRACE directory to VENDOR_DIRS
   VENDOR_DIRS="${GRACE_LIB} ${VENDOR_DIRS}"

])

dnl-------------------------------------------------------------------------dnl
dnl AC_VENDOR_DEFINES
dnl
dnl Run at the end of the environment setup to add defines required by
dnl the vendors.  We do this to allow platform specific mods to the 
dnl vendor defines BEFORE they are output to config.h, etc files.  
dnl
dnl This macro needs to be updated when new vendors are added.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_VENDOR_DEFINES], [dnl

   AC_MSG_CHECKING("vendor defines")
   defines=''

   # ***
   # MPI
   # ***
   if test -n "${vendor_mpi}"; then
       # define mpi include path
       AC_DEFINE_UNQUOTED(MPI_H, ${MPI_H})dnl

       # add to defines
       defines="${defines} ${MPI_H}"
   fi

   # *****
   # SPRNG
   # *****
   if test -n "${vendor_sprng}"; then
       # define sprng include path
       AC_DEFINE_UNQUOTED(SPRNG_H, ${SPRNG_H})dnl 

       # add to defines
       defines="${defines} ${SPRNG_H}"
   fi

   # *****
   # AZTEC
   # *****
   if test -n "${vendor_aztec}"; then 
       # define aztec include paths
       AC_DEFINE_UNQUOTED(AZTEC_H, ${AZTEC_H})dnl
       AC_DEFINE_UNQUOTED(AZTEC_DEFS_H, ${AZTEC_DEFS_H})dnl

       # add to defines
       defines="${defines} ${AZTEC_H} ${AZTEC_DEFS_H}"
   fi

   # *****
   # GSL
   # *****
   if test -n "${vendor_gsl}"; then 
       # define gsl include paths
       AC_DEFINE_UNQUOTED(GSL_H, ${GSL_H})dnl
       AC_DEFINE_UNQUOTED(GSL_DEFS_H, ${GSL_DEFS_H})dnl

       # add to defines
       defines="${defines} ${GSL_H} ${GSL_DEFS_H}"
   fi

   # *****
   # GSLCBLAS
   # *****
   if test -n "${vendor_gslcblas}"; then 
       # define gslcblas include paths
       AC_DEFINE_UNQUOTED(GSLCBLAS_H, ${GSLCBLAS_H})dnl
       AC_DEFINE_UNQUOTED(GSLCBLAS_DEFS_H, ${GSLCBLAS_DEFS_H})dnl

       # add to defines
       defines="${defines} ${GSLCBLAS_H} ${GSLCBLAS_DEFS_H}"
   fi


   # *****
   # TRILINOS
   # *****
   if test -n "${vendor_trilinos}"; then 
       # define trilinos include paths
       AC_DEFINE_UNQUOTED(Trilinos_Util_H, ${Trilinos_Util_H})dnl
       AC_DEFINE_UNQUOTED(AztecOO_H, ${AztecOO_H})dnl
       AC_DEFINE_UNQUOTED(Epetra_MpiComm_H, ${Epetra_MpiComm_H})dnl
       AC_DEFINE_UNQUOTED(Epetra_Map_H, ${Epetra_Map_H})dnl
       AC_DEFINE_UNQUOTED(Epetra_Vector_H, ${Epetra_Vector_H})dnl
       AC_DEFINE_UNQUOTED(Epetra_CrsMatrix_H, ${Epetra_CrsMatrix_H})dnl
       AC_DEFINE_UNQUOTED(Epetra_LinearProblem_H, ${Epetra_LinearProblem_H})dnl
       AC_DEFINE_UNQUOTED(Epetra_IntVector_H, ${Epetra_IntVector_H})dnl
       AC_DEFINE_UNQUOTED(Epetra_Import_H, ${Epetra_Import_H})dnl
       AC_DEFINE_UNQUOTED(Epetra_Export_H, ${Epetra_Export_H})dnl
       AC_DEFINE_UNQUOTED(Epetra_CompObject_H, ${Epetra_CompObject_H})dnl
       AC_DEFINE_UNQUOTED(Epetra_Distributor_H, ${Epetra_Distributor_H})dnl
       AC_DEFINE_UNQUOTED(Epetra_DistObject_H, ${Epetra_DistObject_H})dnl
       AC_DEFINE_UNQUOTED(Epetra_MpiDistributor_H, ${Epetra_MpiDistributor_H})dnl
       AC_DEFINE_UNQUOTED(Epetra_BasicDirectory_H, ${Epetra_BasicDirectory_H})dnl
       AC_DEFINE_UNQUOTED(Epetra_Util_H, ${Epetra_Util_H})dnl
       AC_DEFINE_UNQUOTED(Epetra_Time_H, ${Epetra_Time_H})dnl

       # add to defines
       defines="${defines} ${Trilinos_Util_H} ${AztecOO_H} ${Epetra_MpiComm_H} ${Epetra_Map_H}"
       defines="${defines} ${Epetra_Vector_H} ${Epetra_CrsMatrix_H} ${Epetra_LinearProblem_H}"
       defines="${defines} ${Epetra_IntVector_H} ${Epetra_Import_H} ${Epetra_Export_H}"
       defines="${defines} ${Epetra_CompObject_H} ${Epetra_Distributor_H} ${Epetra_DistObject_H}"
       defines="${defines} ${Epetra_MpiDistributor_H} ${Epetra_BasicDirectory_H}"
       defines="${defines} ${Epetra_Util_H} ${Epetra_Time_H}"
   fi

   # *****
   # GRACE
   # *****
   if test -n "${vendor_grace}"; then
       # define grace include path
       AC_DEFINE_UNQUOTED(GRACE_H, ${GRACE_H})dnl

       # add to defines
       defines="${defines} ${GRACE_H}"
   fi

   if test -n "${defines}"; then
       AC_MSG_RESULT("${defines}")
   else
       AC_MSG_RESULT("no vendors definitions required")
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
   AC_GSLCBLAS_SETUP(pkg)
   AC_TRILINOS_SETUP(pkg)
   AC_LAPACK_SETUP(pkg)
   AC_GANDOLF_SETUP(pkg)
   AC_EOSPAC5_SETUP(pkg)
   AC_GRACE_SETUP(pkg)
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_vendors.m4
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
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_local.m4
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
   dnl DRACO STANDARD HEADERS
   dnl

   dnl defines --enable-draco-stdhdrs
   AC_ARG_ENABLE(draco-stdhdrs,
      [  --enable-draco-stdhdrs  use draco standard headers (off by default)])

   dnl Doxygen options

   AC_ARG_ENABLE(latex-doc,
      [  --enable-latex-doc      built latex docs with doxygen (off by default)])

   dnl end of AC_DRACO_ARGS
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_dracoarg.m4
dnl-------------------------------------------------------------------------dnl


