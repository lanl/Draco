dnl-------------------------------------------------------------------------dnl
dnl ac_vendors.m4
dnl
dnl Macros for each vendor that is used supported by the Draco build
dnl system.
dnl
dnl Thomas M. Evans
dnl 1999/02/04 01:56:22
dnl-------------------------------------------------------------------------dnl
##---------------------------------------------------------------------------##
## $Id$
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## All vendor macros should take the following arguments:
##     pkg      - this vendor is used in the package (default)
##     test     - this vendor is used only in the package's test
##
## Each vendor requires an AC_<VENDOR>_SETUP function.  Additionally,
## any AC_DEFINE or AC_DEFINE_UNQUOTED macros needed by the vendors
## should be done inside of AC_VENDOR_DEFINES.
##---------------------------------------------------------------------------##

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

