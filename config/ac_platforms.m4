dnl-------------------------------------------------------------------------dnl
dnl ac_platforms.m4
dnl
dnl Defines platform-specfic environments, including default vendor
dnl settings for the CCS-4/ASC computer platforms.
dnl
dnl Thomas M. Evans
dnl 2003/04/30 20:29:39
dnl-------------------------------------------------------------------------dnl
##---------------------------------------------------------------------------##
## $Id$
##---------------------------------------------------------------------------##

dnl-------------------------------------------------------------------------dnl
dnl AC_DBS_PLATFORM_ENVIRONMENT
dnl
dnl Configure draco build system platfrom-specfic variables
dnl This function is called within AC_DRACO_ENV
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_PLATFORM_ENVIRONMENT], [dnl

   # we must know the host
   AC_REQUIRE([AC_CANONICAL_HOST])

   # systems setup
   case $host in

   # ***********
   # LINUX SETUP
   # ***********
   *-linux-gnu)

       # print out cpu message
       AC_MSG_CHECKING("host platform cpu")
       AC_MSG_RESULT("${host_cpu}")

       # we do not do any posix source defines unless the user
       # specifically requests them
       if test "${with_posix:=no}" = yes ; then
	   with_posix='199309L'
       fi

       if test "${with_posix}" != no ; then
	   AC_DEFINE(_POSIX_SOURCE)
	   AC_DEFINE_UNQUOTED(_POSIX_C_SOURCE, $with_posix)
       fi

       #
       # setup linux strict if the compiler is KCC (also turn off the
       # warnings about long long being non-standard)
       #
       if test "${CXX}" = KCC && test -n "${STRICTFLAG}" ; then
	   AC_MSG_WARN("Linux KCC strict option set to allow long long type!")
	   STRICTFLAG="--linux_strict -D__KAI_STRICT --diag_suppress 450"
       fi

       #
       # add thread safety if we are using KCC on linux
       #
       if test "${CXX}" = KCC ; then
	   CFLAGS="--thread_safe ${CFLAGS}"
	   CXXFLAGS="--thread_safe ${CXXFLAGS}"
	   ARFLAGS="--thread_safe ${ARFLAGS}"
	   LDFLAGS="--thread_safe ${LDFLAGS}"
       fi

       # determine word sizes
       # AC_DETERMINE_WORD_SIZES

       #
       # setup communication packages
       #
       
       # the default locations for mpi include/lib are:
       #   /usr/local/mpich/include
       #   /usr/local/mpich/lib
       # to make life easy for CCS-2/4 users; needless to say,
       # these can be overridden by --with-mpi-lib and --with-mpi-inc

       # setup for mpi support, on linux vendor and mich are one
       # and the same because there is no vendor for mpi on linux
        
       if test "${with_mpi}" = vendor ; then
	   with_mpi='mpich'
       fi

       if test "${with_mpi}" = mpich ; then
	   
	   # define mpich libraries for v1.2 of mpich
	   if test -n "${MPI_LIB}" ; then
	       linux_mpi_libs="-L${MPI_LIB} -lmpich"
	   elif test -z "${MPI_LIB}" ; then
	       # if /usr/local/mpich/lib exists use it by default;
	       # this is set as the default for the CCS-2/4 network;
	       # it may not be appropriate on other LINUX networks;
	       # in those cases, override with --with-mpi-lib
	       if test -d /usr/local/mpich/lib ; then
	           linux_mpi_libs="-L/usr/local/mpich/lib -lmpich"
	       else
		   linux_mpi_libs="-lmpich"
	       fi
	   fi

	   # define the linux mpi libs
	   AC_VENDORLIB_SETUP(vendor_mpi, ${linux_mpi_libs})

	   # set the default include location on LINUX to
	   # /usr/local/mpich/include; this is specific to the CCS-2/4
	   # LINUX network; to override on other systems use
	   # --with-mpi-inc on the configure line

	   # if MPI_INC is undefined, then put in the explicit default
	   # path if /usr/local/mpich/include exists
	   if test -z "${MPI_INC}" && test -d "/usr/local/mpich/include" ; then
	       MPI_H="\"/usr/local/mpich/include/mpi.h\""
	       AC_DEFINE_UNQUOTED(MPI_H, ${MPI_H})
	   fi

       fi

       # shmem (not available on suns)
       if test "${enable_shmem}" = yes ; then
	   AC_MSG_ERROR("We do not support shmem on linux!")
       fi

       #
       # end of communication package setup
       #  

       # 
       # setup lapack 
       #
       
       # we assume that the vendor option on linux is the install of
       # redhat rpms in /usr/lib; we don't worry about atlas because
       # that has already been defined

       if test "${with_lapack}" = vendor ; then

	   # if an lapack location was defined use it
	   if test -n "${LAPACK_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_lapack, -L${LAPACK_LIB} -llapack -lblas)
	   elif test -z "${LAPACK_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_lapack, -llapack -lblas)
	   fi

       fi 

       # 
       # end of lapack setup
       # 

       #
       # setup eospac
       #
       
       AC_MSG_CHECKING("for extra eospac library requirements.")
       if test -n "${vendor_eospac}"; then
	   extra_eospac_libs="-L/usr/local/lf9560/lib -lfj9i6 -lfj9e6 -lfj9f6 -lfst -lfccx86"
           LIBS="${LIBS} ${extra_eospac_libs}"
           AC_MSG_RESULT("${extra_eospac_libs}")
       else
           AC_MSG_RESULT("none.")
       fi

       #
       # end of eospac
       #

       #
       # add libg2c to LIBS if lapack, gandolf, or pcg is used
       #
       AC_MSG_CHECKING("libg2c requirements")
       if test -n "${vendor_lapack}" || test -n "${vendor_pcg}" ||
	  test -n "${vendor_gandolf}"; then
	   
	   # Add g2c for various compilers
	   if test "${CXX}" = KCC ; then
	       LIBS="${LIBS} --backend -lg2c"
	       AC_MSG_RESULT("--backend -lg2c added to LIBS")
	   elif test "${CXX}" = g++ ; then
	       LIBS="${LIBS} -lg2c"
	       AC_MSG_RESULT("-lg2c added to LIBS")
	   elif test "${CXX}" = icc ; then
               AC_PATH_PROG(GCC_BIN, g++, null)
               GCC_BIN=`dirname ${GCC_BIN}`
               GCC_HOME=`dirname ${GCC_BIN}`
               GCC_LIB_DIR="${GCC_HOME}/lib"
	       LIBS="${LIBS} -L${GCC_LIB_DIR} -lg2c"
	       AC_MSG_RESULT("-lg2c added to LIBS")
	   fi

       else
	   AC_MSG_RESULT("not needed")
       fi

       # set rpath when building shared library executables
       if test "${enable_shared}" = yes ; then

	   # turn off ranlib
	   RANLIB=':'

	   # the g++/icc rpath needs Xlinker in front of it
	   if test "${CXX}" = g++ || test "${CXX}" = icc; then
	       RPATHA="-Xlinker -rpath \${curdir}"
	       RPATHB="-Xlinker -rpath \${curdir}/.."
	       RPATHC="-Xlinker -rpath \${libdir}"
	       RPATH="${RPATHA} ${RPATHB} ${RPATHC} ${RPATH}"
	   else
	       RPATH="-rpath \${curdir}:\${curdir}/..:\${libdir} ${RPATH}"
	   fi

       fi

       # add vendors to rpath
       for vendor_dir in ${VENDOR_DIRS}; 
       do
	   # if we are using gcc/icc then add xlinker
	   if test "${CXX}" = g++ || test "${CXX}" = icc; then
	       RPATH="-Xlinker -rpath ${vendor_dir} ${RPATH}"

	   # else we just add the rpath
	   else
	       RPATH="-rpath ${vendor_dir} ${RPATH}"
	   fi
       done

       # add the intel math library for better performance when
       # compiling with intel
       if test "${CXX}" = icc; then
	   LIBS="$LIBS -limf"
       fi
   ;;

   # *********
   # SGI SETUP
   # *********
   mips-sgi-irix6.*)
   
       # print out cpu message
       AC_MSG_CHECKING("host platform cpu")
       AC_MSG_RESULT("${host_cpu}")

       # posix source defines, by default we set posix on 
       if test "${with_posix:=yes}" = yes ; then
	   with_posix='199309L'
       fi

       if test "${with_posix}" != no ; then
	   AC_DEFINE_UNQUOTED(_POSIX_C_SOURCE, $with_posix)
	   AC_DEFINE(_POSIX_SOURCE)
       fi

       # RANLIB TAG ON SGI
       RANLIB=':'

       # BIT COMPILER FLAGS ON SGI
       if test "${enable_32_bit:=no}" = yes ; then
	   if test "${with_cxx}" = gcc ; then
	       CXXFLAGS="-mabi=n32 ${CXXFLAGS}"
	       CFLAGS="-mabi=n32 ${CFLAGS}"
	       if test "${enable_shared}" = yes ; then
		   ARFLAGS="-mabi=n32 ${ARFLAGS}"
	       fi
	       LDFLAGS="-mabi=n32 ${LDFLAGS}"
	   else
	       CXXFLAGS="-n32 ${CXXFLAGS}"
	       CFLAGS="-n32 ${CFLAGS}"
	       ARFLAGS="-n32 ${ARFLAGS}"
	       LDFLAGS="-n32 ${LDFLAGS}"
	   fi
       else 
	   if test "${with_cxx}" = gcc ; then
	       CXXFLAGS="-mabi=64 ${CXXFLAGS}"
	       CFLAGS="-mabi=64 ${CFLAGS}"
	       if test "${enable_shared}" = yes ; then
		   ARFLAGS="-mabi=64 ${ARFLAGS}"
	       fi
	       LDFLAGS="-mabi=64 ${LDFLAGS}"
	   else
	       CXXFLAGS="-64 ${CXXFLAGS}"
	       CFLAGS="-64 ${CFLAGS}"
	       ARFLAGS="-64 ${ARFLAGS}"
	       LDFLAGS="-64 ${LDFLAGS}"
	   fi
       fi

       # MIPS INSTRUCTIONS ON SGI
       # this is different depending upon the compiler
       if test "${with_cxx}" = kcc ; then
	   CXXFLAGS="-mips${with_mips:=4} --backend -r10000 ${CXXFLAGS}"
	   CFLAGS="-mips${with_mips:=4} -r10000 ${CFLAGS}"
	   ARFLAGS="-mips${with_mips:=4} ${ARFLAGS}"
	   LDFLAGS="-mips${with_mips:=4} ${LDFLAGS}"
       elif test "${with_cxx}" = sgi ; then
	   CXXFLAGS="-mips${with_mips:=4} -r10000 ${CXXFLAGS}"
	   CFLAGS="-mips${with_mips:=4} -r10000 ${CFLAGS}"
	   ARFLAGS="-mips${with_mips:=4} ${ARFLAGS}"
	   LDFLAGS="-mips${with_mips:=4} ${LDFLAGS}"
       elif test "${with_cxx}" = gcc ; then
	   CXXFLAGS="-mips${with_mips:=4} ${CXXFLAGS}"
	   CFLAGS="-mips${with_mips:=4} ${CFLAGS}"
	   if test "${enable_shared}" = yes ; then
	       ARFLAGS="-mips${with_mips:=4} ${ARFLAGS}"
	   fi
	   LDFLAGS="-mips${with_mips:=4} ${LDFLAGS}"
       fi

       #
       # setup communication packages
       #
   
       # setup for mpi support
       # we only support vendor mpi on sgis       
       if test "${with_mpi}" = vendor ; then
	   if test -n "${MPI_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_mpi, -L${MPI_LIB} -lmpi)
	   elif test -z "${MPI_LIB}" ; then
	       if test "${enable_32_bit}" = no ; then
		   AC_VENDORLIB_SETUP(vendor_mpi, -L${MPT_SGI}/usr/lib64 -lmpi)
	       else
		   AC_VENDORLIB_SETUP(vendor_mpi, -L${MPT_SGI}/usr/lib32 -lmpi)
	       fi
	   fi
       elif test "${with_mpi}" = mpich ; then
	   AC_MSG_ERROR("We do not support mpich on the SGI yet!")
       fi

       # MPT (Message Passing Toolkit) for SGI vendor
       # implementation of MPI and SHMEM
       if test -z "${MPI_INC}" &&  test "${with_mpi}" = vendor ; then
	   MPI_INC="${MPT_SGI}/usr/include/"	   
	   MPI_H="\"${MPI_INC}mpi.h\""
	   AC_DEFINE_UNQUOTED(MPI_H, ${MPI_H})
       fi

       # add SGI MPT Specfic options
       if test "${with_mpi}" = vendor ; then
	   # define no C++ bindings
	   AC_DEFINE(MPI_NO_CPPBIND)
       fi

       #
       # end of communication package setup
       #

       #
       # adjust strict flag for KCC
       #

       if test "${CXX}" = KCC && test "${enable_strict_ansi}" = yes ; then
	   
	   # if integer type is long long or vendor mpi is on then we
	   # need to allow long long type
	   if test "${with_mpi}" = vendor; then 
		   AC_MSG_WARN("KCC strict option set to allow long long type")
		   STRICTFLAG="${STRICTFLAG} --diag_suppress 450"
	   fi

       fi

       #
       # setup lapack
       #

       if test "${with_lapack}" = vendor ; then

	   # if an lapack location was defined use it
	   if test -n "${LAPACK_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_lapack, -L${LAPACK_LIB} -lcomplib.sgimath)
	   elif test -z "${LAPACK_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_lapack, -lcomplib.sgimath)
	   fi

       fi

       #
       # end of lapack setup
       #

       #
       # gandolf and eospac requires -lfortran on the link line.
       #

       AC_MSG_CHECKING("libfortran requirements")
       if test -n "${vendor_gandolf}" || test -n "${vendor_eospac}" ; then
          LIBS="${LIBS} -lfortran"
          AC_MSG_RESULT("-lfortran added to LIBS")
       else
	   AC_MSG_RESULT("not needed")
       fi
       
       #
       # end of gandolf/libfortran setup
       #

       # set rpath when building shared library executables
       if test "${enable_shared}" = yes; then

	   # the g++ rpath needs Xlinker in front of it
	   if test "${CXX}" = g++; then
	       RPATHA="-Xlinker -rpath \${curdir}"
	       RPATHB="-Xlinker -rpath \${curdir}/.."
	       RPATHC="-Xlinker -rpath \${libdir}"
	       RPATH="${RPATHA} ${RPATHB} ${RPATHC} ${RPATH}"
	   else
	       RPATH="-rpath \${curdir}:\${curdir}/..:\${libdir} ${RPATH}"
	   fi

       fi

       # add vendors to rpath
       for vendor_dir in ${VENDOR_DIRS}; 
       do
	   # if we are using gcc then add xlinker
	   if test "${CXX}" = g++; then
	       RPATH="-Xlinker -rpath ${vendor_dir} ${RPATH}"

	   # else we just add the rpath
	   else
	       RPATH="-rpath ${vendor_dir} ${RPATH}"
	   fi
       done
   ;;

   # ******************
   # TRU64 COMPAQ SETUP
   # ******************
   alpha*-dec-osf*)
   
       # print out cpu message
       AC_MSG_CHECKING("host platform cpu")
       AC_MSG_RESULT("${host_cpu}")

       # posix source defines, by default we set posix off
       if test "${with_posix:=no}" = yes ; then
	   with_posix='199309L'
       fi

       if test "${with_posix}" != no ; then
	   AC_DEFINE_UNQUOTED(_POSIX_C_SOURCE, $with_posix)
	   AC_DEFINE(_POSIX_SOURCE)
       fi

       #
       # setup communication packages
       #
       
       # setup vendor mpi
       if test "${with_mpi}" = vendor ; then

	   # set up libraries (the headers are already set)
	   if test -n "${MPI_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_mpi, -L${MPI_LIB} -lmpi)
	   elif test -z "${MPI_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_mpi, -lmpi)
	   fi
       
       # setup mpich
       elif test "${with_mpi}" = mpich ; then

	   # set up libraries (the headers are already set)
	   if test -n "${MPI_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_mpi, -L${MPI_LIB} -lmpich)
	   elif test -z "${MPI_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_mpi, -lmpich)
	   fi
   
       fi

       # add COMPAQ ALASKA Specfic options
       if test "${with_mpi}" = vendor ; then
	   # define version check
	   AC_DEFINE(MPI_VERSION_CHECK)
       fi

       #
       # end of communication packages
       #

       #
       # setup lapack
       #

       if test "${with_lapack}" = vendor ; then

	   # if an lapack location was defined use it
	   if test -n "${LAPACK_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_lapack, -L${LAPACK_LIB} -ldxml)
	   elif test -z "${LAPACK_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_lapack, -ldxml)
	   fi

       fi

       #
       # end of lapack setup
       #

       #
       # gandolf, eospac, pcg require -lfor on the link line.
       #

       AC_MSG_CHECKING("libfortran requirements")
       if test -n "${vendor_gandolf}" || test -n "${vendor_eospac}" ||
          test -n "${vendor_pcg}"; then
          LIBS="${LIBS} -lfor"
          AC_MSG_RESULT("-lfor added to LIBS")
       else
	   AC_MSG_RESULT("not needed")
       fi

       #
       # end of gandolf/libfortran setup
       #

       # set rpath when building shared library executables
       if test "${enable_shared}" = yes; then

	   # turn off ranlib
	   RANLIB=':'

	   # the g++ rpath needs Xlinker in front of it
	   if test "${CXX}" = g++; then
	       RPATHA="-Xlinker -rpath \${curdir}"
	       RPATHB="-Xlinker -rpath \${curdir}/.."
	       RPATHC="-Xlinker -rpath \${libdir}"
	       RPATH="${RPATHA} ${RPATHB} ${RPATHC} ${RPATH}"
	   else
	       RPATH="-rpath \${curdir}:\${curdir}/..:\${libdir} ${RPATH}"
	   fi

       fi

       # add vendors to rpath
       for vendor_dir in ${VENDOR_DIRS}; 
       do
	   # if we are using gcc then add xlinker
	   if test "${CXX}" = g++; then
	       RPATH="-Xlinker -rpath ${vendor_dir} ${RPATH}"

	   # else we just add the rpath
	   else
	       RPATH="-rpath ${vendor_dir} ${RPATH}"
	   fi
       done
   ;;

   # *************
   # IBM AIX SETUP
   # *************
   *ibm-aix*)
   
       # print out cpu message
       AC_MSG_CHECKING("host platform cpu")
       AC_MSG_RESULT("${host_cpu}")

       # posix source defines, by default we set posix off
       if test "${with_posix:=no}" = yes ; then
	   with_posix='199309L'
       fi

       if test "${with_posix}" != no ; then
	   AC_DEFINE_UNQUOTED(_POSIX_C_SOURCE, $with_posix)
	   AC_DEFINE(_POSIX_SOURCE)
       fi

       # set up 32 or 64 bit compiling on IBM
       if test "${enable_32_bit:=no}" = yes ; then
	   CXXFLAGS="${CXXFLAGS} -q32"
       elif test "${enable_64_bit:=no}" = yes ; then
	   CXXFLAGS="${CXXFLAGS} -q64"
       fi

       # set up the heap size
       if test "${with_cxx}" = asciwhite ; then
	   LDFLAGS="${LDFLAGS} -bmaxdata=0x80000000"
       fi

       #
       # setup communication packages
       #
       
       # setup vendor mpi
       if test "${with_mpi}" = vendor ; then

	   # set up libraries (the headers are already set)
	   if test -n "${MPI_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_mpi, -L${MPI_LIB} -lmpi)
	   elif test -z "${MPI_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_mpi, -lmpi)
	   fi
       
       # setup mpich
       elif test "${with_mpi}" = mpich ; then

	   # set up libraries (the headers are already set)
	   if test -n "${MPI_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_mpi, -L${MPI_LIB} -lmpich)
	   elif test -z "${MPI_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_mpi, -lmpich)
	   fi
   
       fi

       #
       # end of communication packages
       #

       #
       # setup lapack
       #

       if test "${with_lapack}" = vendor ; then

	   # if an lapack location was defined use it
	   if test -n "${LAPACK_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_lapack, -L${LAPACK_LIB} -ldxml)
	   elif test -z "${LAPACK_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_lapack, -ldxml)
	   fi

       fi

       #
       # end of lapack setup
       #

       #
       # gandolf, eospac, pcg require -lfor on the link line.
       #

       AC_MSG_CHECKING("libfortran requirements")
       if test -n "${vendor_gandolf}" || test -n "${vendor_eospac}" ||
          test -n "${vendor_pcg}"; then
          LIBS="${LIBS} -lfor"
          AC_MSG_RESULT("-lfor added to LIBS") 
       else
	   AC_MSG_RESULT("not needed")
       fi

       #
       # end of gandolf/libfortran setup
       #

       # RPATH is derived from -L, don't need explicit setup

       # do shared specific stuff
       if test "${enable_shared}" = yes ; then
	   # turn off ranlib
	   RANLIB=':'
       fi
   ;;

   # *****************
   # SUN/SOLARIS SETUP
   # *****************
   sparc-sun-solaris2.*)
   
       # print out cpu message
       AC_MSG_CHECKING("host platform cpu")
       AC_MSG_RESULT("${host_cpu}")

       # posix source defines, by default we set poaix on 
       if test "${with_posix:=yes}" = yes ; then
	   with_posix='199309L'
       fi

       if test "${with_posix}" != no ; then
	   AC_DEFINE_UNQUOTED(_POSIX_C_SOURCE, $with_posix)
	   AC_DEFINE(_POSIX_SOURCE)
       fi

       #
       # setup communication packages
       #
   
       # setup for mpi support
       # we only support mpich on sgis       
       if test "${with_mpi}" = vendor ; then
	   AC_MSG_ERROR("We do not support vendor mpi on the SUN yet!")
       elif test "${with_mpi}" = mpich ; then
	   
	   # define sun-required libraries for mpich, v 1.0 (this
	   # needs to be updated for version 1.2)
	   sun_mpi_libs='-lpmpi -lmpi -lsocket -lnsl'
   
	   if test -n "${MPI_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_mpi, -L${MPI_LIB} ${sun_mpi_libs})
	   elif test -z "${MPI_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_mpi, ${sun_mpi_libs})
	   fi

       fi

       # shmem (not available on suns)
       if test "${enable_shmem}" = yes ; then
	   AC_MSG_ERROR("We do not support shmem on suns!")
       fi

       #
       # end of communication package setup
       #

       #
       # setup lapack
       #

       if test "${with_lapack}" = vendor ; then

	   sun_libs='-llapack -lblas -lF77 -lM77 -lsunmath'

	   # if an lapack location was defined use it
	   if test -n "${LAPACK_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_lapack, -L${LAPACK_LIB} ${sun_libs})
	   elif test -z "${LAPACK_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_lapack, ${sun_libs})
	   fi

       fi

       #
       # end of lapack setup
       #

       # set -R when building shared library executables
       if test "${enable_shared}" = yes; then
	   RPATH="-R \${curdir}:\${curdir}/..:\${libdir} ${RPATH}"
	   RANLIB=':'
       fi

       # add vendors to rpath
       for vendor_dir in ${VENDOR_DIRS}; 
       do
	   # if we are using gcc/icc then add xlinker
	   if test "${CXX}" = g++ || test "${CXX}" = icc; then
	       RPATH="-Xlinker -R ${vendor_dir} ${RPATH}"

	   # else we just add the rpath
	   else
	       RPATH="-R ${vendor_dir} ${RPATH}"
	   fi
       done
   ;;

   # *******
   # NOTHING
   # *******
   *)
       AC_MSG_ERROR("Cannot figure out the platform or host!")
   ;;
   esac
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_platforms.m4
dnl-------------------------------------------------------------------------dnl