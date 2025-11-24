!---------------------------------------------*-F90-*---------------------------------------------!
!
! file   FortranCheck/f90sub/drel.f90
! author Kelly Thompson
! date   Tuesday, Jun 12, 2012, 16:03 pm
! brief  Test F90 main linking against C++ library and calling a C++ function.
! note   Copyright (C) 2014-2025 Triad National Security, LLC., All rights reserved.
!-------------------------------------------------------------------------------------------------!

! Ref: http://gcc.gnu.org/onlinedocs/gfortran/Interoperable-Subroutines-and-Functions.html
!      http://fortranwiki.org/fortran/show/Generating+C+Interfaces

subroutine drelf90(nf) bind(c, name="drelf90")

  use iso_c_binding, only: c_size_t, C_NULL_CHAR, c_int
  implicit none
  integer(c_int), intent(out) :: nf

  interface ec_release
    ! include "dsxx//Release.hh"
    subroutine ec_release(release_string, maxlen) bind(C, name="ec_release")
      use iso_c_binding, only: c_char, c_size_t
      implicit none
      character(kind=c_char, len=1), intent(out)  :: release_string
      integer(c_size_t), intent(in), value :: maxlen
    end subroutine ec_release
  end interface ec_release

  interface dsxx_is_big_endian
    ! include "dsxx/Endian.hh"
    function dsxx_is_big_endian() bind(C, name="dsxx_is_big_endian")
      use iso_c_binding, only: c_int
      implicit none
      integer(c_int) :: dsxx_is_big_endian
    end function dsxx_is_big_endian
  end interface dsxx_is_big_endian

  interface dsxx_byte_swap
    ! include "dsxx/Endian.hh"
    subroutine dsxx_byte_swap_int(data) bind(C, name="dsxx_byte_swap_int")
      use iso_c_binding, only: c_int
      implicit none
      integer(c_int), intent(inout) :: data
    end subroutine dsxx_byte_swap_int
    subroutine dsxx_byte_swap_int64_t(data) bind(C, name="dsxx_byte_swap_int64_t")
      use iso_c_binding, only: c_int64_t
      implicit none
      integer(c_int64_t), intent(inout) :: data
    end subroutine dsxx_byte_swap_int64_t
    subroutine dsxx_byte_swap_double(data) bind(C, name="dsxx_byte_swap_double")
      use iso_c_binding, only: c_double
      implicit none
      real(c_double), intent(inout) :: data
    end subroutine dsxx_byte_swap_double
  end interface dsxx_byte_swap

  !----------------------------------------------------------------------
  ! Variable declarations

  integer, parameter :: maxlen = 80
  character(len=maxlen) :: release_string

  integer :: is_big_endian, idata
  real(8) :: ddata

  !----------------------------------------------------------------------
  ! Initialization

  nf = 0 ! init number of failures to zero
  release_string = repeat(' ', maxlen)
  release_string(maxlen:maxlen) = C_NULL_CHAR

  !----------------------------------------------------------------------
  ! Retrieve the version string from dsxx
  call ec_release(release_string, len(release_string, kind=c_size_t))

  print '(a)', trim(release_string)
  if (release_string(1:6) == "Draco-") then
    print '(a)', "Test: passed"
    print '(a)', "     Found 'Draco-' in release string."
  else
    ! LCOV_EXCL_START
    print '(a)', "Test: failed"
    print '(a)', "     Did not find 'Draco-' in the release string."
    nf = nf + 1
    ! LCOV_EXCL_STOP
  end if

  !----------------------------------------------------------------------
  ! Check the dsxx/Endian extern "C" functions...

  is_big_endian = dsxx_is_big_endian()
  if (is_big_endian > 1.or.is_big_endian < 0) then
    ! LCOV_EXCL_START
    print '(a)', "Test: failed"
    print '(a)', "     dsxx_is_big_endian returned an unexpected value."
    nf = nf + 1
    ! LCOV_EXCL_STOP
  end if
  ! note: integers must be signed in F90 (i.e.: we cannot use Z'DEADBEEF')
  idata = 1122867 ! 1122867 = z'00112233'
  call dsxx_byte_swap(idata)
  if (idata /= 857870592) then ! 857870592 = z'33221100'
    ! LCOV_EXCL_START
    print '(a)', "Test: failed"
    print '(a)', "     dsxx_byte_swap(int) returned an unexpected value."
    nf = nf + 1
    ! LCOV_EXCL_STOP
  end if
  ddata = 42.0
  ! Call swap 2x to get initial value
  call dsxx_byte_swap(ddata)
  call dsxx_byte_swap(ddata)
  if (ddata /= 42.0) then
    ! LCOV_EXCL_START
    print '(a)', "Test: failed"
    print '(a)', "     dsxx_byte_swap(double) returned an unexpected value."
    nf = nf + 1
    ! LCOV_EXCL_STOP
  end if

  if (nf > 0) then
    ! LCOV_EXCL_START
    print '(a)', "Test: failed"
    print '(a)', "     Endianess checks had some failures."
    ! LCOV_EXCL_STOP
  else
    print '(a)', "Test: passed"
    print '(a)', "     Endianess checks all pass."
  end if

  !----------------------------------------------------------------------
  ! Summary

  print '(a)', " "
  print '(a)', "*********************************************"
  if (nf > 0) then
    print '(a)', "**** cppmain Test: FAILED." ! LCOV_EXCL_LINE
  else
    print '(a)', "**** cppmain Test: PASSED."
  end if
  print '(a)', "*********************************************"

end subroutine drelf90

!-------------------------------------------------------------------------------------------------!
! end of drel.f90
!-------------------------------------------------------------------------------------------------!
