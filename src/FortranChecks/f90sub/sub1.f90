!----------------------------------------------*-F90-*---------------------------------------------!
! Filename: sub1.f90
! Author:   Kelly Thompson
! Date:     Tue June 12 2012
! Brief:    Fortran subroutine for cross language link checks.
! Note:     Copyright (C) 2014-2025 Triad National Security, LLC., All rights reserved.
!--------------------------------------------------------------------------------------------------!

! http://gcc.gnu.org/onlinedocs/gfortran/Interoperable-Subroutines-and-Functions.html
! http://www.fortran.bcs.org/2002/interop.htm (Example of C calling Fortran)
! http://software.intel.com/sites/products/documentation/hpc/compilerpro/en-us/fortran/lin/compiler_f/bldaps_for/common/bldaps_interopc.htm

! use iso_c_binding

subroutine sub1(alpha, np, nf) bind(c)
  use iso_c_binding, only: c_double, c_size_t
  implicit none
  real(c_double), value, intent(in) :: alpha
  integer(c_size_t), intent(inout) :: np
  integer(c_size_t), intent(inout) :: nf

  ! local variables
  double precision :: small

  !----------------------------------------
  small = 1.0d-13

  write (*, '(a,f5.1,2i3)') "Hello, world.", alpha, np, nf

  if (alpha > 1.0 - small.and.alpha < 1.0 + small) then
    print '(a)', "Test: passed"
    print '(a)', "     alpha == 1.0"
    np = np + 1
  else
    ! LCOV_EXCL_START
    print '(a)', "Test: failed"
    print '(a)', "     alpha != 1.0"
    nf = nf + 1
    ! LCOV_EXCL_STOP
  end if

end subroutine sub1

!-------------------------------------------------------------------------------------------------!
! End sub1.f90
!-------------------------------------------------------------------------------------------------!
