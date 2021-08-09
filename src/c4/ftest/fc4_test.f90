!---------------------------------------------*-F90-*---------------------------------------------!
!
! \file   c4/fc4/fc4_test.f90
! \author Allan Wollaber, Kelly Thompson
! \date   Mon Jul 30 07:06:24 MDT 2012
! \brief  Helper functions for the F90 Draco tests
! \note   Copyright (c) 2016-2021 Triad National Security, LLC., All rights reserved.
!
! This is a modified version of jayenne/src/api/ftest/API_Test.F90.
!-------------------------------------------------------------------------------------------------!

module fc4_test
  use iso_c_binding, only: c_double
  implicit none

  integer, save ::  f90_num_failures = 0

contains

  ! ---------------------------------------------------------------------------
  ! Provide a routine that checks an error code and reports a failure
  ! ---------------------------------------------------------------------------
  subroutine check_fail(ierr, rank)
    implicit none
    integer, intent(in) :: ierr
    integer, intent(in) :: rank

    if (ierr /= 0) then
      write (*, '("**** Test: FAILED on ", I3, " with error ", I3)') rank, ierr
      f90_num_failures = f90_num_failures + 1
    end if
  end subroutine check_fail

  ! ---------------------------------------------------------------------------
  ! Provide a routine that just reports a failure with an error message
  ! ---------------------------------------------------------------------------
  subroutine it_fails(rank, msg)
    implicit none
    integer, intent(in) :: rank
    character(*), intent(in) :: msg

    write (*, '("**** Test: FAILED on ", I3, ":  ", A)') rank, msg

    f90_num_failures = f90_num_failures + 1

  end subroutine it_fails

  ! ---------------------------------------------------------------------------
  ! Provide a routine to report a passing message
  ! ---------------------------------------------------------------------------
  subroutine pass_msg(rank, msg)
    implicit none
    integer, intent(in) :: rank
    character(*), intent(in) :: msg

    write (*, '("**** Test: PASSED on ", I3, ":  ", A)') rank, msg

  end subroutine pass_msg

end module fc4_test

!-------------------------------------------------------------------------------------------------!
! End fc4_test.f90
!-------------------------------------------------------------------------------------------------!
