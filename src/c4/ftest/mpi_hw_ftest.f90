!----------------------------------*-F90-*----------------------------------//
!!
!* \file   c4/ftest/mpi_hw_ftest.f90
!* \author Allan Wollaber, Tom Evans, Kelly Thompson
!* \brief  Infinite medium Draco shunt.
!* \note   Copyright (C) 2020-2026 Triad National Security, LLC., All rights reserved.
!---------------------------------------------------------------------------//

subroutine tst_mpi_hw_f(nf) bind(C, name="tst_mpi_hw")

  use draco_mpi
  use fc4_test

  use iso_c_binding, only: c_int, c_double

  implicit none

  integer(c_int), intent(inout) :: nf
  integer :: ierr

  ! ------------------------------------
  ! Initialize MPI
  ! ------------------------------------
  call fc4_mpi_init(ierr)
  call check_fail(ierr, fc4_rank)

  ! ------------------------------------
  ! Run the problem
  ! ------------------------------------

  if (fc4_rank < fc4_num_ranks) then
    call pass_msg(fc4_rank, "MPI rank index ok")
  else
    call it_fails(fc4_rank, "MPI rank > max")
  end if

  call fc4_mpi_barrier(ierr)
  print *, "Hello from rank ", fc4_rank, "/", fc4_num_ranks

  ! ------------------------------------
  ! Finalize and clean up
  ! ------------------------------------

  call fc4_mpi_finalize(ierr)
  call check_fail(ierr, fc4_rank)

  ! Print the overall test result
  ! call test_report(fc4_rank,nf)

  nf = nf + f90_num_failures

end subroutine tst_mpi_hw_f

!--------------------------------------------------------------------------------------------------!
! End of c4/ftest/mpi_hw_ftest.f90
!--------------------------------------------------------------------------------------------------!
