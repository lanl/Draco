!---------------------------------------------*-F90 *----------------------------------------------!
! \file   c4/fc4/Draco_MPI.F90
! \author Allan Wollaber
! \date   Mon Jul 30 07:06:24 MDT 2012
! \brief  Helper functions to support scalar vs. distributed MPI tests.
! \note   Copyright (c) 2016-2021 Triad National Security, LLC., All rights reserved.
!--------------------------------------------------------------------------------------------------!

module draco_mpi
  use iso_c_binding, only: c_double, c_intptr_t
  implicit none

  integer, public, save :: fc4_rank, fc4_num_ranks

  ! Don't pollute the global namespace with the MPI stuff
  private
#ifdef C4_MPI
  include 'mpif.h'
#endif

  ! Make all the subroutines defined below public, though
  public check_mpi_error

  public fc4_mpi_init
  public fc4_mpi_finalize
  public fc4_mpi_barrier

  interface

     ! Set c4's global 'initialized' boolean.
     subroutine setMpiInit() bind(C, name="setMpiInit")
       implicit none
     end subroutine setMpiInit

  end interface

contains

  ! ---------------------------------------------------------------------------
  ! This subroutine checks for and prints the meaning of MPI errors
  ! ---------------------------------------------------------------------------
  subroutine check_mpi_error(error)
    implicit none
    integer, intent(in)             :: error
#ifdef C4_MPI
    integer                         :: error_string_len, ierror
    character(MPI_MAX_ERROR_STRING) :: error_string
    external mpi_error_string, MPI_Abort

    ! Check and report a nonzero error code
    if (error .ne. 0) then
       call mpi_error_string(error, error_string, error_string_len, ierror)
       write (*, "('*** mpi error = ',i18, ' (', a, ')')") error, trim(error_string)
       call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
    end if
#endif
  end subroutine check_mpi_error

  ! ---------------------------------------------------------------------------
  ! A simple MPI initialization function that also sets rank/num_ranks info.
  ! ---------------------------------------------------------------------------
  subroutine fc4_mpi_init(ierr)
    implicit none
    integer, intent(out) :: ierr
    external mpi_init, mpi_comm_size, mpi_comm_rank

    ierr = 0
    fc4_rank = 0
    fc4_num_ranks = 1
#ifdef C4_MPI
    call mpi_init(ierr)
    call check_mpi_error(ierr)
    call setMpiInit

    call mpi_comm_size(MPI_COMM_WORLD, fc4_num_ranks, ierr)
    call check_mpi_error(ierr)

    call mpi_comm_rank(MPI_COMM_WORLD, fc4_rank, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine fc4_mpi_init

  ! ---------------------------------------------------------------------------
  ! A simple MPI finalize function to wrap the MPI dependencies
  ! ---------------------------------------------------------------------------
  subroutine fc4_mpi_finalize(ierr)
    implicit none
    integer, intent(out) :: ierr
#ifdef C4_MPI
    external mpi_finalize
#endif
    ierr = 0
#ifdef C4_MPI
    call mpi_finalize(ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine fc4_mpi_finalize

  ! ---------------------------------------------------------------------------
  ! Global barrier
  ! ---------------------------------------------------------------------------
  subroutine fc4_mpi_barrier(ierr)
    implicit none
    integer, intent(out) :: ierr
#ifdef C4_MPI
    external mpi_barrier
#endif
    ierr = 0
#ifdef C4_MPI
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine fc4_mpi_barrier

end module draco_mpi

!--------------------------------------------------------------------------------------------------!
! End c4/fc4/Draco_MPI.F90
!--------------------------------------------------------------------------------------------------!
