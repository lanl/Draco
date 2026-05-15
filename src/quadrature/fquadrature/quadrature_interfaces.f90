!---------------------------------------------*-F90-*---------------------------------------------!
!
! \file   quadrature/quadrature_interfaces.f90
! \author Allan Wollaber, Jae Chang
! \brief  Provides F90 quadrature interfaces
! \note   Copyright (C) 2020-2026 Triad National Security, LLC., All rights reserved.
!-------------------------------------------------------------------------------------------------!
!
! Ref: http://gcc.gnu.org/onlinedocs/gfortran/Interoperable-Subroutines-and-Functions.html
!      http://fortranwiki.org/fortran/show/Generating+C+Interfaces
!
!-------------------------------------------------------------------------------------------------!
! This module provides all of the quadrature interface derived types and function signatures.
!-------------------------------------------------------------------------------------------------!

module quadrature_interfaces
  use iso_c_binding, only: c_double, c_int, c_ptr, c_null_ptr
  implicit none

  ! ---------------------------------------------------------------------------------------
  ! This *must* exactly match the layout in Quadrature_Interface.hh Quadrature Data Struct
  ! ---------------------------------------------------------------------------------------
  type, bind(C) :: quadrature_data
    integer(c_int) :: dimension       !< (1, 2, 3)
    integer(c_int) :: type            !< GL, LOBATO, LS, etc
    integer(c_int) :: order           !< 2, 4, 6, ...
    integer(c_int) :: azimuthal_order !< only for PRODUCT_CHEBYSHE_LEGENDRE
    integer(c_int) :: geometry        !< CARTESIAN, AXISYMMETRIC,SPHERICAL
    type(c_ptr)    :: mu              !< directional cosines
    type(c_ptr)    :: eta             !< directional cosines
    type(c_ptr)    :: xi              !< directional cosines
    type(c_ptr)    :: wt              !< ordinate weights
  end type quadrature_data

  interface

    ! ----------------------------------------------------------------
    ! This sets all pointers to NULL and zeros out data
    ! ----------------------------------------------------------------
    subroutine init_quadrature(quad) bind(C, name="init_quadrature")
      import quadrature_data
      implicit none
      ! Argument list
      type(quadrature_data), intent(inout) :: quad
    end subroutine init_quadrature

    ! -------------------------------------------------------------------------------------------
    ! This function inspects the desired quadrature specifications in the enumerated (integer)
    ! struct values and fills in the ordinate data. It assumes there is space to write the data.
    ! -------------------------------------------------------------------------------------------
    subroutine get_quadrature(quad) bind(C, name="get_quadrature")
      import quadrature_data
      implicit none
      ! Argument list
      type(quadrature_data), intent(inout) :: quad
    end subroutine get_quadrature

  end interface

end module quadrature_interfaces ! LCOV_EXCL_LINE

!-------------------------------------------------------------------------------------------------!
! End quadrature_interfaces.f90
!-------------------------------------------------------------------------------------------------!
