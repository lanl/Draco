!---------------------------------------------*-F90-*-----------------------------------------------
!
! file   quadrature/ftest/tstquadrature_interfaces.f90
! author Allan Wollaber
! date   Tuesday, Jun 12, 2012, 16:03 pm
! brief  Test F90 quadrature_data passed into a C++ function
! note   Copyright (C) 2016-2023 Triad National Security, LLC., All rights reserved.
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! This F90 function interface allows the Fortran function test_quadrature_interfaces to directly
! call the extern "C" function rtt_test_quadrature_interfaces (from the library
! rtt_quadrature_test).
!---------------------------------------------------------------------------------------------------
module rtt_test_quadrature_f
  implicit none

  ! Now create an interface to the C routine that accepts that a user-defined type
  interface
    subroutine rtt_test_quadrature_interfaces(quad, err_code) &
      bind(C, name="rtt_test_quadrature_interfaces")
      use iso_c_binding, only: c_int

      use quadrature_interfaces, only: quadrature_data
      implicit none
      type(quadrature_data), intent(in)  :: quad
      integer(c_int), intent(out) :: err_code
    end subroutine rtt_test_quadrature_interfaces
  end interface

end module rtt_test_quadrature_f

!---------------------------------------------------------------------------------------------------
! Test that creates / fills quadrature_data in Fortan and passes to C++
!---------------------------------------------------------------------------------------------------
subroutine test_quadrature_interfaces() bind(c)

  use rtt_test_quadrature_f
  use quadrature_interfaces, only: quadrature_data, get_quadrature, init_quadrature
  use iso_c_binding, only: c_int, c_double, c_loc
  implicit none

  !----------------------------------------------------------------------------------------------
  ! Variable declarations
  type(quadrature_data)      :: quad
  integer(c_int)             :: error_code, sn_size
  real(c_double), allocatable, target, dimension(:) :: q_mu
  real(c_double), allocatable, target, dimension(:) :: q_eta
  real(c_double), allocatable, target, dimension(:) :: q_xi
  real(c_double), allocatable, target, dimension(:) :: q_wt

  !----------------------------------------------------------------------------------------------
  ! Initialization

  call init_quadrature(quad)

  ! Test a Tri Chebyshev Legendre quadrature
  quad % dimension = 2
  quad % type = 1
  quad % order = 4
  quad % geometry = 0
  ! Allocate space to fill in the data
  sn_size = 12
  allocate (q_mu(12), q_eta(12), q_xi(12), q_wt(12))
  quad % mu = c_loc(q_mu)
  quad % eta = c_loc(q_eta)
  quad % xi = c_loc(q_xi)
  quad % wt = c_loc(q_wt)

  ! Fill in the information
  call get_quadrature(quad)

  error_code = -1
  print '(a)', "On the Fortran side, we have created a quadrature_data type"
  print '(a,i1)', "The dimension is ", quad % dimension
  print '(a,i1)', "The type is ", quad % type
  print '(a,i2)', "The order is ", quad % order
  print '(a,i1)', "The geometry is ", quad % geometry
  print '(a,4f7.5)', "The first ordinate is ", q_mu(1), q_eta(1), q_xi(1), q_wt(1)
  print '(a)', "Now calling C to ensure data is passed correctly"
  print '(a)'

  !----------------------------------------------------------------------------------------------
  ! Call the c-function with the derived type and check the error code
  call rtt_test_quadrature_interfaces(quad, error_code)

  deallocate (q_mu, q_eta, q_xi, q_wt)

  if (error_code == 0) then
    print '(a)', "Test: passed"
    print '(a)', "     error code is equal to zero"
  else
    print '(a)', "Test: failed"
    print '(a,i5)', "     error code not equal to zero; it's ", error_code
  endif

  print '(a)', " "
  print '(a)', "*********************************************"
  if (error_code /= 0) then
    print '(a)', "**** ftstquadrature_interface Test: FAILED."
  else
    print '(a)', "**** ftstquadrature_interface Test: PASSED."
  endif
  print '(a)', "*********************************************"

end subroutine test_quadrature_interfaces

!---------------------------------------------------------------------------------------------------
! end of tstquadrature_interfaces.f90
!---------------------------------------------------------------------------------------------------
