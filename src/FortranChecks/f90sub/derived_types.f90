!---------------------------------------------*-F90-*---------------------------------------------!
!
! file   FortranCheck/ftest/derived_types.f90
! author Kelly Thompson, Allan Wollaber
! date   Tuesday, Jun 12, 2012, 16:03 pm
! brief  Test F90 derived types passed into a C++ function
! note   Copyright (C) 2014-2022 Triad National Security, LLC., All rights reserved.
!-------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------!
! Create a module that contains a derived type and an interface to a C function that uses it.
! Note the "import" statement brings the type into local scope.
!-------------------------------------------------------------------------------------------------!

module rtt_test_derived_types
  use iso_c_binding, only: c_double, c_int, c_int64_t, c_ptr, c_null_ptr
  implicit none
  ! Proper enumeration types do not work in F2003, but this is close.  Enumerations with bind(C) are
  ! just integers anyway.  See:
  ! http://www.rhinocerus.net/forum/lang-fortran/92750-no-enumeration-data-type-fortran.html
  integer(c_int), parameter :: GREY = 0, MULTIGROUP = 1, ODF = 2

  ! Create a derived type that contains some "extra" information
  type, bind(C)  :: my_informative_type
    real(c_double)       :: some_double
    integer(c_int)       :: some_int
    integer(c_int64_t)   :: some_large_int
    type(c_ptr)          :: some_pointer = C_NULL_PTR
    integer(c_int)       :: some_enum
  end type my_informative_type

  ! Now create an interface to the C routine that accepts that type
  interface
    subroutine rtt_test_derived_type(dt_in, err_code) bind(C, &
                                                           name="rtt_test_derived_type")
      use iso_c_binding, only: c_int
      import my_informative_type  ! F2003 standard, brings in host scope
      implicit none
      type(my_informative_type), intent(in)  :: dt_in
      integer(c_int), intent(out) :: err_code
    end subroutine rtt_test_derived_type
  end interface

end module rtt_test_derived_types

!---------------------------------------------------------------------------
! Fortran main that uses the above module type and bound c-interface
!---------------------------------------------------------------------------
subroutine test_derived_types() bind(c)

  use rtt_test_derived_types
  use iso_c_binding, only: c_int, c_int64_t, c_double, c_loc, c_ptr
  implicit none

  !----------------------------------------------------------------------
  ! Variable declarations
  type(my_informative_type)  :: mit
  integer(c_int)             :: error_code
  integer(c_int), allocatable, target, dimension(:) :: int_array

  !----------------------------------------------------------------------
  ! Initialization

  mit % some_double = 3.141592654_c_double
  mit % some_int = 137
  mit % some_large_int = 2_c_int64_t**34
  mit % some_enum = MULTIGROUP

  allocate (int_array(10))
  int_array(1) = 2003
  int_array(2) = 2012
  mit % some_pointer = c_loc(int_array)

  error_code = -1

  print '(a,f7.5)', "On Fortran side, derived type contains double = ", mit % some_double
  print '(a,i3)', "integer = ", mit % some_int
  print '(a,i11)', "large integer = ", mit % some_large_int
  print '(a,i4)', "int_array(1) = ", int_array(1)
  print '(a,i4)', "int_array(2) = ", int_array(2)
  print '(a,i1)', "The enumerated type is MULTIGROUP = ", MULTIGROUP
  print '(a)', " "

  !----------------------------------------------------------------------
  ! Call the c-function with the derived type and check the error code
  call rtt_test_derived_type(mit, error_code)

  deallocate (int_array)

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
    print '(a)', "**** test_derived_types Test: FAILED."
  else
    print '(a)', "**** test_derived_types Test: PASSED."
  endif
  print '(a)', "*********************************************"

end subroutine test_derived_types

!---------------------------------------------------------------------------
! end of derived_types.f90
!---------------------------------------------------------------------------
