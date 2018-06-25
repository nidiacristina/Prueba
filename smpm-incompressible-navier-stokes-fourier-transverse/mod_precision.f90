module precision
! Contains kinds for correctly selecting a data type with a specific level of
! precision.

  implicit none
  save

  ! Kind representing a 64-bit floating point value.  This is used to specify
  ! the precision of floating point values regardless of whether they're real or
  ! complex.
  !
  ! XXX: this should use the F2008 features to select 64-bit types directly.
  integer, parameter          :: dp = selected_real_kind( 15, 307 )

  ! Kind representing a 64-bit 2's complement, signed integer.
  !
  ! XXX: this should use the F2008 features to select 64-bit types directly.
  integer, parameter          :: int64 = selected_int_kind( 18 )

  ! NaN for debugging purposes.  This is provided for debugging purposes, and is
  ! not needed for correct solver execution to help identify errant behavior
  ! (e.g. computation that uses uninitialized memory can have said memory
  ! initialized to NaN to confirm the problem).
  !
  ! NOTE: This assumes execution on an IEEE-754 compliant system as the hard
  !       coded constant is a non-signaling NaN as interpreted as a 64-bit
  !       integer.
  real(kind=dp), parameter    :: NaN_value = transfer( -2251799813685248_int64, &
                                                       1.0_dp )

end module precision
