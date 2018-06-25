module timestepping
! Contains the constants for the Adams-Bashforth time-stepping method as well
! as constants related to the time-step.

  use precision, only: dp

  implicit none
  save

  integer       :: time_ndx                        ! Current time-step index.
  real(kind=dp) :: dt                              ! Size of the time step.
  real(kind=dp) :: dt0, dt1, dt2                   ! Current time-step, previous time-step, previous-previous time-step.
  real(kind=dp) :: tend                            ! Final time (physical).
  real(kind=dp) :: cfl_number
  real(kind=dp) :: cfl_max_x, cfl_max_y, cfl_max_z
  real(kind=dp) :: simulation_time
  logical       :: logflag                         ! Global parameter that establishes whether
                                                   ! to write out diagnostic information or not.

  ! The time-stepping constants as in Karniadakis-1991.
  real(kind=dp) :: g0 = 11.0_dp/6.0_dp
  real(kind=dp) :: a0 = 3.0_dp
  real(kind=dp) :: a1 = -3.0_dp/2.0_dp
  real(kind=dp) :: a2 = 1.0_dp/3.0_dp
  real(kind=dp) :: b0 = 3.0_dp
  real(kind=dp) :: b1 = -3.0_dp
  real(kind=dp) :: b2 = 1.0_dp

  contains

  ! Below are functions for the adaptive time-stepping.  These all return the
  ! values of the above coefficients based on the past time-steps and the new
  ! desired timestep.

  ! In all functions below:
  !
  !      dt0 = current time step.
  !      dt1 = previous time step.
  !      dt2 = previous previous time step.

  function get_g0( in_dt0, in_dt1, in_dt2 ) result( g0 )

     use precision, only: dp

     implicit none

     real(kind=dp), intent(in) :: in_dt0
     real(kind=dp), intent(in) :: in_dt1
     real(kind=dp), intent(in) :: in_dt2
     real(kind=dp)             :: g0

     g0 = 1.0_dp + 1.0_dp / (1 + in_dt1 / in_dt0) + 1.0_dp / (1.0_dp + in_dt1 / in_dt0 + in_dt2 / in_dt0)

  end function get_g0

  function get_a0( in_dt0, in_dt1, in_dt2 ) result( a0 )

     use precision, only: dp

     implicit none

     real(kind=dp), intent(in) :: in_dt0
     real(kind=dp), intent(in) :: in_dt1
     real(kind=dp), intent(in) :: in_dt2
     real(kind=dp)             :: a0

     real(kind=dp)             :: dt1, dt2

     dt1 = in_dt1 / in_dt0
     dt2 = in_dt2 / in_dt0
     a0  = (1.0_dp + dt1) * (1.0_dp + dt1 + dt2) / (dt1 * (dt1 + dt2))

  end function get_a0


  function get_a1( in_dt0, in_dt1, in_dt2 ) result( a1 )

     use precision, only: dp

     implicit none

     real(kind=dp), intent(in) :: in_dt0
     real(kind=dp), intent(in) :: in_dt1
     real(kind=dp), intent(in) :: in_dt2
     real(kind=dp)             :: a1

     real(kind=dp)             :: dt1, dt2

     dt1 = in_dt1 / in_dt0
     dt2 = in_dt2 / in_dt0
     a1  = -(1.0_dp + dt1 + dt2) / (dt1 * dt2 * (1.0_dp + dt1))

  end function get_a1


  function get_a2( in_dt0, in_dt1, in_dt2 ) result( a2 )

     use precision, only: dp

     implicit none

     real(kind=dp), intent(in) :: in_dt0
     real(kind=dp), intent(in) :: in_dt1
     real(kind=dp), intent(in) :: in_dt2
     real(kind=dp)             :: a2

     real(kind=dp)             :: dt1, dt2

     dt1 = in_dt1 / in_dt0
     dt2 = in_dt2 / in_dt0
     a2  = (1.0_dp + dt1) / (dt1 * (dt1 + dt2) * (1.0_dp + dt1 + dt2))

  end function get_a2


  function get_b0( in_dt0, in_dt1, in_dt2 ) result( b0 )

     use precision, only: dp

     implicit none

     real(kind=dp), intent(in) :: in_dt0
     real(kind=dp), intent(in) :: in_dt1
     real(kind=dp), intent(in) :: in_dt2
     real(kind=dp)             :: b0

     real(kind=dp)             :: dt1, dt2

     dt1 = in_dt1 / in_dt0
     dt2 = in_dt2 / in_dt0
     b0  = (1.0_dp + dt1) * (1.0_dp + dt1 + dt2) / (dt1 * (dt1 + dt2))

  end function get_b0


  function get_b1( in_dt0, in_dt1, in_dt2 ) result( b1 )

     use precision, only: dp

     implicit none

     real(kind=dp), intent(in) :: in_dt0
     real(kind=dp), intent(in) :: in_dt1
     real(kind=dp), intent(in) :: in_dt2
     real(kind=dp)             :: b1

     real(kind=dp)             :: dt1, dt2

     dt1 = in_dt1 / in_dt0
     dt2 = in_dt2 / in_dt0
     b1  = -(1.0_dp + dt1 + dt2) / (dt1 * dt2)

  end function get_b1


  function get_b2( in_dt0, in_dt1, in_dt2 ) result( b2 )

     use precision, only: dp

     implicit none

     real(kind=dp), intent(in) :: in_dt0
     real(kind=dp), intent(in) :: in_dt1
     real(kind=dp), intent(in) :: in_dt2
     real(kind=dp)             :: b2

     real(kind=dp)             :: dt1, dt2

     dt1 = in_dt1 / in_dt0
     dt2 = in_dt2 / in_dt0
     b2  = (1.0_dp + dt1) / (dt1 * (dt1 + dt2))

  end function get_b2

end module timestepping
