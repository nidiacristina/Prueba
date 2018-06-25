subroutine set_timestepping_coefficients( timestep, dt, dt1, dt2 )
! Updates the timestepping coefficients based on the number of timesteps
! simulated thus far.

  use options, only:                 apply_restart
  use precision, only:               dp
  use timestepping, only:            a0, a1, a2, b0, b1, b2, g0, &
                                     get_a0, get_a1, get_a2, get_b0, get_b1, get_b2, get_g0

  implicit none

  integer, intent(in)       :: timestep
  real(kind=dp), intent(in) :: dt
  real(kind=dp), intent(in) :: dt1
  real(kind=dp), intent(in) :: dt2

  ! compute our timestep coefficients from the most recent step sizes.
  g0 = get_g0( dt, dt1, dt2 )
  a0 = get_a0( dt, dt1, dt2 )
  a1 = get_a1( dt, dt1, dt2 )
  a2 = get_a2( dt, dt1, dt2 )
  b0 = get_b0( dt, dt1, dt2 )
  b1 = get_b1( dt, dt1, dt2 )
  b2 = get_b2( dt, dt1, dt2 )

  ! the first three absolute timesteps, not relative to a restart, in our
  ! simulation have fixed coefficients rather than computed.
  if ( apply_restart .eqv. .false. .and. timestep <= 3 ) then
     if ( timestep == 1 ) then

        ! First order method.
        a0 = 1.0_dp
        a1 = 0.0_dp
        a2 = 0.0_dp
        b0 = 1.0_dp
        b1 = 0.0_dp
        b2 = 0.0_dp
        g0 = 1.0_dp

     else if ( timestep == 2 ) then

        ! Second order method.
        a0 = 2.0_dp
        a1 = -0.5_dp
        a2 = 0.0_dp
        b0 = 2.0_dp
        b1 = -1.0_dp
        b2 = 0.0_dp
        g0 = 1.5_dp

     else if ( timestep == 3 ) then

        ! Third order method.
        a0 = 3.0_dp
        a1 = -1.5_dp
        a2 = 1.0_dp / 3.0_dp
        b0 = 3.0_dp
        b1 = -3.0_dp
        b2 = 1.0_dp
        g0 = 11.0_dp / 6.0_dp

     endif
  endif

end subroutine set_timestepping_coefficients
