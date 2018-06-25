subroutine advance_time()
! Advances time variables in the SMPM main loop.

  use field_variables, only:         Cux0, Cux1, Cux2, Cuy0, Cuy1, Cuy2, Cuz0, Cuz1, Cuz2, &
                                     Nrho0, Nrho1, Nrho2, &
                                     Nux0, Nux1, Nux2, Nuy0, Nuy1, Nuy2, Nuz0, Nuz1, Nuz2, &
                                     rho, rho0, rho1, rho2, &
                                     ux, ux0, ux1, ux2, uy, uy0, uy1, uy2, uz, uz0, uz1, uz2
  use timestepping, only:            dt, dt1, dt2

  implicit none

  ! Advance timestep
  dt2 = dt1
  dt1 = dt

  ! Advance all arrays
  ux2   = ux1
  ux1   = ux0
  ux0   = ux
  uy2   = uy1
  uy1   = uy0
  uy0   = uy
  uz2   = uz1
  uz1   = uz0
  uz0   = uz
  rho2  = rho1
  rho1  = rho0
  rho0  = rho

  Nux2  = Nux1
  Nux1  = Nux0
  Nuy2  = Nuy1
  Nuy1  = Nuy0
  Nuz2  = Nuz1
  Nuz1  = Nuz0

  Cux2  = Cux1
  Cux1  = Cux0
  Cuy2  = Cuy1
  Cuy1  = Cuy0
  Cuz2  = Cuz1
  Cuz1  = Cuz0

  Nrho2 = Nrho1
  Nrho1 = Nrho0

end subroutine advance_time
