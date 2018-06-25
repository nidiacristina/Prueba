subroutine apply_vector_viscous_bc( q, Bq, q_x, q_z, cond, delta )
! Computes the misfit of the penalized boundary conditions in qx and qz for
! the vector viscous equation and adds the misfit to Bqx and Bqz.
!
! cond      - flag indicating Dirichlet or Neumann boundaries (see below).
! (q_x,q_z) - the gradient of q.
! delta     - the viscosity.

  use constants, only:             n, nprocs, nsg, nsubz, pd, rank
  use mesh_deformation_maps, only: eta_x, eta_z, nx, nz, tx, tz, &
                                   tx_x, tx_z, tz_x, tz_z, xi_x, xi_z
  use precision, only:             dp
  use woodbury_matrices, only:     rpk

  implicit none

  real(kind=dp), dimension(1:2*rpk), intent(in)    :: q
  real(kind=dp), dimension(1:2*rpk), intent(inout) :: Bq
  real(kind=dp), dimension(1:2*rpk), intent(in)    :: q_x
  real(kind=dp), dimension(1:2*rpk), intent(in)    :: q_z

  ! boundary condition for each edge of the quad.  edges start at the bottom
  ! and move counter-clockwise (cond(1) = bottom, cond(2) = right, cond(3) = top,
  ! and cond(4) = left.  Dirichlet condition is equal to 1, while Neumann is 2.
  integer, dimension(4), intent(in)                :: cond
  real(kind=dp), intent(in)                        :: delta

  ! Internal variables.
  real(kind=dp), dimension(1:rpk)                  :: qx, qz, Bqx, Bqz
  real(kind=dp), dimension(1:rpk)                  :: qx_x, qx_z, qz_x, qz_z
  integer                                          :: a, z, g
  real(kind=dp), dimension(1:nsg)                  :: rt_x, rt_z, rnx, rnz, qt, qn, qt_n

  real(kind=dp), parameter                         :: alpha = 1.0_dp

  ! constants defining the penalty parameter for the Dirichlet and Neumann
  ! boundaries.
  !
  ! NOTE: these should be parameters, though Fortran does not allow them to be
  !       since they're derived from a module variable (pd).
  real(kind=dp)                                    :: omega, tau_d, tau_n

  ! Split the state vectors into x- and z-velocities.
  qx   =   q(      1:  rpk)
  qz   =   q(rpk + 1:2*rpk)
  qx_x = q_x(      1:  rpk)
  qz_x = q_x(rpk + 1:2*rpk)
  qx_z = q_z(      1:  rpk)
  qz_z = q_z(rpk + 1:2*rpk)
  Bqx  =  Bq(      1:  rpk)
  Bqz  =  Bq(rpk + 1:2*rpk)

  ! Set some constants.
  omega = 2.0_dp / (pd * (pd + 1.0_dp))
  tau_d = delta / alpha / omega**2
  tau_n = delta / omega

  ! Bottom boundary condition.
  if ( cond(1) == 1 ) then
     ! Apply no-slip boundary condition operator.
     a          = 1
     z          = rpk
     g          = n * nsubz
     Bqx(a:z:g) = Bqx(a:z:g) - tau_d * (xi_z(a:z:g) + xi_x(a:z:g))**2 * qx(a:z:g)
     Bqz(a:z:g) = Bqz(a:z:g) - tau_d * (xi_z(a:z:g) + xi_x(a:z:g))**2 * qz(a:z:g)
  else

     ! Apply free-slip boundary condition operator.
     a            = 1
     z            = rpk
     g            = n * nsubz

     ! Compute the normal and tangential components of velocity and the normal
     ! derivative of the tangential component.
       qn(a:z:g) = nx(a:z:g, 1) * qx(a:z:g) + nz(a:z:g, 1) * qz(a:z:g)
       qt(a:z:g) = tx(a:z:g, 1) * qx(a:z:g) + tz(a:z:g, 1) * qz(a:z:g)
     qt_n(a:z:g) = nx(a:z:g, 1) * (tx(a:z:g, 1) * qx_x(a:z:g) + tx_x(a:z:g, 1) * qx(a:z:g) + &
                                   tz(a:z:g, 1) * qz_x(a:z:g) + tz_x(a:z:g, 1) * qz(a:z:g)) + &
                   nz(a:z:g, 1) * (tx(a:z:g, 1) * qx_z(a:z:g) + tx_z(a:z:g, 1) * qx(a:z:g) + &
                                   tz(a:z:g, 1) * qz_z(a:z:g) + tz_z(a:z:g, 1) * qz(a:z:g))
                                !
                                ! This above expression in math is this:
                                !
                                !  nx * Dx( <t,u> ) + nz * Dz( <t,u> )
                                !
                                ! where <a,b> denotes the dot product of two vectors,
                                ! and Dx and Dz are derivatives in x and z.

     ! Compute the residual in (x,z) coordinates.
     rnx(a:z:g)  = nx(a:z:g, 1) * qn(a:z:g)
     rnz(a:z:g)  = nz(a:z:g, 1) * qn(a:z:g)
     rt_x(a:z:g) = tx(a:z:g, 1) * qt_n(a:z:g)
     rt_z(a:z:g) = tz(a:z:g, 1) * qt_n(a:z:g)

     ! Apply the boundary condition.
     Bqx(a:z:g) = Bqx(a:z:g) - tau_d * hypot( xi_z(a:z:g), xi_x(a:z:g) )**2 *  rnx(a:z:g) &
                             - tau_n * hypot( xi_z(a:z:g), xi_x(a:z:g) )    * rt_x(a:z:g)

     Bqz(a:z:g) = Bqz(a:z:g) - tau_d * hypot( xi_z(a:z:g), xi_x(a:z:g) )**2 *  rnz(a:z:g) &
                             - tau_n * hypot( xi_z(a:z:g), xi_x(a:z:g) )    * rt_z(a:z:g)

     ! XXX: Does this reduce to what I expect for the case of undeformed
     !      boundaries?
  endif

  ! Right boundary condition.
  if ( rank == nprocs - 1 ) then
     if ( cond(2) == 1 ) then

        ! Apply no-slip boundary condition operator.
        a        = rpk - (n * nsubz - 1)
        z        = rpk - 0
        Bqx(a:z) = Bqx(a:z) - tau_d * hypot( eta_x(a:z), eta_z(a:z) )**2 * qx(a:z)
        Bqz(a:z) = Bqz(a:z) - tau_d * hypot( eta_x(a:z), eta_z(a:z) )**2 * qz(a:z)

     else

        ! Apply free-slip boundary condition operator.
        a       = rpk - (n * nsubz - 1)
        z       = rpk - 0
        g       = 1

        ! Compute the normal and tangential components of velocity and the
        ! normal derivative of the tangential component.
          qn(a:z:g) = nx(a:z:g, 2) * qx(a:z:g) + nz(a:z:g, 2) * qz(a:z:g)
          qt(a:z:g) = tx(a:z:g, 2) * qx(a:z:g) + tz(a:z:g, 2) * qz(a:z:g)
        qt_n(a:z:g) = nx(a:z:g, 2) * (tx(a:z:g, 2) * qx_x(a:z:g) + tx_x(a:z:g, 2) * qx(a:z:g) + &
                                      tz(a:z:g, 2) * qz_x(a:z:g) + tz_x(a:z:g, 2) * qz(a:z:g)) + &
                      nz(a:z:g, 2) * (tx(a:z:g, 2) * qx_z(a:z:g) + tx_z(a:z:g, 2) * qx(a:z:g) + &
                                      tz(a:z:g, 2) * qz_z(a:z:g) + tz_z(a:z:g, 2) * qz(a:z:g))
                                   !
                                   ! This above expression in math is this:
                                   !
                                   !  nx * Dx( <t,u> ) + nz * Dz( <t,u> )
                                   !
                                   ! where <a,b> denotes the dot product of two vectors,
                                   ! and Dx and Dz are derivatives in x and z.

        ! Compute the residual in (x,z) coordinates.
        rnx(a:z:g)  = nx(a:z:g, 2) * qn(a:z:g)
        rnz(a:z:g)  = nz(a:z:g, 2) * qn(a:z:g)
        rt_x(a:z:g) = tx(a:z:g, 2) * qt_n(a:z:g)
        rt_z(a:z:g) = tz(a:z:g, 2) * qt_n(a:z:g)

        ! Apply the boundary condition.
        Bqx(a:z:g) = Bqx(a:z:g) - tau_d * hypot(eta_z(a:z:g), eta_x(a:z:g))**2 *  rnx(a:z:g) &
                                - tau_n * hypot(eta_z(a:z:g), eta_x(a:z:g))    * rt_x(a:z:g)

        Bqz(a:z:g) = Bqz(a:z:g) - tau_d * hypot(eta_z(a:z:g), eta_x(a:z:g))**2 *  rnz(a:z:g) &
                                - tau_n * hypot(eta_z(a:z:g), eta_x(a:z:g))    * rt_z(a:z:g)
     endif
  endif

  ! Top boundary condition.
  if ( cond(3) == 1 ) then

     ! Apply the no-slip boundary condition operator.
     a          = n * nsubz
     z          = rpk
     g          = n * nsubz
     Bqx(a:z:g) = Bqx(a:z:g) - tau_d * hypot( xi_x(a:z:g), xi_z(a:z:g) )**2 * qx(a:z:g)
     Bqz(a:z:g) = Bqz(a:z:g) - tau_d * hypot( xi_x(a:z:g), xi_z(a:z:g) )**2 * qz(a:z:g)

  else

     ! Apply the free-slip  boundary condition operator.
     a           = n * nsubz
     z           = rpk
     g           = n * nsubz

     ! Compute the normal and tangential components of velocity and the normal
     ! derivative of the tangential component.
       qn(a:z:g) = nx(a:z:g, 3) * qx(a:z:g) + nz(a:z:g, 3) * qz(a:z:g)
       qt(a:z:g) = tx(a:z:g, 3) * qx(a:z:g) + tz(a:z:g, 3) * qz(a:z:g)
     qt_n(a:z:g) = nx(a:z:g, 3) * (tx(a:z:g, 3) * qx_x(a:z:g) + tx_x(a:z:g, 3) * qx(a:z:g) + &
                                   tz(a:z:g, 3) * qz_x(a:z:g) + tz_x(a:z:g, 3) * qz(a:z:g) ) + &
                   nz(a:z:g, 3) * (tx(a:z:g, 3) * qx_z(a:z:g) + tx_z(a:z:g, 3) * qx(a:z:g) + &
                                   tz(a:z:g, 3) * qz_z(a:z:g) + tz_z(a:z:g, 3) * qz(a:z:g) )
                                !
                                ! This above expression in math is this:
                                !
                                !  nx * Dx( <t,u> ) + nz * Dz( <t,u> )
                                !
                                ! where <a,b> denotes the dot product of two vectors,
                                ! and Dx and Dz are derivatives in x and z.

     ! Compute the residual in (x,z) coordinates.
     rnx(a:z:g)  = nx(a:z:g, 3) * qn(a:z:g)
     rnz(a:z:g)  = nz(a:z:g, 3) * qn(a:z:g)
     rt_x(a:z:g) = tx(a:z:g, 3) * qt_n(a:z:g)
     rt_z(a:z:g) = tz(a:z:g, 3) * qt_n(a:z:g)

     ! Apply the boundary condition.
     Bqx(a:z:g) = Bqx(a:z:g) - tau_d * hypot( xi_z(a:z:g), xi_x(a:z:g) )**2 *  rnx(a:z:g) &
                             - tau_n * hypot( xi_z(a:z:g), xi_x(a:z:g) )    * rt_x(a:z:g)

     Bqz(a:z:g) = Bqz(a:z:g) - tau_d * hypot( xi_z(a:z:g), xi_x(a:z:g) )**2 *  rnz(a:z:g) &
                             - tau_n * hypot( xi_z(a:z:g), xi_x(a:z:g) )    * rt_z(a:z:g)

  endif

  ! Left boundary condition.
  if ( rank == 0 ) then
     if ( cond(4) == 1 ) then

        ! Apply no-slip boundary condition operator.
        a        = 1
        z        = n * nsubz
        Bqx(a:z) = Bqx(a:z) - tau_d * hypot( eta_x(a:z), eta_z(a:z) )**2 * qx(a:z)
        Bqz(a:z) = Bqz(a:z) - tau_d * hypot( eta_x(a:z), eta_z(a:z) )**2 * qz(a:z)

     else

        ! Apply free-slip boundary condition operator.
        a       = 1
        z       = n * nsubz
        g       = 1

        ! Compute the normal and tangential components of velocity and the
        ! normal derivative of the tangential component.
          qn(a:z:g) = nx(a:z:g, 4) * qx(a:z:g) + nz(a:z:g, 4) * qz(a:z:g)
          qt(a:z:g) = tx(a:z:g, 4) * qx(a:z:g) + tz(a:z:g, 4) * qz(a:z:g)
        qt_n(a:z:g) = nx(a:z:g, 4) * (tx(a:z:g, 4) * qx_x(a:z:g) + tx_x(a:z:g, 4) * qx(a:z:g) + &
                                      tz(a:z:g, 4) * qz_x(a:z:g) + tz_x(a:z:g, 4) * qz(a:z:g)) + &
                      nz(a:z:g, 4) * (tx(a:z:g, 4) * qx_z(a:z:g) + tx_z(a:z:g, 4) * qx(a:z:g) + &
                                      tz(a:z:g, 4) * qz_z(a:z:g) + tz_z(a:z:g, 4) * qz(a:z:g))
                                   !
                                   ! This above expression in math is this:
                                   !
                                   !  nx * Dx( <t,u> ) + nz * Dz( <t,u> )
                                   !
                                   ! where <a,b> denotes the dot product of two vectors,
                                   ! and Dx and Dz are derivatives in x and z.

        ! Compute the residual in (x,z) coordinates.
        rnx(a:z:g)  = nx(a:z:g, 4) * qn(a:z:g)
        rnz(a:z:g)  = nz(a:z:g, 4) * qn(a:z:g)
        rt_x(a:z:g) = tx(a:z:g, 4) * qt_n(a:z:g)
        rt_z(a:z:g) = tz(a:z:g, 4) * qt_n(a:z:g)

        ! Apply the boundary condition.
        Bqx(a:z:g) = Bqx(a:z:g) - tau_d * hypot( eta_z(a:z:g), eta_x(a:z:g) )**2 *  rnx(a:z:g) &
                                - tau_n * hypot( eta_z(a:z:g), eta_x(a:z:g) )    * rt_x(a:z:g)

        Bqz(a:z:g) = Bqz(a:z:g) - tau_d * hypot( eta_z(a:z:g), eta_x(a:z:g) )**2 *  rnz(a:z:g) &
                                - tau_n * hypot( eta_z(a:z:g), eta_x(a:z:g) )    * rt_z(a:z:g)
     endif
  endif

  ! Reassemble components into a state vector.
  Bq(    1:  rpk) = Bqx
  Bq(rpk+1:2*rpk) = Bqz

end subroutine apply_vector_viscous_bc
