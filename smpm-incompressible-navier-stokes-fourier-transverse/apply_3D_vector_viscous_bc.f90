subroutine apply_3D_vector_viscous_bc( q, Bq, q_x, q_z, cond, delta )
! Computes the misfit of the penalized boundary conditions in qx and qz for
! the vector viscous equation and adds the misfit to Bqx and Bqz.
!
! cond          - flag indicating Dirichlet or Neumann boundaries (see below).
! (q_x,q_y,q_z) - the gradient of q.
! delta         - the viscosity.

  use constants, only:             n, nprocs, nsuby, nsubz, pd, rank
  use mesh_deformation_maps, only: eta_x, eta_z, nx, nz, xi_x, xi_z
  use precision, only:             dp
  use woodbury_matrices, only:     rpk

  implicit none

  real(kind=dp), dimension(1:3*rpk, 1:nsuby), intent(in)    :: q
  real(kind=dp), dimension(1:3*rpk, 1:nsuby), intent(inout) :: Bq
  real(kind=dp), dimension(1:3*rpk, 1:nsuby), intent(in)    :: q_x
  real(kind=dp), dimension(1:3*rpk, 1:nsuby), intent(in)    :: q_z

  ! Boundary condition for each edge of the quad.  edges start at the bottom
  ! and move counter-clockwise (cond(1) = bottom, cond(2) = right, cond(3) = top,
  ! and cond(4) = left.  Dirichlet condition is equal to 1, while Neumann is 2.
  integer, dimension(4), intent(in)                         :: cond
  real(kind=dp), intent(in)                                 :: delta

  ! Internal variables.
  integer                                                   :: a, z, g, ii

  real(kind=dp), dimension(1:2*rpk, 1:nsuby)                :: qxz, Bqxz, qxz_x, qxz_z
  real(kind=dp), dimension(1:rpk, 1:nsuby)                  :: Bqy

  real(kind=dp), parameter                                  :: alpha = 1.0_dp
  ! Constants defining the penalty parameter for the Dirichlet and Neumann
  ! boundaries.
  !
  ! NOTE: these should be parameters, though Fortran does not allow them to be
  !       since they're derived from a module variable (pd).
  real(kind=dp)                                             :: omega, tau_d, tau_n

  ! NOTE:
  !
  !   q_x := x derivative of [ux, uy, uz] = [ux_x, uy_x, uz_x]
  !   q_y := y derivative of [ux, uy, uz] = [ux_y, uy_y, uz_y]
  !   q_z := z derivative of [ux, uy, uz] = [ux_z, uy_z, uz_z]

  ! Set up inputs for the 2D x/z vector viscous boundary conditions.

  ! Package the x velocity-related components we'll need.
  a              = 1
  z              = rpk
  Bqxz(a:z, :)   =  Bq(a:z, :)
  qxz(a:z, :)    =   q(a:z, :)
  qxz_x(a:z, :)  = q_x(a:z, :)
  qxz_z(a:z, :)  = q_z(a:z, :)

  ! Package the z velocity-related components we'll need (we'll skip the
  ! y-velocity components for now).
  a              = rpk + 1
  z              = 2 * rpk
  Bqxz(a:z, :)   =  Bq(a + rpk:z + rpk, :)
  qxz(a:z, :)    =   q(a + rpk:z + rpk, :)
  qxz_x(a:z, :)  = q_x(a + rpk:z + rpk, :)
  qxz_z(a:z, :)  = q_z(a + rpk:z + rpk, :)

  ! Compute the x/z velocity coupled vector viscous boundary conditions one
  ! plane at a time.
  do ii = 1, nsuby
     call apply_vector_viscous_bc( qxz(:, ii), Bqxz(:, ii), qxz_x(:, ii), qxz_z(:, ii), cond, delta )
  enddo

  ! Update the output with the x/z velocity boundary conditions.

  ! Update the x velocity-related components.
  a          = 1
  z          = rpk
  Bq(a:z, :) = Bq(a:z, :) + Bqxz(a:z, :) ! XXX: might be a minus.

  ! Update z velocity-related components.
  a                      = rpk + 1
  z                      = 2 * rpk
  Bq(a + rpk:z + rpk, :) = Bq(a + rpk:z + rpk, :) + Bqxz(a:z, :) ! XXX: might be a minus.

  ! Compute the boundary condition on the y-velocity.

  ! Extract the y-component of the operator (it is the middle third).
  Bqy = Bq(rpk + 1:2 * rpk, 1:nsuby)

  ! Set some constants.
  omega = 2.0_dp / (pd * (pd + 1.0_dp))
  tau_d = delta / alpha / omega**2
  tau_n = delta / omega

  ! Apply the y-velocity boundary condition one plane at a time.
  do ii = 1, nsuby

     ! Bottom boundary condition.
     if ( cond(1) == 1 ) then
        ! Do nothing, these surfaces have no boundary condition since the
        ! normal vector is orthogonal to the y-velocity.
     else

        ! Apply free-slip boundary condition operator.
        a            = 1
        z            = rpk
        g            = n * nsubz

        ! Compute the boundary condition.
        Bqy(a:z:g, ii) = Bqy(a:z:g, ii) - tau_n * hypot( xi_z(a:z:g), xi_x(a:z:g) ) * &
                         (nx(a:z:g, 1) * q_x(rpk + 1:2 * rpk, ii) + nz(a:z:g, 1) * q_z(rpk + 1:2 * rpk, ii))

     endif

     ! Right boundary condition.
     if ( rank == nprocs - 1 ) then
        if ( cond(2) == 1 ) then
           ! Do nothing, these surfaces have no boundary condition since the
           ! normal vector is orthogonal to the y-velocity.
        else

           ! Apply free-slip boundary condition operator.
           a       = rpk - (n * nsubz - 1)
           z       = rpk - 0
           g       = 1

           ! Apply the free-slip boundary condition operator in y.
           Bqy(a:z, ii) = Bqy(a:z, ii) - tau_n * hypot( eta_z(a:z), eta_x(a:z) ) * &
                          (nx(a:z, 2) * q_x(rpk + 1:2 * rpk, ii) + nz(a:z, 2) * q_z(rpk + 1:2 * rpk, ii))

        endif
     endif

     ! Top boundary condition.
     if ( cond(3) == 1 ) then
        ! Do nothing, these surfaces have no boundary condition since the
        ! normal vector is orthogonal to the y-velocity.
     else

        ! Apply the free-slip  boundary condition operator.
        a           = n * nsubz
        z           = rpk
        g           = n * nsubz

        ! Compute the boundary condition.
        Bqy(a:z:g, ii) = Bqy(a:z:g, ii) - tau_n * hypot( xi_z(a:z:g), xi_x(a:z:g) ) * &
                         (nx(a:z:g, 3) * q_x(rpk + 1:2 * rpk, ii) + nz(a:z:g, 3) * q_z(rpk + 1:2 * rpk, ii))

     endif

     ! Left boundary condition.
     if ( rank == 0 ) then
        if ( cond(4) == 1 ) then
           ! Do nothing.
        else

           ! Apply free-slip boundary condition operator.
           a       = 1
           z       = n * nsubz
           g       = 1

           ! Apply the free-slip boundary condition operator in y.
           Bqy(a:z, ii) = Bqy(a:z, ii) - tau_n * hypot( eta_z(a:z), eta_x(a:z) ) * &
                          (nx(a:z, 4) * q_x(rpk + 1:2 * rpk, ii) + nz(a:z, 4) * q_z(rpk + 1:2 * rpk, ii))

        endif
     endif
  enddo

  ! Return the y-velocity operator component into the output vector.
  Bq(rpk+1:2*rpk, 1:nsuby) = Bqy

end subroutine apply_3D_vector_viscous_bc
