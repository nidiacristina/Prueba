subroutine apply_bc( q, Bq, cond, delta, q_x, q_z )
! Computes the misfit of the penalized boundary conditions in q and adds the
! misfit to Bq.
!
! cond      - flag indicating Dirichlet or Neumann boundaries (see below).
! (q_x,q_z) - the gradient of q.
! delta     - the viscosity.

  use constants, only:             n, nprocs, nsubz, pd, rank
  use mesh_deformation_maps, only: eta_x, eta_z, x_eta_n, x_xi_n, xi_x, xi_z, z_eta_n, z_xi_n
  use precision, only:             dp
  use woodbury_matrices, only:     rpk

  implicit none

  real(kind=dp), dimension(rpk), intent(in)    :: q
  real(kind=dp), dimension(rpk), intent(inout) :: Bq

  ! Boundary condition for each edge of the quad.  edges start at the bottom
  ! and move counter-clockwise (cond(1) = bottom, cond(2) = right, cond(3) = top,
  ! and cond(4) = left.  Dirichlet condition is equal to 1, while Neumann is 2.
  integer, dimension(4), intent(in)            :: cond
  real(kind=dp), intent(in)                    :: delta
  real(kind=dp), dimension(rpk), intent(in)    :: q_x
  real(kind=dp), dimension(rpk), intent(in)    :: q_z

  integer                                      :: a, z, gap

  real(kind=dp), parameter                     :: alpha = 1.0_dp
  ! Constants defining the penalty parameter for the Dirichlet and Neumann
  ! boundaries.
  !
  ! NOTE: these should be parameters, though Fortran does not allow them to be
  !       since they're derived from a module variable (pd).
  real(kind=dp)                                :: omega, tau_d, tau_n

  ! Set some constants.
  omega = 2.0_dp / (pd * (pd + 1.0_dp))
  tau_d = delta / alpha / omega**2
  tau_n = delta / omega

  ! Bottom boundary condition.
  if ( cond(1) == 1 ) then
     ! Apply Dirichlet boundary condition operator.
     a           = 1
     z           = rpk
     gap         = n * nsubz
     !Bq(a:z:gap) = Bq(a:z:gap) - tau_d * (xi_z(a:z:gap) + xi_x(a:z:gap))**2 * q(a:z:gap)
     Bq(a:z:gap) = Bq(a:z:gap) - tau_d * hypot( xi_z(a:z:gap), xi_x(a:z:gap) )**2 * q(a:z:gap)
  else
     ! Apply Neumann boundary condition operator.
     a           = 1
     z           = rpk
     gap         = n * nsubz
!     Bq(a:z:gap) = Bq(a:z:gap) - &
!                   tau_n * xi_z(a:z:gap) * ( z_eta(a:z:gap) * q_x(a:z:gap)) / sqrt( z_eta(a:z:gap)**2 + x_eta(a:z:gap)**2 ) - &
!                   tau_n * xi_z(a:z:gap) * (-x_eta(a:z:gap) * q_z(a:z:gap)) / sqrt( z_eta(a:z:gap)**2 + x_eta(a:z:gap)**2 ) - &
!                   tau_n * xi_x(a:z:gap) * ( z_eta(a:z:gap) * q_x(a:z:gap)) / sqrt( z_eta(a:z:gap)**2 + x_eta(a:z:gap)**2 ) - &
!                   tau_n * xi_x(a:z:gap) * (-x_eta(a:z:gap) * q_z(a:z:gap)) / sqrt( z_eta(a:z:gap)**2 + x_eta(a:z:gap)**2 )
!                   !     penalty coef.   ! normal vector    ! D( field )    ! normalization factor for normal vector
     Bq(a:z:gap) = Bq(a:z:gap) - hypot( xi_z(a:z:gap), xi_x(a:z:gap) ) * &
                   tau_n * (z_eta_n(a:z:gap) * q_x(a:z:gap) - x_eta_n(a:z:gap) * q_z(a:z:gap))
                   !     penalty coef.   ! normal vector    ! D( field )    ! normalization factor for normal vector
  endif

  ! Right boundary condition.
  if ( rank == nprocs - 1 ) then
     if ( cond(2) == 1 ) then
        ! Apply Dirichlet boundary condition operator.
        a       = rpk - (n * nsubz - 1)
        z       = rpk - 0
        !Bq(a:z) = Bq(a:z) - tau_d * (eta_x(a:z) + eta_z(a:z))**2 * q(a:z)
        Bq(a:z) = Bq(a:z) - tau_d * hypot( eta_x(a:z), eta_z(a:z) )**2 * q(a:z)
     else
        ! Apply Neumann boundary condition operator.
        a       = rpk - (n * nsubz - 1)
        z       = rpk - 0
!        Bq(a:z) = Bq(a:z) - &
!                  tau_n * eta_x(a:z) * ( z_xi(a:z) * q_x(a:z)) / sqrt( z_xi(a:z)**2 + x_xi(a:z)**2 ) - &
!                  tau_n * eta_x(a:z) * (-x_xi(a:z) * q_z(a:z)) / sqrt( z_xi(a:z)**2 + x_xi(a:z)**2 ) - &
!                  tau_n * eta_z(a:z) * ( z_xi(a:z) * q_x(a:z)) / sqrt( z_xi(a:z)**2 + x_xi(a:z)**2 ) - &
!                  tau_n * eta_z(a:z) * (-x_xi(a:z) * q_z(a:z)) / sqrt( z_xi(a:z)**2 + x_xi(a:z)**2 )
        Bq(a:z) = Bq(a:z) - hypot( eta_x(a:z), eta_z(a:z) ) * &
                  tau_n * (z_xi_n(a:z) * q_x(a:z) - x_xi_n(a:z) * q_z(a:z))
     endif
  endif

  ! Top boundary condition.
  if ( cond(3) == 1 ) then
     ! Apply Dirichlet boundary condition operator.
     a           = n * nsubz
     z           = rpk
     gap         = n * nsubz
     !Bq(a:z:gap) = Bq(a:z:gap) - tau_d * (xi_x(a:z:gap) + xi_z(a:z:gap))**2 * q(a:z:gap)
     Bq(a:z:gap) = Bq(a:z:gap) - tau_d * hypot( xi_x(a:z:gap), xi_z(a:z:gap) )**2 * q(a:z:gap)
  else
     ! Apply Neumann boundary condition operator.
     a           = n * nsubz
     z           = rpk
     gap         = n * nsubz
!     Bq(a:z:gap) = Bq(a:z:gap) - &
!                   tau_n * xi_z(a:z:gap) * (-z_eta(a:z:gap) * q_x(a:z:gap)) / sqrt( z_eta(a:z:gap)**2 + x_eta(a:z:gap)**2 ) - &
!                   tau_n * xi_z(a:z:gap) * ( x_eta(a:z:gap) * q_z(a:z:gap)) / sqrt( z_eta(a:z:gap)**2 + x_eta(a:z:gap)**2 ) - &
!                   tau_n * xi_x(a:z:gap) * (-z_eta(a:z:gap) * q_x(a:z:gap)) / sqrt( z_eta(a:z:gap)**2 + x_eta(a:z:gap)**2 ) - &
!                   tau_n * xi_x(a:z:gap) * ( x_eta(a:z:gap) * q_z(a:z:gap)) / sqrt( z_eta(a:z:gap)**2 + x_eta(a:z:gap)**2 )
     Bq(a:z:gap) = Bq(a:z:gap) - hypot( xi_z(a:z:gap), xi_x(a:z:gap) ) * &
                   tau_n * (-z_eta_n(a:z:gap) * q_x(a:z:gap) + x_eta_n(a:z:gap) * q_z(a:z:gap))

  endif

  ! Left boundary condition.
  if ( rank == 0 ) then
     if ( cond(4) == 1 ) then
        ! Apply Dirichlet boundary condition operator.
        a       = 1
        z       = n * nsubz
        !Bq(a:z) = Bq(a:z) - tau_d * (eta_x(a:z) + eta_z(a:z))**2 * q(a:z)
        Bq(a:z) = Bq(a:z) - tau_d * hypot( eta_x(a:z), eta_z(a:z) )**2 * q(a:z)
     else
        ! Apply Neumann boundary condition operator.
        a       = 1
        z       = n * nsubz
!        Bq(a:z) = Bq(a:z) - &
!                  tau_n * eta_x(a:z) * (-z_xi(a:z) * q_x(a:z)) / sqrt( z_xi(a:z)**2 + x_xi(a:z)**2 ) - &
!                  tau_n * eta_x(a:z) * ( x_xi(a:z) * q_z(a:z)) / sqrt( z_xi(a:z)**2 + x_xi(a:z)**2 ) - &
!                  tau_n * eta_z(a:z) * (-z_xi(a:z) * q_x(a:z)) / sqrt( z_xi(a:z)**2 + x_xi(a:z)**2 ) - &
!                  tau_n * eta_z(a:z) * ( x_xi(a:z) * q_z(a:z)) / sqrt( z_xi(a:z)**2 + x_xi(a:z)**2 )
        Bq(a:z) = Bq(a:z) - hypot( eta_x(a:z), eta_z(a:z) ) * &
                  tau_n * (-z_xi_n(a:z) * q_x(a:z) + x_xi_n(a:z) * q_z(a:z))

        ! Let me explain what this preceeding equation is.
        !
        !    eta_x(a:z) <-- This is the metric term in the penalty coefficent,
        !                   the analog of 2/L.
        !    eta_z(a:z) <-- This is the metric term in the penalty coefficent,
        !                   the analog of 2/L.
        !    -z_xi(a:z) <-- This is nx.
        !     x_xi(a:z) <-- This is nz.
        !   sqrt( ... ) <-- This is norm(n), the normalizing factor to make
        !                   the normal vector have unit length.
        !
        ! nx = -z_xi / sqrt( z_xi**2 + x_xi**2 )
        ! nz =  x_xi / sqrt( z_xi**2 + x_xi**2 )
        ! n-dot-grad(q) = nx * q_x + nz * q_z
        !               = nx * ( eta_x * q_eta + xi_x * q_xi ) + nz * ( eta_z * q_eta + xi_z * q_xi )
        !               = ( nx * eta_x + nz * eta_z ) * q_eta  + ( nx * xi_x + nz * xi_z  ) * q_xi
     endif
  endif

end subroutine apply_bc
