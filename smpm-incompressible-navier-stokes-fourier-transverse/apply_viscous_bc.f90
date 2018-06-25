subroutine apply_viscous_bc( q, q_b, q_flag, delta )
! Takes the boundary condition value q_b and adds the misfit of q and q_b to
! q.
!
! q_flag - array of length(4) indicating Dirichlet or Neumann (see below).
! delta  - viscosity * dt / g0.

  use constants, only:             n, nprocs, nsuby, nsubz, pd, rank
  use mesh_deformation_maps, only: eta_x, eta_z, xi_x, xi_z
  use precision, only:             dp
  use woodbury_matrices, only:     rpk

  implicit none

  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(inout) :: q      ! The field variable.
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(in)    :: q_b    ! The boundary values (only those indices of
                                                                    ! this array on the boundary will be accessed).
  integer, dimension(1:4), intent(in)                     :: q_flag ! q_flag(1) is the bottom; work counterclockwise. 
                                                                    ! 1 = Dirichlet; 2 = Neumann.
  real(kind=dp), intent(in)                               :: delta
  real(kind=dp), parameter                                :: alpha = 1.0_dp
  ! constants defining the penalty parameter for the Dirichlet and Neumann
  ! boundaries.
  !
  ! NOTE: these should be parameters, though Fortran does not allow them to be
  !       since they're derived from a module variable (pd).
  real(kind=dp)                                           :: omega, tau_d, tau_n
  integer                                                 :: a, z, gap, ii

  ! Set some constants.
  omega = 2.0_dp / (pd * (pd + 1))
  tau_d = delta / alpha / omega**2
  tau_n = delta / omega

  ! Apply the boundary conditions one transverse plane at a time.
  do ii = 1, nsuby

     ! Apply the boundary condition to the left boundary.
     if ( rank == 0 ) then
        a      = 1
        z      = n * nsubz
        if ( q_flag(4) == 1 ) then
           q(a:z,ii) = q(a:z,ii) + tau_d * hypot(eta_x(a:z), eta_z(a:z))**2 * q_b(a:z,ii)
        else
           q(a:z,ii) = q(a:z,ii) + tau_n * hypot(eta_x(a:z), eta_z(a:z)) * q_b(a:z,ii)
        endif
     endif

     ! Apply the boundary condition to the right boundary.
     if ( rank == nprocs - 1 ) then
        a      = rpk - (n * nsubz - 1)
        z      = rpk - 0
        if ( q_flag(2) == 1 ) then
           q(a:z,ii) = q(a:z,ii) + tau_d * hypot(eta_x(a:z), eta_z(a:z))**2 * q_b(a:z,ii)
        else
           q(a:z,ii) = q(a:z,ii) + tau_n * hypot(eta_x(a:z), eta_z(a:z)) * q_b(a:z,ii)
        endif
     endif

     ! Apply the boundary condition to the top boundary.
     a          = n * nsubz
     z          = rpk
     gap        = n * nsubz
     if ( q_flag(3) == 1 ) then
        q(a:z:gap,ii) = q(a:z:gap,ii) + tau_d * hypot(xi_x(a:z:gap), xi_z(a:z:gap))**2 * q_b(a:z:gap,ii)
     else
        q(a:z:gap,ii) = q(a:z:gap,ii) + tau_n * hypot(xi_x(a:z:gap), xi_z(a:z:gap)) * q_b(a:z:gap,ii)
     endif

     ! Apply the boundary condition to the bottom boundary.
     a          = 1
     z          = rpk
     gap        = n * nsubz
     if ( q_flag(1) == 1 ) then
        q(a:z:gap,ii) = q(a:z:gap,ii) + tau_d * hypot(xi_z(a:z:gap), xi_x(a:z:gap))**2 * q_b(a:z:gap,ii)
     else
        q(a:z:gap,ii) = q(a:z:gap,ii) + tau_n * hypot(xi_z(a:z:gap), xi_x(a:z:gap)) * q_b(a:z:gap,ii)
     endif

  enddo

end subroutine apply_viscous_bc
