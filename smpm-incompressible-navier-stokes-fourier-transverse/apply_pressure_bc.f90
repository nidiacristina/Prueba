subroutine apply_pressure_bc( p_rhs )
! Computes the misfit of the penalized boundary conditions on the
! right-hand-side of the pressure Poisson equation add adds it to the
! right-hand-side.

  use constants, only:             n, nprocs, nsuby, nsubz, nu, pd, rank
  use field_variables, only:       Cux0, Cux1, Cux2, Cuz0, Cuz1, Cuz2, &
                                   Nux0, Nux1, Nux2, Nuz0, Nuz1, Nuz2
  use mesh_deformation_maps, only: eta_x, eta_z, nx, nz, xi_x, xi_z
  use precision, only:             dp
  use timestepping, only:          b0, b1, b2
  use woodbury_matrices, only:     rpk

  implicit none

  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(inout) :: p_rhs
  real(kind=dp)                                           :: tau, omega
  integer                                                 :: a, z, gap, ii

  ! Set some constants.
  omega = 2.0_dp / (pd * (pd + 1.0_dp))
  tau   = 1.0_dp / omega

  ! Apply the boundary conditions one wavenumber at a time.
  do ii = 1, nsuby

     ! Left boundary.
     if ( rank == 0 ) then
        a              = 1
        z              = n * nsubz
        p_rhs(a:z, ii) = p_rhs(a:z, ii) + (tau * hypot( eta_z(a:z), eta_x(a:z) ) * &
                              (nx(a:z, 4) * ((b0 * Nux0(a:z, ii) + b1 * Nux1(a:z, ii) + b2 * Nux2(a:z, ii)) - &
                                        nu * (b0 * Cux0(a:z, ii) + b1 * Cux1(a:z, ii) + b2 * Cux2(a:z, ii))) + &
                               nz(a:z, 4) * ((b0 * Nuz0(a:z, ii) + b1 * Nuz1(a:z, ii) + b2 * Nuz2(a:z, ii)) - &
                                        nu * (b0 * Cuz0(a:z, ii) + b1 * Cuz1(a:z, ii) + b2 * Cuz2(a:z, ii)))) )
     endif

     ! Right boundary.
     if ( rank == nprocs - 1 ) then
        a              = rpk - (n * nsubz - 1)
        z              = rpk - 0
        p_rhs(a:z, ii) = p_rhs(a:z, ii) + (tau * hypot( eta_z(a:z), eta_x(a:z) ) * &
                              (nx(a:z, 2) * ((b0 * Nux0(a:z, ii) + b1 * Nux1(a:z, ii) + b2 * Nux2(a:z, ii)) - &
                                        nu * (b0 * Cux0(a:z, ii) + b1 * Cux1(a:z, ii) + b2 * Cux2(a:z, ii))) + &
                               nz(a:z, 2) * ((b0 * Nuz0(a:z, ii) + b1 * Nuz1(a:z, ii) + b2 * Nuz2(a:z, ii)) - &
                                        nu * (b0 * Cuz0(a:z, ii) + b1 * Cuz1(a:z, ii) + b2 * Cuz2(a:z, ii)))) )
     endif

     ! Top boundary.
     a                  = n * nsubz
     z                  = rpk
     gap                = n * nsubz
     p_rhs(a:z:gap, ii) = p_rhs(a:z:gap, ii) + (tau * hypot( xi_z(a:z:gap), xi_x(a:z:gap) ) * &
                               (nx(a:z:gap, 3) * ((b0 * Nux0(a:z:gap, ii) + b1 * Nux1(a:z:gap, ii) + b2 * Nux2(a:z:gap, ii)) - &
                                             nu * (b0 * Cux0(a:z:gap, ii) + b1 * Cux1(a:z:gap, ii) + b2 * Cux2(a:z:gap, ii))) + &
                                nz(a:z:gap, 3) * ((b0 * Nuz0(a:z:gap, ii) + b1 * Nuz1(a:z:gap, ii) + b2 * Nuz2(a:z:gap, ii)) - &
                                             nu * (b0 * Cuz0(a:z:gap, ii) + b1 * Cuz1(a:z:gap, ii) + b2 * Cuz2(a:z:gap, ii)))) )

     ! Bottom boundary.
     a                  = 1
     z                  = rpk
     gap                = n * nsubz
     p_rhs(a:z:gap, ii) = p_rhs(a:z:gap, ii) + (tau * hypot( xi_z(a:z:gap), xi_x(a:z:gap) ) * &
                               (nx(a:z:gap, 1) * ((b0 * Nux0(a:z:gap, ii) + b1 * Nux1(a:z:gap, ii) + b2 * Nux2(a:z:gap, ii)) - &
                                             nu * (b0 * Cux0(a:z:gap, ii) + b1 * Cux1(a:z:gap, ii) + b2 * Cux2(a:z:gap, ii))) + &
                                nz(a:z:gap, 1) * ((b0 * Nuz0(a:z:gap, ii) + b1 * Nuz1(a:z:gap, ii) + b2 * Nuz2(a:z:gap, ii) ) - &
                                             nu * (b0 * Cuz0(a:z:gap, ii) + b1 * Cuz1(a:z:gap, ii) + b2 * Cuz2(a:z:gap, ii)))) )
  enddo

end subroutine apply_pressure_bc
