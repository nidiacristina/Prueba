subroutine apply_3D_patching( q, Bq, delta, fac, q_x, q_z )
! Computes the misfit of the penalized patching conditions in q and adds the
! misfit to Bq.
!
! cond      - flag indicating Dirichlet or Neumann boundaries (see below).
! (q_x,q_z) - the gradient of q.
! delta     - the viscosity.
! fac       - the multiplicative factor on the penalty coefficient.

  use constants, only:             n, nprocs, nsuby, nsubz, pd, rank
  use mesh_deformation_maps, only: eta_x, eta_z, x_eta_n, x_xi_n, xi_x, xi_z, z_eta_n, z_xi_n
  use precision, only:             dp
  use woodbury_matrices, only:     dimA, numA_per_rank, rpk, s

  implicit none

  real(kind=dp), dimension(rpk, 1:nsuby), intent(in)    :: q
  real(kind=dp), dimension(rpk, 1:nsuby), intent(inout) :: Bq
  real(kind=dp), intent(in)                             :: delta
  real(kind=dp), intent(in)                             :: fac
  real(kind=dp), dimension(rpk, 1:nsuby), intent(in)    :: q_x, q_z

  integer                                               :: aL, zL     ! start, and stop index.
  integer                                               :: aR, zR     ! start, and stop index.
  integer                                               :: aB, zB, gB ! start, stop, and gap index.
  integer                                               :: aT, zT, gT ! start, stop, and gap index.
  real(kind=dp)                                         :: tau, omega, kappa, odb
  integer                                               :: ii, jj, a, z
  real(kind=dp), dimension(1:3*s, 1:nsuby)              :: flux_L, flux_R

  real(kind=dp), parameter                              :: alpha = 1.0_dp
  real(kind=dp), parameter                              :: beta = 1.0_dp

  ! Set some constants.
  omega = 2.0_dp / (pd * (pd + 1.0_dp))
  kappa = omega * alpha / beta
  odb   = omega * delta * beta
  if ( odb > 0.0_dp ) then
     tau = (1.0_dp / odb) * (delta + (2.0_dp * kappa) - (2 * sqrt( kappa**2 + (delta * kappa) )))
  else
     tau = 0.0_dp
  endif

  ! Apply the patching conditions one transverse plane at a time.
  do jj = 1, nsuby

     ! Loop over internal vertical interfaces.
     do ii = 1, numA_per_rank - 1

        ! Get the indices of vertical line of grid points on either side of
        ! this interface.
        aL   = dimA * ii - (n * nsubz - 1)
        zL   = dimA * ii - 0
        aR   = zL + 1
        zR   = zL + n * nsubz

        ! NOTE: "L" means left side of the interface, "R" means right side of
        !       the interface.

        ! Impose the C0 continuity condition on each side of the interface.
        Bq(aL:zL, jj) = Bq(aL:zL, jj) - fac * tau * hypot( eta_x(aL:zL), eta_z(aL:zL) ) * &
                        (q(aL:zL, jj) - q(aR:zR, jj))
        Bq(aR:zR, jj) = Bq(aR:zR, jj) - fac * tau * hypot( eta_x(aR:zR), eta_z(aR:zR) ) * &
                        (q(aR:zR, jj) - q(aL:zL, jj))

        ! Impose the C1 continuity condition on each side of the interface.
        Bq(aL:zL, jj) = Bq(aL:zL, jj) - fac * delta * tau * hypot( eta_x(aL:zL), eta_z(aL:zL) ) * &
                       ((z_xi_n(aL:zL) * q_x(aL:zL, jj) - x_xi_n(aL:zL) * q_z(aL:zL, jj)) - &
                        (z_xi_n(aL:zL) * q_x(aR:zR, jj) - x_xi_n(aL:zL) * q_z(aR:zR, jj)))
        Bq(aR:zR, jj) = Bq(aR:zR, jj) - fac * delta * tau * hypot( eta_x(aR:zR), eta_z(aR:zR) ) * &
                       ((-z_xi_n(aR:zR) * q_x(aR:zR, jj) + x_xi_n(aR:zR) * q_z(aR:zR, jj)) - &
                        (-z_xi_n(aR:zR) * q_x(aL:zL, jj) + x_xi_n(aR:zR) * q_z(aL:zL, jj)))
     enddo

     ! Loop over horizontal interfaces.
     do ii = 1, nsubz - 1

        ! Get the indices of horizontal line of grid points on either side of
        ! this interface.
        aB = n * ii
        zB = n * ii + nsubz * n * (numA_per_rank * n - 1) ! mz * n * ( mx * n - 1 ) gets you
                                                          ! to the last column of the grid.
        gB = n * nsubz
        aT = aB + 1
        zT = zB + 1
        gT = n * nsubz

        ! NOTE: "B" means bottom side of the interface, "T" means top side of
        !       the interface.

        ! Impose the C0 continuity condition on each side of the interface.
        Bq(aB:zB:gB, jj) = Bq(aB:zB:gB, jj) - fac * tau * hypot( xi_x(aB:zB:gB), xi_z(aB:zB:gB) ) * &
                           (q(aB:zB:gB, jj) - q(aT:zT:gT, jj))
        Bq(aT:zT:gT, jj) = Bq(aT:zT:gT, jj) - fac * tau * hypot( xi_x(aT:zT:gT), xi_z(aT:zT:gT) ) * &
                           (q(aT:zT:gT, jj) - q(aB:zB:gB, jj))

        ! Impose the C1 continuity condition on each side of the interface.
        Bq(aB:zB:gB, jj) = Bq(aB:zB:gB, jj) - fac * delta * tau * hypot( xi_x(aB:zB:gB), xi_z(aB:zB:gB) ) * &
                           ((-z_eta_n(aB:zB:gB) * q_x(aB:zB:gB, jj) + x_eta_n(aB:zB:gB) * q_z(aB:zB:gB, jj)) - &
                            (-z_eta_n(aB:zB:gB) * q_x(aT:zT:gT, jj) + x_eta_n(aB:zB:gB) * q_z(aT:zT:gT, jj)))
        Bq(aT:zT:gT, jj) = Bq(aT:zT:gT, jj) - fac * delta * tau * hypot(xi_x(aT:zT:gT), xi_z(aT:zT:gT)) * &
                           ((z_eta_n(aT:zT:gT) * q_x(aT:zT:gT, jj) - x_eta_n(aT:zT:gT) * q_z(aT:zT:gT, jj)) - &
                            (z_eta_n(aT:zT:gT) * q_x(aB:zB:gB, jj) - x_eta_n(aT:zT:gT) * q_z(aB:zB:gB, jj)))
     enddo
  enddo

  ! Compute the inter-rank interfaces.

  ! Pack up the fluxes for communication.
  a                          = 1
  z                          = s
  flux_L(    1:s,   1:nsuby) =   q(a:z, :)
  flux_L(  s+1:2*s, 1:nsuby) = q_x(a:z, :)
  flux_L(2*s+1:3*s, 1:nsuby) = q_z(a:z, :)

  a                          = rpk - s + 1
  z                          = rpk
  flux_R(    1:s,   1:nsuby) =   q(a:z, :)
  flux_R(  s+1:2*s, 1:nsuby) = q_x(a:z, :)
  flux_R(2*s+1:3*s, 1:nsuby) = q_z(a:z, :)

  ! Send and receive fluxes.
  call sync_ranks( flux_L, flux_R, 3 * s * nsuby )

  ! Apply patching conditions on inter-rank interfaces.

  ! Apply to the left side of this rank.
  if ( rank > 0 ) then

     ! Apply patching conditions one transverse plane at a time.
     do jj = 1, nsuby

        ! Set the bounds for this block.
        a = 1
        z = s

        ! Apply C0 condition.
        Bq(a:z, jj) = Bq(a:z, jj) - fac * tau * hypot( eta_x(a:z), eta_z(a:z) ) * &
                      (q(a:z, jj) - flux_L(1:s, jj))

        ! Apply C1 condition.
        Bq(a:z, jj) = Bq(a:z, jj) - fac * delta * tau * hypot( eta_x(a:z), eta_z(a:z) ) * &
                      ((-z_xi_n(a:z) * q_x(a:z, jj)        + x_xi_n(a:z) * q_z(a:z, jj)) - &
                       (-z_xi_n(a:z) * flux_L(s+1:2*s, jj) + x_xi_n(a:z) * flux_L(2*s+1:3*s, jj)))
     enddo

  endif

  ! Apply to the right side of this rank.
  if ( rank < nprocs - 1 ) then

     ! Set the bounds for this block.
     a = rpk - s + 1
     z = rpk

     ! apply patching conditions one transverse plane at a time.
     do jj = 1, nsuby

        ! Apply C0 condition.
        Bq(a:z, jj) = Bq(a:z, jj) - fac * tau * hypot( eta_x(a:z), eta_z(a:z) ) * &
                      (q(a:z, jj) - flux_R(1:s, jj))

        ! Apply C1 condition.
        Bq(a:z, jj) = Bq(a:z, jj) - fac * delta * tau * hypot( eta_x(a:z), eta_z(a:z) ) * &
                      ((z_xi_n(a:z) * q_x(a:z, jj)        - x_xi_n(a:z) * q_z(a:z, jj)) - &
                       (z_xi_n(a:z) * flux_R(s+1:2*s, jj) - x_xi_n(a:z) * flux_R(2*s+1:3*s, jj)))
     enddo

  endif

end subroutine apply_3D_patching
