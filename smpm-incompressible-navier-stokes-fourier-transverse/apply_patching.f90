subroutine apply_patching( q, Bq, delta, fac, q_x, q_z )
! Computes the misfit of the penalized patching conditions in q and adds the
! misfit to Bq.
!
! cond      - flag indicating Dirichlet or Neumann boundaries (see below).
! (q_x,q_z) - the gradient of q.
! delta     - the viscosity.
! fac       - the multiplicative factor on the penalty coefficient.

  use constants, only:             n, nprocs, nsubz, pd, rank
  use mesh_deformation_maps, only: eta_x, eta_z, x_eta_n, x_xi_n, xi_x, xi_z, z_eta_n, z_xi_n
  use precision, only:             dp
  use woodbury_matrices, only:     dimA, numA_per_rank, rpk, s

  implicit none

  real(kind=dp), dimension(rpk), intent(in)    :: q
  real(kind=dp), dimension(rpk), intent(inout) :: Bq
  real(kind=dp), intent(in)                    :: delta
  real(kind=dp), intent(in)                    :: fac
  real(kind=dp), dimension(rpk), intent(in)    :: q_x, q_z

  integer                                      :: aL, zL     ! start, and stop index.
  integer                                      :: aR, zR     ! start, and stop index.
  integer                                      :: aB, zB, gB ! start, stop, and gap index.
  integer                                      :: aT, zT, gT ! start, stop, and gap index.
  real(kind=dp)                                :: tau, omega, kappa, odb
  integer                                      :: ii, a, z
  real(kind=dp), dimension(1:3*s)              :: flux_L, flux_R

  real(kind=dp), parameter                     :: alpha = 1.0_dp
  real(kind=dp), parameter                     :: beta = 1.0_dp

  ! Set some constants.
  omega = 2.0_dp / (pd * (pd + 1.0_dp))
  kappa = omega * alpha / beta
  odb   = omega * delta * beta
  if ( odb > 0.0_dp ) then
     tau = (1.0_dp / odb) * (delta + (2.0_dp * kappa) - (2.0_dp * sqrt( kappa**2 + (delta * kappa) )))
  else
     tau = 0.0_dp
  endif

  ! Loop over internal interfaces first.

  ! Loop over vertical interfaces.
  do ii = 1, numA_per_rank - 1

     ! Get the indices of vertical line of grid points on either side of this
     ! interface.
     aL   = dimA * ii - (n * nsubz - 1)
     zL   = dimA * ii - 0
     aR   = zL + 1
     zR   = zL + n * nsubz

     ! NOTE: "L" means left side of the interface, "R" means right side of the
     !       interface.

     ! Impose the C0 continuity condition on each side of the interface.
     Bq(aL:zL) = Bq(aL:zL) - fac * tau * hypot( eta_x(aL:zL), eta_z(aL:zL) ) * &
                 (q(aL:zL) - q(aR:zR))
     Bq(aR:zR) = Bq(aR:zR) - fac * tau * hypot( eta_x(aR:zR), eta_z(aR:zR) ) * &
                 (q(aR:zR) - q(aL:zL))

     ! Impose the C1 continuity condition on each side of the interface.
     Bq(aL:zL) = Bq(aL:zL) - fac * delta * tau * hypot( eta_x(aL:zL), eta_z(aL:zL) ) * &
                 ((z_xi_n(aL:zL) * q_x(aL:zL) - x_xi_n(aL:zL) * q_z(aL:zL)) - &
                  (z_xi_n(aL:zL) * q_x(aR:zR) - x_xi_n(aL:zL) * q_z(aR:zR)))
     Bq(aR:zR) = Bq(aR:zR) - fac * delta * tau * hypot( eta_x(aR:zR), eta_z(aR:zR) ) * &
                 ((-z_xi_n(aR:zR) * q_x(aR:zR) + x_xi_n(aR:zR) * q_z(aR:zR)) - &
                  (-z_xi_n(aR:zR) * q_x(aL:zL) + x_xi_n(aR:zR) * q_z(aL:zL)))
  enddo

  ! Loop over horizontal interfaces.
  do ii = 1, nsubz - 1

     ! Get the indices of horizontal line of grid points on either side of
     ! this interface.
     aB = n * ii
     zB = n * ii + nsubz * n * (numA_per_rank * n - 1) ! NOTE: mz * n * (mx * n - 1) gets you
                                                       !       to the last column of the grid.
     gB = n * nsubz
     aT = aB + 1
     zT = zB + 1
     gT = n * nsubz

     ! NOTE: "B" means bottom side of the interface, "T" means top side of the
     !       interface.

     ! Impose the C0 continuity condition on each side of the interface.
     Bq(aB:zB:gB) = Bq(aB:zB:gB) - fac * tau * hypot( xi_x(aB:zB:gB), xi_z(aB:zB:gB) ) * &
                    (q(aB:zB:gB) - q(aT:zT:gT))
     Bq(aT:zT:gT) = Bq(aT:zT:gT) - fac * tau * hypot( xi_x(aT:zT:gT), xi_z(aT:zT:gT) ) * &
                    (q(aT:zT:gT) - q(aB:zB:gB))

     ! Impose the C1 continuity condition on each side of the interface.
     Bq(aB:zB:gB) = Bq(aB:zB:gB) - fac * delta * tau * hypot( xi_x(aB:zB:gB), xi_z(aB:zB:gB) ) * &
                    ((-z_eta_n(aB:zB:gB) * q_x(aB:zB:gB) + x_eta_n(aB:zB:gB) * q_z(aB:zB:gB)) - &
                     (-z_eta_n(aB:zB:gB) * q_x(aT:zT:gT) + x_eta_n(aB:zB:gB) * q_z(aT:zT:gT)))
     Bq(aT:zT:gT) = Bq(aT:zT:gT) - fac * delta * tau * hypot( xi_x(aT:zT:gT), xi_z(aT:zT:gT) ) * &
                    ((z_eta_n(aT:zT:gT) * q_x(aT:zT:gT) - x_eta_n(aT:zT:gT) * q_z(aT:zT:gT)) - &
                     (z_eta_n(aT:zT:gT) * q_x(aB:zB:gB) - x_eta_n(aT:zT:gT) * q_z(aB:zB:gB)))
  enddo

  ! Compute the inter-rank interfaces.

  ! Pack up the fluxes for communication.
  a                 = 1
  z                 = s
  flux_L(    1:s  ) =   q(a:z)
  flux_L(  s+1:2*s) = q_x(a:z)
  flux_L(2*s+1:3*s) = q_z(a:z)

  a                 = rpk - s + 1
  z                 = rpk
  flux_R(    1:s  ) =   q(a:z)
  flux_R(  s+1:2*s) = q_x(a:z)
  flux_R(2*s+1:3*s) = q_z(a:z)

  ! Send and receive fluxes.
  call sync_ranks( flux_L, flux_R, 3 * s )

  ! Apply patching conditions on inter-rank interfaces.

  ! Apply to the left side of this rank.
  if ( rank > 0 ) then

     ! Set the bounds for this block.
     a = 1
     z = s

     ! Apply C0 condition.
     Bq(a:z) = Bq(a:z) - fac * tau * hypot( eta_x(a:z), eta_z(a:z) ) * &
               (q(a:z) - flux_L(1:s))

     ! Apply C1 condition.
     Bq(a:z) = Bq(a:z) - fac * delta * tau * hypot( eta_x(a:z), eta_z(a:z) ) * &
               ((-z_xi_n(a:z) * q_x(a:z)        + x_xi_n(a:z) * q_z(a:z)) - &
                (-z_xi_n(a:z) * flux_L(s+1:2*s) + x_xi_n(a:z) * flux_L(2*s+1:3*s)))

  endif

  ! Apply to the right side of this rank.
  if ( rank < nprocs - 1 ) then

     ! Set the bounds for this block.
     a = rpk - s + 1
     z = rpk

     ! Apply C0 condition.
     Bq(a:z) = Bq(a:z) - fac * tau * hypot( eta_x(a:z), eta_z(a:z) ) * &
               (q(a:z) - flux_R(1:s))

     ! Apply C1 condition.
     Bq(a:z) = Bq(a:z) - fac * delta * tau * hypot( eta_x(a:z), eta_z(a:z) ) * &
               ((z_xi_n(a:z) * q_x(a:z)        - x_xi_n(a:z) * q_z(a:z)) - &
                (z_xi_n(a:z) * flux_R(s+1:2*s) - x_xi_n(a:z) * flux_R(2*s+1:3*s)))

  endif

end subroutine apply_patching
