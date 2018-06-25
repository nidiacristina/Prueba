subroutine assemble_capacitance_block( iiC, block_number, iiAE, tau_metric_nbr, x_xi_nbr, z_xi_nbr )
! Assemble the capacitance block with internal block number block_number
! (i.e. the block_number-th block on this rank) and return it to the calling
! routine.
!
! iiAE is this current block's A\E.  iiC = B * iiAE.
!
! (q_x,q_z) - the gradient of q.
! delta     - the viscosity.
! fac       - the multiplicative factor on the penalty coefficient.

  use options, only:           facrobin_PPE
  use constants, only:         nsubx, pd, rank
  use precision, only:         dp
  use woodbury_matrices, only: dimA, numA_per_rank, rpk, s

  implicit none

  real(kind=dp), dimension(1:4*s, 1:2*s), intent(out) :: iiC
  integer, intent(in)                                 :: block_number   ! Within this rank.
  real(kind=dp), dimension(1:dimA, 1:2*s), intent(in) :: iiAE
  real(kind=dp), dimension(1:2*s), intent(in)         :: tau_metric_nbr !
  real(kind=dp), dimension(1:2*s), intent(in)         :: x_xi_nbr       !
  real(kind=dp), dimension(1:2*s), intent(in)         :: z_xi_nbr       ! Metric terms from neighboring blocks.
                                                                        ! First half is left, second half is right.

  real(kind=dp), dimension(1:4*s)                     :: BAEj
  real(kind=dp), dimension(1:dimA)                    :: AEj, AEj_x, AEj_z
  real(kind=dp), dimension(1:rpk)                     :: iiAEj, iiAEj_x, iiAEj_Z
  real(kind=dp)                                       :: fac, delta
  real(kind=dp), dimension(1:2*s)                     :: x_xi_norm, z_xi_norm

  real(kind=dp)                                       :: tau, omega, kappa, odb
  integer                                             :: jj, global_block_number
  integer                                             :: a, z

  real(kind=dp), parameter                            :: alpha = 1.0_dp
  real(kind=dp), parameter                            :: beta = 1.0_dp

  ! Set some constants.
  delta = 1.0_dp
  fac   = facrobin_PPE
  omega = 2.0_dp / (pd * (pd + 1.0_dp))
  kappa = omega * alpha / beta
  odb   = omega * delta * beta
  tau   = (1.0_dp / odb) * (delta + (2.0_dp * kappa) - (2.0_dp * sqrt( kappa**2 + (delta * kappa) )))

  ! Rename the incoming metric terms.
  !
  ! NOTE: This still exists as a hold over from when the metric terms were
  !       not normal vectors.
  x_xi_norm = x_xi_nbr
  z_xi_norm = z_xi_nbr

  ! Get some row ranges we'll need.

  ! Get the ownership range of rows in this capacitance block (for gradient
  ! computation).
  a = (block_number - 1) * dimA + 1
  z = (block_number - 1) * dimA + dimA

  ! Get the ownership range of rows in the capacitance matrix for this
  ! capacitance block.
  global_block_number = rank * numA_per_rank + block_number

  ! Loop over the columns of the A\E block.
  iiC = 0.0_dp
  do jj = 1, 2 * s

     ! Initialize the arrays for this block and compute their gradient.
     iiAEj(a:z) = 0.0_dp
     iiAEj(a:z) = iiAE(:, jj)
     call compute_gradient( iiAEj, iiAEj_x, iiAEj_z )
     AEj   = iiAEj(a:z)
     AEj_x = iiAEj_x(a:z)
     AEj_z = iiAEj_z(a:z)
     BAEj  = 0.0_dp

     ! Compute the fluxes we would have sent to the left block.
     if ( global_block_number > 1 ) then

        ! Apply the C0 continuity conditions.
        BAEj(1:s) = BAEj(1:s) - fac * tau * tau_metric_nbr(1:s) * (-1.0_dp * AEj(1:s))

        ! Apply the C1 continuity conditions.
        BAEj(1:s) = BAEj(1:s) - fac * tau * tau_metric_nbr(1:s) * delta * &
                    (-1.0_dp) * (z_xi_norm(1:s) * AEj_x(1:s) - x_xi_norm(1:s) * AEj_z(1:s))

     endif

     ! Compute the fluxes we would have sent to the right block.
     if ( global_block_number < nsubx ) then

        ! Apply the C0 continuity conditions.
        BAEj(3*s+1:4*s) = BAEj(3*s+1:4*s) - fac * tau * tau_metric_nbr(s+1:2*s) * (-1.0_dp * AEj(dimA-s+1:dimA))

        ! Apply the C1 continuity conditions.
        !
        ! NOTE: the -1 in front of z_xi_norm accounts for the normal vector
        !       of the _other_ block, not this one.
        BAEj(3*s+1:4*s) = BAEj(3*s+1:4*s) - fac * tau * tau_metric_nbr(s+1:2*s) * delta * &
                          (-1.0_dp) * (-z_xi_norm(s+1:2*s) * AEj_x(dimA-s+1:dimA) + x_xi_norm(s+1:2*s) * AEj_z(dimA-s+1:dimA))

     endif

     ! Store this column in the capacitance block.
     iiC(:, jj) = BAEj

  enddo

end subroutine assemble_capacitance_block
