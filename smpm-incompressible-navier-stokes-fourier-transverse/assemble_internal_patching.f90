subroutine assemble_internal_patching( q, Bq, delta, fac, q_x, q_z, sub_block_number )
! Computes the misfit of the penalized patching conditions in q and adds the
! misfit to Bq, for assembly of the poisson matrix.  This function only
! computes the internal conditions that contribute to the block-diagonal
! component A of the poisson matrix.
!
! (q_x,q_z)        - the gradient of q.
! delta            - the viscosity.
! fac              - the multiplicative factor on the penalty coefficient.
! sub_block_number - which sub_block this is within the A block for this
!                    rank.

  use constants, only:             n, nprocs, nsubz, pd, rank
  use mesh_deformation_maps, only: eta_x, eta_z, x_eta_n, xi_x, xi_z, x_xi_n, z_xi_n, z_eta_n
  use precision, only:             dp
  use woodbury_matrices, only:     dimA, numA_per_rank

  implicit none

  real(kind=dp), dimension(dimA), intent(in)    :: q
  real(kind=dp), dimension(dimA), intent(inout) :: Bq
  real(kind=dp), intent(in)                     :: delta
  real(kind=dp), intent(in)                     :: fac
  real(kind=dp), dimension(dimA), intent(in)    :: q_x
  real(kind=dp), dimension(dimA), intent(in)    :: q_z
  integer, intent(in)                           :: sub_block_number

  integer                                       :: a, z, ma, mz, iistart, iiend
  integer                                       :: aB, zB, gB ! start, stop, and gap index.
  integer                                       :: aT, zT, gT ! start, stop, and gap index.
  integer                                       :: maB, mzB   ! start and stop index.
  integer                                       :: maT, mzT   ! start and stop index.
  real(kind=dp)                                 :: tau, omega, kappa, odb
  integer                                       :: ii

  real(kind=dp), parameter                      :: alpha = 1.0_dp
  real(kind=dp), parameter                      :: beta = 1.0_dp

  ! Get the ownership range of this sub-block.
  iistart = (sub_block_number - 1) * dimA + 1
  iiend   = (sub_block_number - 1) * dimA + dimA

  ! Set some constants.
  omega = 2.0_dp / (pd * (pd + 1.0_dp))
  kappa = omega * alpha / beta
  odb   = omega * delta * beta
  tau   = (1.0_dp / odb) * (delta + (2.0_dp * kappa) - (2.0_dp * sqrt( kappa**2 + (delta * kappa) )))

  ! Add the contributions of the vertical interfaces to the left and right of
  ! this sub-block.

  ! Left side of this block (right side of that interface).
  if ( (rank > 0) .or. (sub_block_number > 1)  ) then

     ! Get the indices of the left column of points in this block.
     a = 1
     z = n * nsubz

     ! Get the indices of the metric terms within the global grid.
     ma = (sub_block_number - 1) * dimA + a
     mz = (sub_block_number - 1) * dimA + z

     ! Impose the C0 continuity condition.
     Bq(a:z) = Bq(a:z) - fac * tau * hypot( eta_x(ma:mz), eta_z(ma:mz) ) * q(a:z)

     ! Impose the C1 continuity condition.
     Bq(a:z) = Bq(a:z) - fac * delta * tau * hypot( eta_x(ma:mz), eta_z(ma:mz) ) * &
               ((-z_xi_n(ma:mz) * q_x(a:z) + x_xi_n(ma:mz) * q_z(a:z)))
  endif

  ! Right side of this block (left side of that interface).
  if ( (rank < nprocs - 1) .or. (sub_block_number < numA_per_rank )) then

     ! Get the indices of the rght column of points in this block.
     a = n * n * nsubz - (n * nsubz - 1)
     z = n * n * nsubz

     ! Get the indices of the metric terms within the global grid.
     ma = (sub_block_number - 1) * dimA + a
     mz = (sub_block_number - 1) * dimA + z

     ! Impose the C0 continuity condition.
     Bq(a:z) = Bq(a:z) - fac * tau * hypot( eta_x(ma:mz), eta_z(ma:mz) ) * q(a:z)

     ! Impose the C1 continuity condition.
     Bq(a:z) = Bq(a:z) - fac * delta * tau * hypot( eta_x(ma:mz), eta_z(ma:mz) ) * &
               ((z_xi_n(ma:mz) * q_x(a:z) - x_xi_n(ma:mz) * q_z(a:z)) )
  endif

  ! Loop over horizontal interfaces.
  do ii = 1, nsubz - 1

     ! Get the indices of horizontal line of grid points on either side of
     ! this interface.
     aB   = n * ii
     zB   = n * ii + nsubz * n * (1 * n - 1) ! mz * n * (mx * n - 1) gets you
                                             ! to the last column of the grid.
     gB   = n * nsubz
     aT   = aB + 1
     zT   = zB + 1
     gT   = n * nsubz

     ! Get the indices of the metric terms on this block.
     maB = (sub_block_number - 1) * dimA + aB
     mzB = (sub_block_number - 1) * dimA + zB
     maT = (sub_block_number - 1) * dimA + aT
     mzT = (sub_block_number - 1) * dimA + zT

     ! NOTE: "B" means bottom side of the interface, "T" means top side of the
     !       interface.

     ! Impose the C0 continuity condition on each side of the interface.
     Bq(aB:zB:gB) = Bq(aB:zB:gB) - fac * tau * hypot( xi_x(maB:mzB:gB), xi_z(maB:mzB:gB) ) * (q(aB:zB:gB) - q(aT:zT:gT))
     Bq(aT:zT:gT) = Bq(aT:zT:gT) - fac * tau * hypot( xi_x(maT:mzT:gT), xi_z(maT:mzT:gT) ) * (q(aT:zT:gT) - q(aB:zB:gB))

     ! Impose the C1 continuity condition on each side of the interface.
     Bq(aB:zB:gB) = Bq(aB:zB:gB) - fac * delta * tau * hypot( xi_x(maB:mzB:gB), xi_z(maB:mzB:gB) ) * &
                    ((-z_eta_n(maB:mzB:gB) * q_x(aB:zB:gB) + x_eta_n(maB:mzB:gB) * q_z(aB:zB:gB)) - &
                     (-z_eta_n(maB:mzB:gB) * q_x(aT:zT:gT) + x_eta_n(maB:mzB:gB) * q_z(aT:zT:gT)))

     Bq(aT:zT:gT) = Bq(aT:zT:gT) - fac * delta * tau * hypot( xi_x(maT:mzT:gT), xi_z(maT:mzT:gT) ) * &
                    ((z_eta_n(maT:mzT:gT) * q_x(aT:zT:gT) - x_eta_n(maT:mzT:gT) * q_z(aT:zT:gT)) - &
                     (z_eta_n(maT:mzT:gT) * q_x(aB:zB:gB) - x_eta_n(maT:mzT:gT) * q_z(aB:zB:gB)))

  enddo

end subroutine assemble_internal_patching
