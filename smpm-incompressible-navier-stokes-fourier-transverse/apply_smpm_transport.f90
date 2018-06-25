subroutine apply_smpm_transport( Aq, q, ux_in, uy_in, uz_in )
! Computes the advective operator on q in skew-symmetric form, without the
! divergence term.

  use constants, only:             n, nprocs, nsuby, nsubz, pd, rank
  use field_variables, only:       ubc
  use mesh_deformation_maps, only: eta_x, eta_z, x_eta_n, x_xi_n, xi_x, xi_z, z_eta_n, z_xi_n
  use precision, only:             dp
  use woodbury_matrices, only:     dimA, numA_per_rank, rpk, s

  implicit none

  ! Input variables.
  real(kind=dp), intent(out), dimension(1:rpk, 1:nsuby) :: Aq
  real(kind=dp), intent(in),  dimension(1:rpk, 1:nsuby) :: q
  real(kind=dp), intent(in),  dimension(1:rpk, 1:nsuby) :: ux_in
  real(kind=dp), intent(in),  dimension(1:rpk, 1:nsuby) :: uy_in
  real(kind=dp), intent(in),  dimension(1:rpk, 1:nsuby) :: uz_in

  ! Intermediate variables for vector calculus.
  real(kind=dp), dimension(1:rpk, 1:nsuby)              :: q_x, q_y, q_z
  real(kind=dp), dimension(1:s  , 1:nsuby)              :: qleft, qright
  real(kind=dp), dimension(1:rpk)                       :: uDOTn

  ! Internal variables.
  integer                                               :: ii, jj, kk, ll, ndx
  integer                                               :: top, bottom, left, right
  real(kind=dp)                                         :: tau, omega

  ! Set some constants for the penalty parameters.
  omega = 2.0_dp / (pd * (pd + 1))
  tau   = 1.0_dp / omega

  ! Do some calculus to get derivatives we'll need.
  call compute_3D_gradient( q, q_x, q_y, q_z )

  ! Compute the advection operator in conservative form.
  Aq = -1.0_dp * (ux_in * q_x + uy_in * q_y + uz_in * q_z)

  ! Loop over the transverse planes, computing the inter-subdomain continuity conditions.
  do ll = 1, nsuby

     ! Apply the patching conditions in the vertical.

     ! Compute u-dot-n on the horizontal interfaces (normal vectors are always
     ! upward pointing).
     uDOTn = -z_eta_n * (ux_in(:, ll) + ubc(:, ll)) + x_eta_n * uz_in(:, ll)

     ! C0 continuity over n*numA_per_rank copies of a 1D multidomain grid with
     ! nsubz subdomains.
     do jj = 1, n * numA_per_rank
        do kk = 1, (nsubz - 1)

           ! Get the correct indices for the patching condition.
           bottom = (jj - 1) * (nsubz * n) + (kk * n)
           top    = (jj - 1) * (nsubz * n) + (kk * n) + 1

           ! Fix the fact that my normal vectors are strictly upward-pointing.
           uDOTn(top) = -1 * uDOTn(top)

           ! Check for inflow below the interface.
           if ( uDOTn(bottom) < 0.0_dp ) then

              ! Apply patching condition to the bottom subdomain.
              Aq(bottom, ll) = Aq(bottom, ll) - &
                               tau * hypot( xi_x(bottom), xi_z(bottom) ) * &
                               abs( uDOTn(bottom) ) * (q(bottom, ll) - q(top, ll))                               
           else

              ! Apply patching condition to the top subdomain.
              Aq(top, ll) = Aq(top, ll) - &
                            tau * hypot( xi_x(top), xi_z(top) ) * &
                            abs( uDOTn(top) ) * (q(top, ll) - q(bottom, ll))
           endif
        enddo
     enddo

     ! Compute u-dot-n on the vertical (normal vectors are always left
     ! pointing).
     uDOTn = -z_xi_n * (ux_in(:, ll) + ubc(:, ll)) + x_xi_n * uz_in(:, ll) ! Always left-pointing.

     ! C0 continuity of the internal vertical interfaces.
     do ii = 1, numA_per_rank - 1
        do kk = 1, nsubz * n

           ! Get the correct indices for the patching condition.
           left  = ii * dimA - n * nsubz + kk
           right = ii * dimA + kk

           ! Correct the right's u_dot_n so that the vector points outward.
           uDOTn(left) = - 1 * uDOTn(left) ! The left side of an interface is the right side of a subdomain.

           ! Check for inflow at the left side of the interface.
           if ( uDOTn(left) < 0 ) then
              Aq(left, ll) = Aq(left, ll) - &
                             tau * hypot( eta_x(left), eta_z(left) ) * abs( uDOTn(left) ) * (q(left, ll) - q(right, ll))
           else
              Aq(right,ll ) = Aq(right, ll) - &
                              tau * hypot( eta_x(right), eta_z(right) ) * abs( uDOTn(right) ) * (q(right, ll) - q(left, ll))
           endif

           uDOTn(left) = -1 * uDOTn(left)

        enddo
     enddo
  enddo

  ! Pack up fluxes to send to the left and the right (need to send/recv uDOTn
  ! and q).
  qleft(1:s, :)  = q(1:s, :)
  qright(1:s, :) = q(rpk - s + 1:rpk, :)

  ! Pass and receive fluxes from neighboring ranks.
  call sync_ranks( qleft, qright, s * nsuby )

  ! Apply the patching conditions on the inter-rank blocks, one transverse
  ! plane at a time.
  do ll = 1, nsuby

     ! Compute u-dot-n on the vertical (normal vectors are always left
     ! pointing).
     uDOTn = -z_xi_n * (ux_in(:, ll) + ubc(:, ll)) + x_xi_n * uz_in(:, ll) ! Always left-pointing.

     ! If we're not the first rank, apply the inter-rank patching condition on
     ! the left.
     if ( rank .ne. 0 ) then

        ! Loop over the boundary points on the left side of this rank.
        do kk = 1, s
           ! Check for inflow on this rank's side of the interface.
           if ( uDOTn(kk) < 0 ) then
              Aq(kk, ll) = Aq(kk, ll) - tau * hypot( eta_x(kk), eta_z(kk) ) * abs( uDOTn(kk) ) * (q(kk, ll) - qleft(kk, ll))
           endif
        enddo

     endif

     ! If we're not the last rank, apply the inter-rank patching condition on the
     ! right.
     if ( rank .ne. nprocs - 1 ) then

        ! Loop over the boundary points on the right side of this rank.
        do kk = 1, s
           ! Get the index within this rank of this grid point.
           ndx = rpk - s + kk

           ! Since normal vectors always point "left", flip them so that we get a
           ! negative sign for inflow.
           uDOTn(ndx) = - 1.0_dp * uDOTn(ndx)

           ! Check for inflow on this rank's side of the interface.
           if ( uDOTn(ndx) < 0 ) then
              Aq(ndx, ll) = Aq(ndx, ll) - tau * hypot( eta_x(ndx), eta_z(ndx) ) * abs( uDOTn(ndx) ) * (q(ndx, ll) - qright(kk, ll))
           endif

           ! Flip the normal vector back.
           uDOTn(ndx) = -1.0_dp * uDOTn(ndx)

        enddo
     endif
  enddo

end subroutine apply_smpm_transport
