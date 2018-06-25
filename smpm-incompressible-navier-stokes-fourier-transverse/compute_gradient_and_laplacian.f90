subroutine compute_gradient_and_laplacian( u, Lu, u_x, u_z )
! Computes the Laplacian Lu of the scalar field u.  In computing the
! Laplacian, the gradient operator is implicitly computed, so this subroutine
! also returns the gradient in components (u_x,u_z).

  use constants, only:             n, nsubz, rank
  use legendre, only:              D, D2
  use mesh_deformation_maps, only: detJ, d_eta_to_laplacian, d_etaeta_to_laplacian, &
                                   d_xi_to_laplacian, d_xixi_to_laplacian, d_xieta_to_laplacian, &
                                   x_eta, x_xi, z_eta, z_xi
  use precision, only:             dp
  use woodbury_matrices, only:     numA_per_rank, rpk

  implicit none

  ! Input variables.
  real(kind=dp), dimension(rpk), intent(in)  :: u
  real(kind=dp), dimension(rpk), intent(out) :: Lu
  real(kind=dp), dimension(rpk), intent(out) :: u_x
  real(kind=dp), dimension(rpk), intent(out) :: u_z

  ! Internal variables.
  integer                                    :: ii, iistart, iiend
  real(kind=dp), dimension(1:rpk)            :: iiu, iiLxz

  ! Block variables.
  real(kind=dp), dimension(1:rpk)            :: iiu_xi, iiu_eta, iiu_xixi, iiu_etaeta, iiu_xieta

  ! Apply the Laplacian.

  ! Get the start and stop indices.
  ii      = rank
  iistart = 1
  iiend   = rpk

  ! Extract this portion of the data.
  iiu     = u

  ! Apply the xi-terms ("z" direction).

  ! 2nd derivative in xi. (this is a folded iterated matvec).
  call DGEMM( 'N', 'N', n, rpk / n, n, 1.0_dp, D2, n, iiu, n, 0.0_dp, iiu_xixi, n )

  ! 1st derivative in xi. (this is a folded iterated matvec).
  call DGEMM( 'N', 'N', n, rpk / n, n, 1.0_dp, D, n, iiu, n, 0.0_dp, iiu_xi, n )

  ! 1st derivative in xi for the cross term.
  iiLxz = iiu_xi

  ! Apply the eta-terms ("x" direction).

  ! Shuffle into eta-first indexing.
  call perfect_shuffle( n * numA_per_rank, n * nsubz, iiu )

  ! 2nd derivative in eta. (this is a folded iterated matvec).
  call DGEMM( 'N', 'N', n, rpk / n, n, 1.0_dp, D2, n, iiu, n, 0.0_dp, iiu_etaeta, n )
  call perfect_shuffle( n * nsubz, n * numA_per_rank, iiu_etaeta )

  ! 1st derivative in eta. (this is a folded iterated matvec).
  call DGEMM( 'N', 'N', n, rpk / n, n, 1.0_dp, D, n, iiu, n, 0.0_dp, iiu_eta, n )
  call perfect_shuffle( n * nsubz, n * numA_per_rank, iiu_eta )

  ! 1st derivative in eta for the cross term. (this is a folded iterated
  ! matvec).
  call perfect_shuffle( n * numA_per_rank, n * nsubz, iiLxz )
  call DGEMM( 'N', 'N', n, rpk / n, n, 1.0_dp, D, n, iiLxz, n, 0.0_dp, iiu_xieta, n )
  call perfect_shuffle( n * nsubz, n * numA_per_rank, iiu_xieta )

  ! Construct this block's contribution in the deformed derivative.
  Lu =      iiu_eta    * d_eta_to_laplacian(iistart:iiend)
  Lu = Lu + iiu_xi     * d_xi_to_laplacian(iistart:iiend)
  Lu = Lu + iiu_etaeta * d_etaeta_to_laplacian(iistart:iiend)
  Lu = Lu + iiu_xixi   * d_xixi_to_laplacian(iistart:iiend)
  Lu = Lu + iiu_xieta  * d_xieta_to_laplacian(iistart:iiend)

  ! Compute this block's local gradient.
  u_x = (z_xi(iistart:iiend) * iiu_eta - z_eta(iistart:iiend) * iiu_xi) / detJ(iistart:iiend)
  u_z = (x_eta(iistart:iiend) * iiu_xi  -  x_xi(iistart:iiend) * iiu_eta) / detJ(iistart:iiend)

end subroutine compute_gradient_and_laplacian
