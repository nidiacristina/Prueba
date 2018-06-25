subroutine compute_divergence( vx, vz, divv )
! Computes the divergence of vector field (vx,vz).

  use constants, only:             n, nsubz
  use legendre, only:              D
  use mesh_deformation_maps, only: eta_x, eta_z, xi_x, xi_z
  use precision, only:             dp
  use woodbury_matrices, only:     dimblock, numA_per_rank, rpk

  implicit none

  real(kind=dp), dimension(1:rpk), intent(in)  :: vx
  real(kind=dp), dimension(1:rpk), intent(in)  :: vz
  real(kind=dp), dimension(1:rpk), intent(out) :: divv

  ! Internal variables.
  integer                                      :: iistart, iiend
  real(kind=dp), dimension(1:rpk)              :: iivx, iivz, iivx_xi, iivx_eta, iivz_xi, iivz_eta

  ! Get the start and stop indices for this block.
  iistart = 1
  iiend   = rpk

  ! Get the data for this block.
  iivx    = vx(iistart:iiend)
  iivz    = vz(iistart:iiend)

  ! Calculate the local derivatives in xi.
  call DGEMM( 'n', 'n', n, dimblock / n, n, 1.0_dp, D, n, iivx, n, 0.0_dp, iivx_xi, n )
  call DGEMM( 'n', 'n', n, dimblock / n, n, 1.0_dp, D, n, iivz, n, 0.0_dp, iivz_xi, n )

  ! Permute into eta-first indexing.
  call perfect_shuffle( n * numA_per_rank, n * nsubz, iivx )
  call perfect_shuffle( n * numA_per_rank, n * nsubz, iivz )

  ! Calculate the local derivatives in eta.
  call DGEMM( 'n', 'n', n, dimblock / n, n, 1.0_dp, D, n, iivx, n, 0.0_dp, iivx_eta, n )
  call DGEMM( 'n', 'n', n, dimblock / n, n, 1.0_dp, D, n, iivz, n, 0.0_dp, iivz_eta, n )

  ! Permate back into xi-first indexing.
  call perfect_shuffle( n * nsubz, n * numA_per_rank, iivx_eta )
  call perfect_shuffle( n * nsubz, n * numA_per_rank, iivz_eta )

  ! Compute the local divergence.
  divv = eta_x * iivx_eta + xi_x * iivx_xi + xi_z * iivz_xi + eta_z * iivz_eta

end subroutine compute_divergence
