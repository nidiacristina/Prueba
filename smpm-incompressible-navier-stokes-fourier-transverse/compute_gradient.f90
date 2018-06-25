subroutine compute_gradient( v, vx, vz )
! Computes the gradient of the scalar field v and returns components (vx,vz).
!
! This computes the gradient on the local block and returns the gradient on
! just that piece.

  use constants, only:             n, nsubz
  use legendre, only:              D
  use mesh_deformation_maps, only: detJ, x_eta, x_xi, z_eta, z_xi
  use precision, only:             dp
  use woodbury_matrices, only:     dimblock, numA_per_rank, rpk

  implicit none

  real(kind=dp), dimension(1:rpk), intent(in)  :: v
  real(kind=dp), dimension(1:rpk), intent(out) :: vx
  real(kind=dp), dimension(1:rpk), intent(out) :: vz

  ! Internal variables.
  integer                                      :: iistart, iiend
  real(kind=dp), dimension(1:dimblock)         :: iiv, iiv_eta, iiv_xi

  ! Compute the divergence block-by-block.

  ! Get the start and stop indices for this block.
  iistart = 1
  iiend   = dimblock

  ! Get the data for this block.
  iiv     = v(iistart:iiend)

  ! Calculate the local derivative in xi.
  call DGEMM( 'n', 'n', n, dimblock / n, n, 1.0_dp, D, n, iiv, n, 0.0_dp, iiv_xi, n )

  ! Permute into eta-first indexing.
  call perfect_shuffle( n * numA_per_rank, n * nsubz, iiv )

  ! Calculate the local derivatives in eta.
  call DGEMM( 'n', 'n', n, dimblock / n, n, 1.0_dp, D, n, iiv, n, 0.0_dp, iiv_eta, n )

  ! Permate back into xi-first indexing.
  call perfect_shuffle( n * nsubz, n * numA_per_rank, iiv_eta )

  ! Compute the local gradient.
  vx = ( z_xi(iistart:iiend) * iiv_eta - z_eta(iistart:iiend) * iiv_xi)  / detJ(iistart:iiend)
  vz = (x_eta(iistart:iiend) * iiv_xi  -  x_xi(iistart:iiend) * iiv_eta) / detJ(iistart:iiend)

end subroutine compute_gradient
