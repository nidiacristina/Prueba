subroutine apply_smpm_poisson( Ax, x )
! Applies the Poisson-Neumann operator on x.

  use options, only:           facrobin_PPE
  use precision, only:         dp
  use woodbury_matrices, only: rpk

  implicit none

  real(kind=dp), dimension(1:rpk), intent(out) :: Ax
  real(kind=dp), dimension(1:rpk), intent(in)  :: x
  real(kind=dp), dimension(1:rpk)              :: q_x, q_z

  integer, dimension(1:4), parameter           :: cond_lhsgmres = 2
  real(kind=dp), parameter                     :: delta_poisson = 1.0_dp

  ! Apply the Laplacian.
  call compute_gradient_and_laplacian( x, Ax, q_x, q_z )

  ! Apply the boundary conditions.
  call apply_bc( x, Ax, cond_lhsgmres, delta_poisson, q_x, q_z )

  ! Apply the patching conditions.
  call apply_patching( x, Ax, delta_poisson, facrobin_PPE, q_x, q_z )

end subroutine apply_smpm_poisson
