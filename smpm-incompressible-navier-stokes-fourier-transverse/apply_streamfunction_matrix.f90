subroutine apply_streamfunction_matrix( Aq, q )
! This function assumes xi-first indexing.  This function is identical to
! apply_smpm_poisson but with Dirichlet boundary conditions instead of
! Neumann.

  use precision, only:         dp
  use woodbury_matrices, only: rpk

  implicit none

  real(kind=dp), dimension(1:rpk), intent(out) :: Aq
  real(kind=dp), dimension(1:rpk), intent(in)  :: q
  integer, dimension(1:4)                      :: cond_lhsgmres = (/ 1, 1, 1, 1 /)
  real(kind=dp)                                :: delta_poisson = 1.0_dp
  real(kind=dp), dimension(1:rpk)              :: q_x, q_z

  ! Apply the Laplacian.
  call compute_gradient_and_laplacian( q, Aq, q_x, q_z )

  ! Apply the boundary conditions.
  call apply_bc( q, Aq, cond_lhsgmres, delta_poisson, q_x, q_z )

  ! Apply the patching conditions.
  call apply_patching( q, Aq, delta_poisson, 10000.0_dp, q_x, q_z )

end subroutine apply_streamfunction_matrix
