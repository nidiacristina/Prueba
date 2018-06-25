subroutine apply_smpm_diffusion( Ax, x )
! Applies the diffusive operator on x.

  use constants, only:         nu_d
  use options, only:           bc_flag_diffusion, facrobin
  use precision, only:         dp
  use timestepping, only:      dt, g0
  use woodbury_matrices, only: rpk

  implicit none

  real(kind=dp), dimension(1:rpk), intent(out) :: Ax
  real(kind=dp), dimension(1:rpk), intent(in)  :: x
  real(kind=dp)                                :: delta
  real(kind=dp), dimension(1:rpk)              :: rho_x, rho_z

  ! Set the viscosity to the value specified by the input file.
  delta = nu_d * dt / g0

  ! Apply the Laplacian.
  call compute_gradient_and_laplacian( x, Ax, rho_x, rho_z )
  Ax = delta * Ax - x

  ! Apply the boundary condition.
  call apply_bc( x, Ax, bc_flag_diffusion, delta, rho_x, rho_z )

  ! Apply the patching condition.
  call apply_patching( x, Ax, delta, facrobin, rho_x, rho_z )

end subroutine apply_smpm_diffusion
