subroutine apply_smpm_viscous( Aq, q )
! Applies the Spectral Multidomain Penalty method operator for the viscous
! vector equation in x/z to a vector q and returns the matrix-vector product
! Aq.

  use constants, only:         nu
  use options, only:           bc_flag_viscous, facrobin
  use precision, only:         dp
  use timestepping, only:      dt, g0
  use woodbury_matrices, only: rpk

  implicit none

  real(kind=dp), dimension(1:2*rpk), intent(out) :: Aq
  real(kind=dp), dimension(1:2*rpk), intent(in)  :: q
  real(kind=dp)                                  :: delta
  real(kind=dp), dimension(1:2*rpk)              :: qx, qz
  real(kind=dp), dimension(1:2*rpk)              :: q_x, q_z
  integer                                        :: x1, xN, z1, zN

  ! Get the ranges of the x- and z-components of these packed arrays.
  x1 = 1
  xN = rpk
  z1 = rpk + 1
  zN = 2 * rpk

  ! Set the viscosity to the value specified by the input file.
  delta = nu * dt / g0 ! gamma0 because of KIO time-splitting.

  ! Apply the Laplacian to each component of velocity.
  call compute_gradient_and_laplacian( q(x1:xN), Aq(x1:xN), qx(x1:xN), qx(z1:zN) )
  call compute_gradient_and_laplacian( q(z1:zN), Aq(z1:zN), qz(x1:xN), qz(z1:zN) )

  ! Subtract the identity from the pure Laplacian.
  Aq = delta * Aq - q

  ! Apply the patching condition to each component of velocity.
  call apply_patching( q(x1:xN), Aq(x1:xN), delta, facrobin, qx(x1:xN), qx(z1:zN) )
  call apply_patching( q(z1:zN), Aq(z1:zN), delta, facrobin, qz(x1:xN), qz(z1:zN) )

  ! Reassemble components into a state vector.
  q_x(      1:    rpk) = qx(x1:xN)
  q_x(rpk + 1:2 * rpk) = qz(x1:xN)
  q_z(      1:    rpk) = qx(z1:zN)
  q_z(rpk + 1:2 * rpk) = qz(z1:zN)

  ! Compute the boundary condition application (this is where the ux and uz
  ! velocities interact).
  call apply_vector_viscous_bc( q, Aq, q_x, q_z, bc_flag_viscous, delta )

end subroutine apply_smpm_viscous
