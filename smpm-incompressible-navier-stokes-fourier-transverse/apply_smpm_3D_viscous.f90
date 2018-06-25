subroutine apply_smpm_3D_viscous( Aq, q )
! Applies the Spectral Multidomain Penalty method operator for the viscous
! vector equation in x/z to a vector q and returns the matrix-vector product
! Aq.

  use constants, only:         nsuby, nu
  use options, only:           bc_flag_viscous, facrobin
  use precision, only:         dp
  use timestepping, only:      dt, g0
  use woodbury_matrices, only: rpk

  implicit none

  real(kind=dp), dimension(1:3*rpk, 1:nsuby), intent(out) :: Aq
  real(kind=dp), dimension(1:3*rpk, 1:nsuby), intent(in)  :: q
  real(kind=dp)                                           :: delta
  real(kind=dp), dimension(1:3*rpk, 1:nsuby)              :: ux_grad, uy_grad, uz_grad
  real(kind=dp), dimension(1:3*rpk, 1:nsuby)              :: q_x, q_y, q_z
  integer                                                 :: x1, xN
  integer                                                 :: y1, yN
  integer                                                 :: z1, zN

  ! Get the ranges of the x-, y-, and z-components of these packed arrays.
  x1 = 1
  xN = rpk
  y1 = rpk + 1
  yN = 2 * rpk
  z1 = 2 * rpk + 1
  zN = 3 * rpk

  ! Set the viscosity to the value specified by the input file.
  delta = nu * dt / g0 ! gamma0 because of KIO time-splitting.

  ! Apply the Laplacian to each component of velocity.
  call compute_3D_gradient_and_laplacian( q(x1:xN, :), Aq(x1:xN, :), ux_grad(x1:xN, :), ux_grad(y1:yN, :), ux_grad(z1:zN, :) )
  call compute_3D_gradient_and_laplacian( q(y1:yN, :), Aq(y1:yN, :), uy_grad(y1:yN, :), uy_grad(y1:yN, :), uy_grad(z1:zN, :) )
  call compute_3D_gradient_and_laplacian( q(z1:zN, :), Aq(z1:zN, :), uz_grad(x1:xN, :), uz_grad(y1:yN, :), uz_grad(z1:zN, :) )

  ! Note:
  !
  !   ux_grad := gradient of the x-velocity.
  !   uy_grad := gradient of the y-velocity.
  !   uz_grad := gradient of the z-velocity.

  ! Subtract the identity from the pure Laplacian.
  Aq = delta * Aq - q

  ! Apply the patching condition to each component of velocity.
  call apply_3D_patching( q(x1:xN, :), Aq(x1:xN, :), delta, facrobin, ux_grad(x1:xN, :), ux_grad(z1:zN, :) )
  call apply_3D_patching( q(y1:yN, :), Aq(y1:yN, :), delta, facrobin, uy_grad(y1:yN, :), uy_grad(z1:zN, :) )
  call apply_3D_patching( q(z1:zN, :), Aq(z1:zN, :), delta, facrobin, uz_grad(z1:zN, :), uz_grad(z1:zN, :) )

  ! Reassemble components into a state vector.
  q_x(x1:xN, :) = ux_grad(x1:xN, :)
  q_x(y1:yN, :) = uy_grad(x1:xN, :)
  q_x(z1:zN, :) = uz_grad(x1:xN, :)

  q_y(x1:xN, :) = ux_grad(y1:yN, :)
  q_y(y1:yN, :) = uy_grad(y1:yN, :)
  q_y(z1:zN, :) = uz_grad(y1:yN, :)

  q_z(x1:xN, :) = ux_grad(z1:zN, :)
  q_z(y1:yN, :) = uy_grad(z1:zN, :)
  q_z(z1:zN, :) = uz_grad(z1:zN, :)

  ! NOTE:
  !
  !   q_x := x derivative of [ux, uy, uz] = [ux_x, uy_x, uz_x]
  !   q_y := y derivative of [ux, uy, uz] = [ux_y, uy_y, uz_y]
  !   q_z := z derivative of [ux, uy, uz] = [ux_z, uy_z, uz_z]

  ! Compute the 3D vector viscous boundary conditions.
  call apply_3D_vector_viscous_bc( q, Aq, q_x, q_z, bc_flag_viscous, delta )

end subroutine apply_smpm_3D_viscous
