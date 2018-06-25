subroutine apply_3D_bc( q, Bq, cond, delta, q_x, q_z )
! Computes the misfit of the penalized boundary conditions in q and adds the
! misfit to Bq.
!
! cond      - flag indicating Dirichlet or Neumann boundaries (see below).
! (q_x,q_z) - the gradient of q.
! delta     - the viscosity.

  use constants, only:             nsuby
  use precision, only:             dp
  use woodbury_matrices, only:     rpk

  implicit none

  real(kind=dp), dimension(rpk, 1:nsuby), intent(in)    :: q
  real(kind=dp), dimension(rpk, 1:nsuby), intent(inout) :: Bq

  ! Boundary condition for each edge of the quad.  edges start at the bottom
  ! and move counter-clockwise (cond(1) = bottom, cond(2) = right, cond(3) = top,
  ! and cond(4) = left.  Dirichlet condition is equal to 1, while Neumann is 2.
  integer, dimension(4), intent(in)                     :: cond
  real(kind=dp), intent(in)                             :: delta
  real(kind=dp), dimension(rpk, 1:nsuby), intent(in)    :: q_x
  real(kind=dp), dimension(rpk, 1:nsuby), intent(in)    :: q_z

  integer                                               :: ii

  ! Loop over transverse planes, applying boundary conditions on each.
  do ii = 1, nsuby
     call apply_bc( q(:, ii), Bq(:, ii), cond, delta, q_x(:, ii), q_z(:, ii) )
  enddo

end subroutine apply_3D_bc
