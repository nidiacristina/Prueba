subroutine compute_3D_gradient( v, vx, vy, vz )
! Computes the gradient of the scalar field v and returns components
! (vx,vy,vz).

  use constants, only:             nsuby
  use precision, only:             dp
  use woodbury_matrices, only:     rpk

  implicit none

  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(in)  :: v
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out) :: vx
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out) :: vy
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out) :: vz

  ! Internal variables.
  integer                                               :: ii

  ! Compute the derivatives in x and z.
  do ii = 1, nsuby
     call compute_gradient( v(:, ii), vx(:, ii), vz(:, ii) )
  enddo

  ! Compute the derivative in y. 
  ! Check if we have a transverse domain first.
  if (nsuby > 1 ) then
     call compute_y_derivative( v, vy )
  else
     vy = 0.0_dp
  endif

end subroutine compute_3D_gradient
