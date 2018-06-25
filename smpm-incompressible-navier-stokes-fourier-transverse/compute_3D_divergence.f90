subroutine compute_3D_divergence( vx, vy, vz, divv )
! Computes the divergence of the real vector field (vx,vy,vz).

  use constants, only:             nsuby
  use precision, only:             dp
  use woodbury_matrices, only:     rpk

  implicit none

  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(in)  :: vx
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(in)  :: vy
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(in)  :: vz
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out) :: divv

  real(kind=dp), dimension(1:rpk, 1:nsuby)              :: Dvy

  integer                                               :: ii

  ! Get the 2D divergence D_x * vx + D_z * vz.
  do ii = 1, nsuby
     call compute_divergence( vx(:, ii), vz(:, ii), divv(:, ii) )
  enddo

  ! Get the derivative in y. Check first if we have a transverse domain.
  if (nsuby > 1) then
     call compute_y_derivative( vy, Dvy )
  else
     Dvy = 0.0_dp
  endif
  divv = divv + Dvy

end subroutine compute_3D_divergence
