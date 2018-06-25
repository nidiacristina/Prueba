subroutine compute_double_curl( vx, vy, vz, CCvx, CCvy, CCvz )
! Computes the double curl:
!
!     curl( curl( v ) ) = gradient( divergence( v ) ) - laplacian( v )
!
! Where v is a vector field with components (vx,vz).

  use constants,         only: nsuby
  use precision, only:         dp
  use woodbury_matrices, only: rpk

  implicit none

  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(in)  :: vx
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(in)  :: vy
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(in)  :: vz
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out) :: CCvx
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out) :: CCvy
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out) :: CCvz

  ! Internal variables.
  real(kind=dp), dimension(1:rpk, 1:nsuby)              :: Lvx, Lvy, Lvz
  real(kind=dp), dimension(1:rpk, 1:nsuby)              :: divv, Gdivvx, Gdivvy, Gdivvz

  ! Compute the vector Laplacian.  In cartesian coordinates it's component-wise
  ! scalar Laplacian.
  call compute_3D_laplacian( vx, Lvx )
  call compute_3D_laplacian( vy, Lvy )
  call compute_3D_laplacian( vz, Lvz )

  ! Compute the divergence of the vector field.
  call compute_3D_divergence( vx, vy, vz, divv )

  ! Compute the gradient of the divergence of the vector field.
  call compute_3D_gradient( divv, Gdivvx, Gdivvy, Gdivvz )

  ! Compute the double curl.
  CCvx = Gdivvx - Lvx
  CCvy = Gdivvy - Lvy
  CCvz = Gdivvz - Lvz

end subroutine compute_double_curl
