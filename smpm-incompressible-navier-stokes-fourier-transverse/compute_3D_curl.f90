subroutine compute_3D_curl( vx, vy, vz, curli, curlj, curlk )
! Computes the curl of the vector field (vx,vy,vz).

  use constants, only: nsuby
  use precision, only: dp
  use woodbury_matrices, only: rpk

  implicit none

  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(in)  :: vx
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(in)  :: vy
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(in)  :: vz
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out) :: curli
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out) :: curlj
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out) :: curlk

  ! Internal variables
  real(kind=dp), dimension(1:rpk, 1:nsuby)              :: vx_x, vx_y, vx_z
  real(kind=dp), dimension(1:rpk, 1:nsuby)              :: vy_x, vy_y, vy_z
  real(kind=dp), dimension(1:rpk, 1:nsuby)              :: vz_x, vz_y, vz_z

  ! Set the curl to zero.
  curli = 0.0_dp
  curlj = 0.0_dp
  curlk = 0.0_dp

  ! Compute the gradients we'll need for the curl computation.
  call compute_3D_gradient( vx, vx_x, vx_y, vx_z )
  call compute_3D_gradient( vy, vy_x, vy_y, vy_z )
  call compute_3D_gradient( vz, vz_x, vz_y, vz_z )

  ! Assemble the curl.
  curli = vz_y - vy_z
  curlj = vx_z - vz_x ! -(vz_x - vx_z)
  curlk = vy_x - vx_y

end subroutine compute_3D_curl
