subroutine compute_curl( vx, vz, curlv )
! Computes the curl of the vector field (vx,vz).

  use precision, only:         dp
  use woodbury_matrices, only: rpk

  implicit none

  real(kind=dp), dimension(1:rpk), intent(in)  :: vx
  real(kind=dp), dimension(1:rpk), intent(in)  :: vz
  real(kind=dp), dimension(1:rpk), intent(out) :: curlv

  ! Internal variables
  real(kind=dp), dimension(1:rpk)              :: vz_x, vx_z, vx_x, vz_z

  ! Set the curl to zero.
  curlv = 0.0_dp

  ! Compute the gradients we'll need for the curl computation.
  call compute_gradient( vx, vx_x, vx_z )
  call compute_gradient( vz, vz_x, vz_z )

  ! Assemble the curl.
  curlv = vz_x - vx_z

end subroutine compute_curl
