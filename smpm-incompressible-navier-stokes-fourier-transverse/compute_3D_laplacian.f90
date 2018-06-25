subroutine compute_3D_laplacian( u, Lu )
! Computes the Laplacian Lu of the scalar field u.  In computing the
! Laplacian, the gradient operator is implicitly computed, so this subroutine
! also returns the gradient in components (u_x,u_y,u_z).

  use constants, only:             nky, nsuby
  use precision, only:             dp
  use woodbury_matrices, only:     rpk
  use transverse, only:            qy

  implicit none

  ! Input variables.
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out) :: Lu
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(in)  :: u

  ! Internal variables.
  integer                                               :: ii
  complex(kind=dp)                                      :: i
  complex(kind=dp), dimension(1:rpk, 1:nky)             :: Fu
  real(kind=dp), dimension(1:rpk, 1:nsuby)              :: u_yy

  ! Compute the 2D gradient and Laplacian plane-by-plane.
  do ii = 1, nsuby
     call compute_laplacian( u(:, ii), Lu(:, ii) )
  enddo

  ! Compute the transverse gradient.
  ! Check if we have a transverse domain first.
  if ( nsuby > 1 ) then

     ! Set the complex unit.
     i = (0.0_dp, 1.0_dp)

     ! Compute the transverse first and second derivatives.
     call apply_fft( u, Fu )
     do ii = 1, nky
        Fu(:, ii) = i * i * qy(ii) * qy(ii) * Fu(:, ii)  
     enddo
     call apply_ifft( Fu, u_yy )

  else

     u_yy = 0.0_dp

  endif

  ! Add the second derivative to the Laplacian.
  Lu = Lu + u_yy

end subroutine compute_3D_laplacian

