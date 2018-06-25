subroutine compute_3D_gradient_and_laplacian( u, Lu, u_x, u_y, u_z )
! Computes the Laplacian Lu of the scalar field u.  In computing the
! Laplacian, the gradient operator is implicitly computed, so this subroutine
! also returns the gradient in components (u_x,u_y,u_z).

  use constants, only:             nky, nsuby
  use precision, only:             dp
  use transverse, only:            qy
  use woodbury_matrices, only:     rpk

  implicit none

  ! Input variables.
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(in)  :: u
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out) :: Lu
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out) :: u_x
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out) :: u_y
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out) :: u_z

  ! Internal variables.
  integer                                               :: ii
  complex(kind=dp)                                      :: i
  complex(kind=dp), dimension(1:rpk, 1:nky)             :: Fu
  real(kind=dp), dimension(1:rpk, 1:nsuby)              :: u_yy

  ! Compute the 2D gradient and Laplacian plane-by-plane.
  do ii = 1, nsuby
     call compute_gradient_and_laplacian( u(:, ii), Lu(:, ii), u_x(:, ii), u_z(:, ii) )
  enddo

  ! Compute the transverse gradient
  ! Check if we have a tranvserse domain first.
  if ( nsuby > 1 ) then

     ! Set the complex unit.
     i = cmplx( 0.0_dp, 1.0_dp, kind=dp )

     ! Compute the transverse first and second derivatives.
     call apply_fft( u, Fu )
     do ii = 1, nky
        Fu(:, ii) = i * qy(ii) * Fu(:, ii)
     enddo
     call apply_ifft( Fu, u_y )
     do ii = 1, nky
        Fu(:, ii) = i * qy(ii) * Fu(:, ii)
     enddo
     call apply_ifft( Fu, u_yy )

  else

     ! Set transverse array to 0
     u_y  = 0.0_dp
     u_yy = 0.0_dp

  endif

  ! Add the second derivative to the Laplacian.
  Lu = Lu + u_yy

end subroutine compute_3D_gradient_and_laplacian
