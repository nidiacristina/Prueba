subroutine compute_y_derivative( v, Dyv )
! Computes the transverse derivative of a 3D scalar field.

  use constants, only:          nky, nsuby
  use precision, only:          dp
  use transverse, only:         qy
  use woodbury_matrices, only:  rpk

  implicit none

  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(in)  :: v
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out) :: Dyv

  complex(kind=dp), dimension(1:rpk, 1:nky)             :: Fv
  complex(kind=dp)                                      :: i
  integer                                               :: ii

  ! Set the complex unit.
  i = cmplx( 0.0_dp, 1.0_dp, kind=dp )

  ! Get the Fourier transform.
  call apply_fft( v, Fv )

  ! Multiply by wavenumbers.
  do ii = 1, nky
     Fv(:, ii) = i * qy(ii) * Fv(:, ii)
  enddo

  ! Get the inverse Fourier transform.
  call apply_ifft( Fv, Dyv )

end subroutine compute_y_derivative
