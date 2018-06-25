subroutine validate_fftw( tolerance, success_flag )
! Checks to make sure FFTW is working properly.

  use constants, only:         nky, nsuby
  use parallel_linear_algebra
  use precision, only:         dp
  use woodbury_matrices, only: rpk

  implicit none

  real(kind=dp), intent(in)                      :: tolerance
  logical, intent(out)                           :: success_flag

  ! Buffers for our random numbers and their Fourier spectra.
  real(kind=dp), allocatable, dimension(:, :)    :: realdata, ifftdata
  complex(kind=dp), allocatable, dimension(:, :) :: fftdata
 
  real(kind=dp)                                  :: error

  character(len=64)                              :: caststr

  allocate( realdata(1:rpk, 1:nsuby), fftdata(1:rpk, 1:nky), ifftdata(1:rpk, 1:nsuby) )

  call notify( 'Testing FFTW.')

  ! FFT a buffer of random numbers and verify that we get the same values when
  ! we do an inverse FFT.
  call random_number( realdata )
  call  apply_fft( realdata, fftdata )
  call apply_ifft( fftdata, ifftdata )

  ! Compute error
  error = pmaxval( reshape( abs( realdata - ifftdata ), (/rpk * nsuby/) ) )

  write( caststr, '(A,D17.10)' ) 'Error in FFT/IFFT pair:', error
  call notify( caststr )

  success_flag = ((error < tolerance))


  deallocate( realdata, fftdata, ifftdata )

end subroutine validate_fftw

