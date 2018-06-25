subroutine apply_ifft( q, Fq )
! Applies the Poisson-Neumann operator on x.

  use, intrinsic :: iso_c_binding
  use constants,         only: nky, nsuby
  use precision,         only: dp
  use transverse,        only: complex_buff_FFT, planIFFT, real_buff_FFT
  use woodbury_matrices, only: rpk

  implicit none

  include 'fftw3.f03'

  complex(kind=dp), dimension(1:rpk, 1:nky), intent(in)    :: q
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out)    :: Fq

  integer                                                  :: ii

  ! Loop over all the points taking FFTS.
  do ii = 1, rpk

     ! NOTE: FFTW computes unnormalized FFTs so we have to scale the output
     !       by the number of points.
     complex_buff_FFT = q(ii, :)
     call FFTW_EXECUTE_DFT_C2R( planIFFT, complex_buff_FFT, real_buff_FFT )
     Fq(ii, :) = real_buff_FFT / nsuby

  enddo

end subroutine apply_ifft
