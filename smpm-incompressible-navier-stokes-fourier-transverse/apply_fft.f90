subroutine apply_fft( q, Fq )
! Applies the Poisson-Neumann operator on x.

  use, intrinsic :: iso_c_binding
  use constants, only:         nky, nsuby
  use precision, only:         dp
  use transverse, only:        complex_buff_FFT, planFFT, real_buff_FFT
  use woodbury_matrices, only: rpk

  implicit none

  include 'fftw3.f03'

  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(in)   :: q
  complex(kind=dp), dimension(1:rpk, 1:nky), intent(out) :: Fq

  integer                                                :: ii

  ! Loop over all the points taking FFTs.
  do ii = 1, rpk

     real_buff_FFT = q(ii, :)
     call FFTW_EXECUTE_DFT_R2C( planFFT, real_buff_FFT, complex_buff_FFT )
     Fq(ii, :) = complex_buff_FFT

  enddo

end subroutine apply_fft
