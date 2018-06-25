subroutine setup_transverse
! Constructs the transverse grid, the transverse wavenumbers, and the
! tranvserse FFTs.

  use, intrinsic :: iso_c_binding

  use constants, only:  nky, nsuby
  use precision, only:  dp
  use transverse, only: real_buff_FFT, complex_buff_FFT, ky, Ly, cy, &
                        planFFT, planIFFT, qy, delta_y

  implicit none

  include 'fftw3.f03'

  integer                       :: ii, ndx
  real(kind=dp)                 :: pi

  ! Get pi.
  pi = 2.0_dp * acos( 0.0_dp )

  ! Build the transverse wavenumber grid.
  if (nsuby > 1) then
     delta_y = cy(1, 2) - cy(1, 1)
     Ly      = cy(1, nsuby) - cy(1, 1) + delta_y
     ndx     = 1

     ! Build the transverse wavenumber grid from [-Ny/2 , Ny/2 -1]
     do ii = -nsuby / 2, nsuby / 2 - 1
        ky(ndx) = cmplx( 2 * pi * ii / Ly, kind=dp )
        ndx     = ndx + 1
     enddo

     ! Build the transverse wavenumber grid for the wavenumbers we're actually
     ! calculating.
     do ii = 1, nky
        qy(ii) = cmplx( 2.0_dp * pi * (ii - 1) / Ly, kind=dp )
     enddo

  else

     ! Set critical parameters in order to collapse 3D SMPM to 2D SMPM
     cy      = 0.0_dp
     ky(1)   = cmplx( 0, kind=dp )
     qy(1)   = cmplx( 0.0_dp, kind=dp )
     delta_y = 0.0_dp
  end if



  ! Set up plans for the ffts, reversing the order.
  planFFT  = fftw_plan_dft_r2c_1d( int( nsuby, kind=C_INT ), real_buff_FFT, complex_buff_FFT, FFTW_ESTIMATE )
  planIFFT = fftw_plan_dft_c2r_1d( int( nsuby, kind=C_INT ), complex_buff_FFT, real_buff_FFT, FFTW_ESTIMATE )

end subroutine setup_transverse
