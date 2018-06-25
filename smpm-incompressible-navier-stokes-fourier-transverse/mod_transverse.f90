module transverse
! Contains variables related to the transverse directions (FFTs, transverse
! grid, and wavenumbers).

  use, intrinsic :: iso_c_binding

  use precision, only: dp

  implicit none
  save

  include 'fftw3.f03'

  real(kind=dp)                               :: Ly       ! Transverse dimension.
  real(kind=dp)                               :: delta_y  ! Transverse spacing
  real(kind=dp), allocatable, dimension(:, :) :: cy       ! Transverse grid.
  complex(kind=dp), allocatable, dimension(:) :: ky       ! All transverse wavenumbers.
  real(kind=dp), allocatable, dimension(:)    :: filter_y ! Transverse filter function
  complex(kind=dp), allocatable, dimension(:) :: qy       ! Just the transverse wavenumbers we calculate the solution on.


  ! FFTW plans and buffers.
  ! XXX: need to review real_buf_FFT and complex_buf_FFT
  type(C_PTR)                                 :: planFFT, planIFFT
  real(C_DOUBLE), pointer                     :: real_buff_FFT(:)
  complex(C_DOUBLE_COMPLEX), pointer          :: complex_buff_FFT(:)
  type(C_PTR)                                 :: p_real_buff_FFT, p_complex_buff_FFT

end module transverse
