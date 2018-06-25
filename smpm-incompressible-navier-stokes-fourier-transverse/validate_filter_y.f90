subroutine validate_filter_y( tolerance, success_flag )
! Validates the laplacian function in 2D and 3D using a sinusoidal 
! with a user-specified wavenumber. 

  use constants, only:               nky, nsuby
  use parallel_linear_algebra
  use precision, only:               dp
  use transverse, only:              cy
  use woodbury_matrices, only:       rpk

  implicit none

  real(kind=dp), intent(in)                      :: tolerance
  logical, intent(out)                           :: success_flag

  real(kind=dp)                                  :: error

  real(kind=dp), allocatable, dimension(:,:)     :: phi, phi_check, phi_filtered
  real(kind=dp), allocatable, dimension(:)       :: noise_real, noise_imag

  complex(kind=dp), allocatable, dimension(:)    :: noise
  complex(kind=dp), allocatable, dimension(:,:)  :: Fphi

  character(len=256)                             :: caststr

  real(kind=dp)                                  :: pi
  integer                                        :: ii

  allocate(Fphi(1:rpk,1:nky))
  allocate(phi(1:rpk,1:nsuby))
  allocate(phi_check(1:rpk,1:nsuby))
  allocate(phi_filtered(1:rpk,1:nsuby))
  allocate(noise(1:nky))
  allocate(noise_real(1:nky), noise_imag(1:nky))

  ! Define Pi
  pi = 2.0_dp * acos( 0.0_dp )

  ! Setup Filter 
  call setup_filter_y()

  ! Initialize the function
  phi = 0.0_dp

  ! Set function
  phi = exp(-(cy - pi)**2)

  ! Filter the function
  call apply_filter_y(phi, phi_filtered)

  ! Store the filtered function for error calculation
  phi_check = phi_filtered

  ! Compute inverse fourier transform of unfiltered function
  call apply_fft(phi, Fphi)

  ! Generate White Noise
  call random_number( noise_real )
  call random_number( noise_imag )

  ! Construct complex noise array
  do ii = 1,nky
     noise(ii) = cmplx( noise_real(ii), noise_imag(ii), kind=dp)
  enddo

  ! Add noise to the highest wavenumbers of the unfiltered function
  Fphi(:,nky)   = Fphi(:,nky)   + noise(nky)

  ! Bring the unfiltered function back to physical space
  call apply_ifft(Fphi,phi)

  ! Filter the noisy unfiltered function
  call apply_filter_y( phi,phi )

  error = pmaxval( reshape( abs( phi_check - phi ), (/rpk * nsuby/) )  )
  write( caststr, '(A,D17.10)' ) 'Error in filter:', error
  call notify( caststr )

  success_flag = ((error < tolerance))


end subroutine validate_filter_y
