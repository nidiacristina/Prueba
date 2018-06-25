subroutine apply_filter_y( q, q_filtered )
! Applies the exponential filter in the transverse direction, similarly 
! to the scheme found in smpm-fourier-fourier
!
!
! April 2017
! Gustavo Rivera

  use constants, only:          nky, nsuby
  use precision, only:          dp
  use transverse, only:         filter_y
  use woodbury_matrices, only:  rpk

  implicit none

  ! Input/Output buffers
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(in)  :: q
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out) :: q_filtered

  ! Local Buffers
  complex(kind=dp), dimension(1:rpk, 1:nky)             :: Fq
  !real(kind=dp), dimension(1:rpk, 1:nky)                :: Fq_real, Fq_imag

  ! Local variables.
  integer                                               :: ii, kk
  real(kind=dp)                                         :: cosc, sinc
  real(kind=dp)                                         :: filter_y_magnitude
  real(kind=dp)                                         :: Fq_mag, Fq_real, Fq_imag  
  real(kind=dp)                                         :: phase

  ! Initialize Buffers
  Fq                 = cmplx( 0.0_dp,0.0_dp,kind=dp )

  ! Initialize Variables
  filter_y_magnitude = 0.0_dp 
  Fq_real            = 0.0_dp
  Fq_imag            = 0.0_dp
  Fq_mag             = 0.0_dp
  cosc               = 0.0_dp
  sinc               = 0.0_dp
  phase              = 0.0_dp

  ! Check if there is a transverse direction
  if (nsuby .GT. 1) then

     ! Fourier Transform the Array
     call apply_fft( q, Fq )  

     ! Begin loop over wavenumbers
     do kk = 1,nky
 
        ! Begin loop over physical domain
        do ii = 1,rpk

           ! Pass complex Fouerir Coefficients to real variables
           Fq_real = real(  Fq( ii,kk ), kind=dp )
           Fq_imag = aimag( Fq( ii,kk ) )

           ! Compute Magnitude of Fourier Coefficients
           Fq_mag = sqrt( Fq_real**2.0_dp + Fq_imag**2.0_dp )

           ! Check for zeroth mode
           if ( Fq_mag .NE. 0.0_dp ) then
            
              ! Normalize Fourier Coefficients
              cosc = Fq_real/Fq_mag
              sinc = Fq_imag/Fq_mag
              if (sinc .GE. 0.0_dp) phase =  acos(cosc)
              if (sinc .LT. 0.0_dp) phase = -acos(cosc)

              ! Incorporate magnitude to transverse filter
              filter_y_magnitude = Fq_mag * filter_y(kk)
           
              ! Apply filter function to Fourier Coefficients
              Fq_real = filter_y_magnitude * cos(phase)
              Fq_imag = filter_y_magnitude * sin(phase)
    
           else

              Fq_real = 0.0_dp
              Fq_imag = 0.0_dp
 
           endif !End if statement for zeroth mode

           ! Reconstruct Function in Fourier Space
           Fq( ii,kk ) = cmplx( Fq_real,Fq_imag, kind=dp )
        
           ! Reset variables
           filter_y_magnitude = 0.0_dp 
           Fq_imag            = 0.0_dp
           Fq_real            = 0.0_dp
           Fq_mag             = 0.0_dp
           cosc               = 0.0_dp
           sinc               = 0.0_dp
           phase              = 0.0_dp

        enddo !End Loop over physical domain

     enddo !End loop over fourier space
 
     ! Inverse Fourier Transform the array
     call apply_ifft( Fq,q_filtered )

  else

     q_filtered = q

  endif !End if statement for transverse direction check


end subroutine apply_filter_y
