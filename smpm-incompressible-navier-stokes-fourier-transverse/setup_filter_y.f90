subroutine setup_filter_y
! Generates the vector for spectral filtering in the transverse,
! meaning the y-direction.

  use constants, only: nsuby, nky
  use transverse, only: filter_y, Ly, qy
  use options, only: filter_order_y
  use precision, only: dp

  implicit none


  ! Local variables.
  integer                         :: ii
  real(kind=dp)                   :: alpha, amp
  real(kind=dp)                   :: akc, akmax
  real(kind=dp)                   :: kcf
  real(kind=dp)                   :: pi
  real(kind=dp)                   :: theta
  real(kind=dp)                   :: wavenumber 

  ! Initialize the filter buffer
  filter_y = 0.0_dp

  ! Set fraction of maximum wave# above which filter is applied
  kcf = 0.0_dp

  ! Define Pi
  pi = 2.0_dp * acos( 0.0_dp )

  ! Check presence of transverse domain
  if (nsuby > 1) then

     ! Maximum Wavenumber Value
     alpha = 2.0_dp * pi / Ly

     ! Highest Resolved Wavenumber
     akmax = (real(nky, kind=dp) - 1) * alpha

     ! Set the Filter amplitude
     amp = -log( 1.0d-16 )

     ! Set the Filter Lag
     akc = akmax * kcf
 
     ! Loop across each wavenumber
     do ii = 1,nky

        ! Set the Wavenumber
        wavenumber = real(qy(ii), kind=dp)

        ! Construct Filter
        if ( wavenumber .ge. akc ) then
           theta = ( wavenumber - akc ) / (akmax - akc)
           filter_y(ii) = exp( -amp * theta**filter_order_y )
        else
           filter_y(ii) = 1.0_dp
        endif

        ! Reset buffers
        theta = 0.0_dp
        wavenumber = 0.0_dp

     enddo

  else

      ! If no transverse domain, leave filter at 1
      filter_y(1) = 1.0_dp

  endif

end subroutine setup_filter_y
