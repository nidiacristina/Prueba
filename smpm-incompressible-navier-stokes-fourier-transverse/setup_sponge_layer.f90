subroutine setup_sponge_layer
! This subroutine generates the sponge layer at the beginning (x = 0) and end
! (x = Lx) of the domain.  It is based on the raycoeff() function in Peter
! Diamessis' Fourier-Fourier code.  The method implemented here is detailed
! in Amar's thesis. 

  use constants, only:                nsg, nsuby, rank, root
  use geom, only:                     cx
  use parallel_linear_algebra, only:  pmaxval, pminval
  use precision, only:                dp
  use sponge, only:                   left_fraction, raycoeff, &
                                      right_fraction, sponge_layer_location, &
                                      time_scale_in_x
  use woodbury_matrices, only:        rpk

  implicit none

  ! Loop and range indices.
  integer                                 :: kk, ll

  ! Left and right thicknesses of the sponge layer.
  real(kind=dp)                           :: left_thickness, right_thickness

  ! General use variables.
  real(kind=dp)                           :: pi

  ! Sponge layer variables.
  real(kind=dp)                           :: eta0, eta, xx
  real(kind=dp)                           :: beginning, domain_length, end

  ! Full domain buffer to output the sponge layer
  real(kind=dp), dimension(1:nsg,1:nsuby) :: raycoeff_full

  ! Get pi.
  pi = 2.0_dp * acos( 0.0_dp )

  ! NOTE: These maximum and minimum are based on a coordinate system sitting
  !       at the free surface.  We have to make sure that if this is not the
  !       case the values are revised.

  ! Throw some error messages to ensure that the sponge layer has been properly requested
  !
  ! Check for location of sponge
  if (( sponge_layer_location(2) == 0 ) .or. ( sponge_layer_location(4) == 0)) then
     call notify( ' No Sponge Locations Specified. Stopping Run ')
     stop
  endif

  ! Check for time scale of the sponge
  if (time_scale_in_x .eq. 0.0_dp) then
     call notify( ' No Sponge Timescale Specified. Stopping Run ')
     stop
  endif

  ! Determine the extremities of the domain in the x-coordinate
  end           = pmaxval( reshape( cx , (/rpk * nsuby/) ))
  beginning     = pminval( reshape( cx , (/rpk * nsuby/) ))
  domain_length = end - beginning

  ! Set up the sponge layer thickness at left and right domain edges.
  left_thickness  = left_fraction  * domain_length
  right_thickness = right_fraction * domain_length

  ! Loop over the transverse planes.
  do ll = 1, nsuby

     do kk = 1, rpk

       ! Grab the coordinate.
       xx = cx(kk,ll)

       ! Handle the left wall...
       if ( xx .lt. beginning + left_thickness .and. sponge_layer_location(4) == 1 ) then
          eta0             = beginning + left_thickness
          eta              = pi/2.0_dp * (xx - eta0) / left_thickness
          raycoeff(kk, ll) = time_scale_in_x * sin( eta )**2
       endif

       ! ... and the right wall.
       if ( xx .gt. end - right_thickness .and. sponge_layer_location(2) == 1 ) then
          eta0             = end - right_thickness
          eta              = pi/2.0_dp * (xx - eta0) / right_thickness
          raycoeff(kk, ll) = time_scale_in_x * sin( eta )**2
       endif

     enddo

  enddo

  ! Write Output to file for debug
  call gather_3D_array( raycoeff, raycoeff_full)

  if (rank .eq. root) then
     open(65, file='sponge.txt')
     write(65,*) raycoeff_full
     close(65)
  endif

end subroutine setup_sponge_layer
