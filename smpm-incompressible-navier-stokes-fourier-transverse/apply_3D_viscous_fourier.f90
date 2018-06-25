subroutine apply_3D_viscous_fourier( Aq, q )
! Applies the Poisson-Neumann operator on x in fourier space using real data type
! (aka real and imaginary parts stored separately on real datatype arrays)

  use constants,only:          nky, nu
  use options, only:           bc_flag_viscous, facrobin
  use precision, only:         dp
  use timestepping, only:      dt, g0
  use transverse,only:         qy
  use woodbury_matrices, only: rpk

  implicit none

  real(kind=dp), dimension(1:3*rpk, 1:nky), intent(out)   :: Aq
  real(kind=dp), dimension(1:3*rpk, 1:nky), intent(in)    :: q
  real(kind=dp), dimension(1:3*rpk)                       :: real_buff
  real(kind=dp), dimension(1:3*rpk)                       :: q_x, q_z
  integer                                                 :: ii, jj, a, z

  real(kind=dp), dimension(1:2*rpk)                       :: qx_and_qz
  real(kind=dp), dimension(1:2*rpk)                       :: Aqx_and_Aqz
  real(kind=dp), dimension(1:2*rpk)                       :: real_dudx
  real(kind=dp), dimension(1:2*rpk)                       :: real_dudz

  real(kind=dp)                                           :: delta

  real(kind=dp), parameter                                :: delta_poisson = 1.0_dp

  ! Set the viscosity to the value specified by the input file.
  delta = nu * dt / g0 ! gamma0 because of KIO time-splitting.
 
  ! Loop over each wavenumber.
  do ii = 1, nky

     ! Loop over the three components of velocity.
     do jj = 1, 3

        ! Set array bounds for this component.
        a = (jj - 1) * rpk + 1
        z = (jj - 1) * rpk + rpk

        ! Apply the shifted Laplacian to the real and imaginary parts
        ! separately.
        call compute_gradient_and_laplacian(  q(a:z, ii) , real_buff(a:z), q_x(a:z), q_z(a:z) )
        Aq(a:z, ii) = delta * real_buff(a:z) - delta * real( qy(ii)**2, kind=dp ) * q(a:z, ii) - q(a:z, ii)

        ! Apply the patching conditions.
        real_buff(a:z) =   Aq(a:z, ii)
        call apply_patching( q(a:z, ii) , real_buff(a:z), delta, facrobin, q_x(a:z), q_z(a:z) )
        Aq(a:z, ii) = real_buff(a:z)

     enddo

     ! Set up arrays to hold the real/imaginary parts of the x and z
     ! derivatives of the ux and uz velocity fields.

     ! Real part of the x derivative of ux and uz.
     real_dudx(    1:rpk    ) = q_x(       1:rpk   )
     real_dudx( rpk+1:2*rpk ) = q_x( 2*rpk+1:3*rpk )

     ! Real part of the z derivative of ux and uz.
     real_dudz(     1:rpk   ) = q_z(       1:rpk )
     real_dudz( rpk+1:2*rpk ) = q_z( 2*rpk+1:3*rpk )

     ! Set up an array that contains ux and uz.
     qx_and_qz(     1:rpk   ) = q(       1:rpk, ii   )
     qx_and_qz( rpk+1:2*rpk ) = q( 2*rpk+1:3*rpk, ii )

     ! Set up an array to hold the viscous operator applied to the x and z
     ! components.
     Aqx_and_Aqz(     1:rpk   ) =  Aq(       1:rpk, ii   )
     Aqx_and_Aqz( rpk+1:2*rpk ) =  Aq( 2*rpk+1:3*rpk, ii )

     ! Apply the vector viscous boundary conditions in (ux,uz).
     call apply_vector_viscous_bc( qx_and_qz, Aqx_and_Aqz, real_dudx, real_dudz, bc_flag_viscous, delta )

     ! Store the result in the output vector in the right place.
     Aq(       1:rpk, ii   ) = Aqx_and_Aqz(     1:rpk )
     Aq( 2*rpk+1:3*rpk, ii ) = Aqx_and_Aqz( rpk+1:2*rpk )

  enddo

end subroutine apply_3D_viscous_fourier


