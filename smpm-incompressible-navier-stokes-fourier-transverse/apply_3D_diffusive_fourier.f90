subroutine apply_3D_diffusive_fourier( Aq, q )
! Applies the Poisson-Neumann operator on x in fourier space using real data type
! (aka real and imaginary parts stored separately on real datatype arrays)

  use constants,only:                nu_d, nky
  use options, only:                 bc_flag_diffusion, facrobin
  use parallel_linear_algebra, only: pdot_product, pmaxval, pnorm2, znorm2
  use precision, only:               dp
  use timestepping, only:            dt, g0
  use transverse,only:               qy
  use woodbury_matrices, only:       rpk

  implicit none

  real(kind=dp), dimension(1:rpk, 1:nky), intent(out)  :: Aq
  real(kind=dp), dimension(1:rpk, 1:nky), intent(in)   :: q
  real(kind=dp), dimension(1:rpk)                      :: real_buff, Rq_x, Rq_z
  real(kind=dp)                                        :: delta

  integer                                              :: ii

  ! Set the viscosity to the value specified by the input file.
  delta = nu_d * dt / g0 ! gamma0 because of KIO time-splitting.

  do ii = 1, nky

     ! Apply the shifted Laplacian to the real and imaginary parts separately.
     call compute_gradient_and_laplacian(  q(:, ii) , real_buff, Rq_x, Rq_z )
     Aq(:, ii) = delta * real_buff - ( delta * real( qy(ii)**2, kind=dp ) * q(:, ii) ) - q(:, ii)

     ! Set up an array to hold the diffusive operator applied to the x and z
     ! components.
     real_buff =   Aq(:, ii)

     ! Apply the patching conditions.
     call apply_patching( q(:, ii) , real_buff, delta, facrobin, Rq_x, Rq_z )

     ! Apply the diffusive boundary conditions in rho for both real and
     ! imaginary components.
     call apply_bc(  q(:, ii) , real_buff, bc_flag_diffusion, delta,  Rq_x, Rq_z )

     ! Store the result in the output vector in the right place.
     Aq(:, ii)  = real_buff

  enddo

end subroutine apply_3D_diffusive_fourier
