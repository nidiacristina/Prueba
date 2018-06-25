subroutine apply_3D_poisson( Ax, x )
! Applies the Poisson-Neumann operator on x 

  use constants,only:          nky
  use options, only:           bc_flag_lhsgmres, facrobin_PPE
  use precision, only:         dp
  use transverse,only:         qy
  use woodbury_matrices, only: rpk

  implicit none

  real(kind=dp), dimension(1:rpk, 1:nky), intent(out)     :: Ax
  real(kind=dp), dimension(1:rpk, 1:nky), intent(in)      :: x
  real(kind=dp), dimension(1:rpk, 1:nky)                  :: q_x, q_z

  integer                                                 :: ii

  real(kind=dp), parameter                                :: delta_poisson = 1.0_dp

 
  ! Apply each wavenumber.
  do ii = 1, nky

     ! Apply the shifted Laplacian.
     call compute_gradient_and_laplacian( x(:, ii), Ax(:, ii), q_x(:, ii), q_z(:, ii) )
     Ax(:, ii) = Ax(:, ii) - real( qy(ii)**2, kind=dp ) * x(:, ii)

     ! Apply the boundary conditions.
     call apply_bc( x(:, ii), Ax( :, ii), bc_flag_lhsgmres, delta_poisson, q_x(:, ii), q_z(:, ii) )

     ! Apply the patching conditions.
     call apply_patching( x(:, ii), Ax( :, ii), delta_poisson, facrobin_PPE, q_x(:, ii), q_z(:, ii) )
  
  enddo

end subroutine apply_3D_poisson
