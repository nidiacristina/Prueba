subroutine make_streamfunction( vx, vz, psi )
! Builds a streamfunction from velocity vectors.

  use constants, only:               nsg, nsuby
  use legendre, only:                filter_xz
  use options, only:                 gmres_maxit_poisson, gmres_tol_poisson
  use parallel_linear_algebra, only: pnorm2
  use precision, only:               dp
  use woodbury_matrices, only:       rpk

  implicit none

  real(kind=dp), dimension(1:rpk,1:nsuby), intent(in)  :: vx
  real(kind=dp), dimension(1:rpk,1:nsuby), intent(in)  :: vz
  real(kind=dp), dimension(1:rpk,1:nsuby), intent(out) :: psi

  real(kind=dp), dimension(1:rpk,1:nsuby)              :: Ax
  real(kind=dp), allocatable, dimension(:,:)           :: gmres_rhs
  real(kind=dp)                                        :: gmres_tol
  integer                                              :: gmres_maxit, gmres_restart
  character(len=64)                                    :: caststr

  external                                             :: apply_streamfunction_matrix

  ! Build the vorticity field.
  allocate( gmres_rhs(1:rpk,1:nsuby) )
  call compute_curl( vx, vz, gmres_rhs )
  gmres_rhs = -1 * gmres_rhs
  psi       = 1.0_dp

  ! Set some constants for the GMRES solver.
  gmres_tol     = gmres_tol_poisson
  gmres_maxit   = min( gmres_maxit_poisson, nsg - 1 )
  gmres_restart = gmres_maxit

  ! Filter the vorticity.
  call apply_filter_xz( gmres_rhs, filter_xz )

  ! Solve the Poisson-Dirichlet problem for the streamfunction.
  call compute_gmres_householder( psi, gmres_rhs, rpk, gmres_tol, gmres_maxit, &
                                  gmres_restart, apply_streamfunction_matrix )

  ! Check the residual.
  call apply_streamfunction_matrix( Ax, psi )
  write( caststr, '(D17.10)' )  pnorm2( reshape( Ax - gmres_rhs ,(/nsuby * rpk/) ) )

  ! Filter the streamfunction.
  call apply_filter_xz( psi, filter_xz  )

  ! Notify the user of some solver statistics.
  write( caststr, '(I10) ' ) gmres_maxit
  call notify( '      Poisson GMRES iterations         : ' // caststr )
  write( caststr, '(D17.10)' ) gmres_tol
  call notify( '      Poisson GMRES residual           : ' // caststr )
  call notify( '      Computed error in GMRES          : ' // caststr )

end subroutine make_streamfunction
