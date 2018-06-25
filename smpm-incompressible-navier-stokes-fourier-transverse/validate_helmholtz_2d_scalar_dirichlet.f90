subroutine validate_helmholtz_2d_scalar_dirichlet( tolerance, success_flag )
! Validates the GMRES implementation in 2D with dirichlet bcs using the 
! 2D Helmholtz Equation, which is given as,
!
!        d * Lap( phi ) - phi = f 
! 
!  where f is the RHS containing the density field, phi is the 
!  unknown scalar field such as the density and d is the coefficient containing 
!  the timestep and diffusion information. 
!
!  April 2017
!  Greg Thomsen & Gustavo Rivera

  use constants, only:               nsuby, nsg, nu_d
  use geom, only:                    cx, cz
  use options, only:                 gmres_maxit_viscous, gmres_restart_viscous, gmres_tol_viscous
  use parallel_linear_algebra, only: pmaxval, pminval
  use precision, only:               dp
  use timestepping, only:            g0, dt
  use woodbury_matrices, only:       rpk

  implicit none

  real(kind=dp), intent(in)                :: tolerance
  logical, intent(out)                     :: success_flag

  ! General use variables.
  real(kind=dp)                            :: pi, mm, kk

  real(kind=dp)                            :: linf_laplacian, linf_operator, &
                                              linf_solution

  ! Variables for GMRES.
  real(kind=dp)                            :: gmres_tol
  integer                                  :: gmres_diffusion_iterations
  integer                                  :: gmres_restart

  ! Variables for testing real/complex GMRES implementation.
  real(kind=dp), allocatable, dimension(:) :: lap, phi0, dummy1, phi, rhs, dummy2, dummy3

  character(len=256)                       :: caststr

  real(kind=dp)                            :: Lx, Lz

  ! Declare the functions used with GMRES to solve the viscous problems.
  external                                 :: apply_smpm_diffusion

  ! Allocate Arrays
  allocate( rhs(1:rpk) )
  allocate( dummy1(1:rpk), dummy2(1:rpk), dummy3(1:rpk) )
  allocate( lap(1:rpk) )
  allocate( phi0(1:rpk), phi(1:rpk) )

  ! Get pi.
  pi = 2.0_dp * acos( 0.0_dp )

  ! Get the domain lengths
  Lx = pmaxval( reshape( cx, (/rpk * nsuby/) ) ) - pminval( reshape( cx, (/rpk * nsuby/) ) )
  Lz = pmaxval( reshape( cz, (/rpk * nsuby/) ) ) - pminval( reshape( cz, (/rpk * nsuby/) ) )

  ! Set the timestepping coefficient value.
  g0 = 1.0_dp

  ! Generate phi using a single wavenumber in the transverse.
  kk = 2.0_dp * pi / Lx
  mm = 2.0_dp * pi / Lz
  phi0 = sin( kk * cx(:, 1) ) * sin( mm * cz(:, 1) )

  ! Compute the Laplacian.
  call compute_laplacian( phi0, lap )

  ! Check the Laplacian error.
  dummy1         = -(kk**2.0_dp + mm**2.0_dp) * sin( kk * cx(:, 1) ) * sin( mm * cz(:, 1) )
  linf_laplacian = pmaxval( abs( dummy1 - lap) )
  write( caststr, '(A,D17.10)' ) 'Linf Laplacian in 2D:', linf_laplacian
  call notify( caststr )

  ! Set GMRES parameters.
  gmres_tol                  = gmres_tol_viscous
  gmres_diffusion_iterations = min( gmres_maxit_viscous, nsg - 1 )
  gmres_restart              = min( gmres_restart_viscous, gmres_maxit_viscous )
  phi                        = 0.0_dp * phi0 * (1 - exp( -nu_d * dt ))

  ! Manufacture a RHS.
  call apply_smpm_diffusion( rhs, phi0 )

  ! Check error of RHS.
  dummy2        = nu_d * dt / g0 * dummy1  - phi0
  linf_operator = pmaxval( abs( dummy2 - rhs ) )
  write( caststr, '(A,D17.10)' ) 'Linf Linear Operator in 2D:', linf_operator
  call notify( caststr )

  ! We are solving the purely 2D Helmholtz (i.e. apply_smpm_diffusion as the operator) 
  rhs = rhs / g0
  call compute_gmres_householder( phi, rhs, rpk, gmres_tol, &
                                  gmres_diffusion_iterations, gmres_restart, &
                                  apply_smpm_diffusion )

  ! Check the solution's error.
  linf_solution = pmaxval( abs( phi - phi0 ) )
  write( caststr, '(A,D17.10)' ) 'Linf Solution Error in 2D:', linf_solution
  call notify( caststr )


  success_flag = (linf_solution  < tolerance)

   ! Reset the timestepping coefficient
   g0 = 11.0_dp/6.0_dp


end subroutine validate_helmholtz_2d_scalar_dirichlet
