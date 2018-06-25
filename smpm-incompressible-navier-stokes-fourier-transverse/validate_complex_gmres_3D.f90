subroutine validate_complex_gmres_3D( output_prefix, tolerance, success_flag )
! Validates the complex-valued GMRES routine that supports 3D.

  use constants, only:               nky, nsuby, nu_d
  use geom, only:                    cx, cz
  use parallel_linear_algebra, only: pdot_product, pmaxval, pnorm2, znorm2
  use precision, only:               dp
  use timestepping, only:            g0, dt
  use transverse, only:              cy
  use woodbury_matrices, only:       rpk

  implicit none

  character(len=*), intent(in)                   :: output_prefix
  real(kind=dp), intent(in)                      :: tolerance
  logical, intent(out)                           :: success_flag

  ! General use variables.
  real(kind=dp)                                  :: pi, mm, ll, kk
  real(kind=dp)                                  :: rhs_error

  real(kind=dp)                                  :: linf_laplacian, linf_operator, &
                                                    linf_diffusion

  ! Variables for GMRES.
  real(kind=dp)                                  :: gmres_tol
  integer                                        :: gmres_diffusion_iterations

  ! Variables for testing real/comples GMRES implementation.
  real(kind=dp), allocatable, dimension(:, :)    :: lap, phi0, phi, rhs, dummy1, dummy2, dummy3
  complex(kind=dp), allocatable, dimension(:, :) :: Fphi0, Frhs

  character(len=256)                             :: caststr

  ! Declare the functions used with GMRES to solve the viscous problems.
  external                                       :: apply_3D_diffusive_fourier

  ! Allocate arrays.
  allocate( phi0(1:rpk, 1:nsuby), phi(1:rpk, 1:nsuby) )
  allocate( dummy1(1:rpk, 1:nsuby), dummy2(1:rpk, 1:nsuby), dummy3(1:rpk, 1:nsuby) )
  allocate( rhs(1:rpk, 1:nsuby) )
  allocate( Fphi0(1:rpk, 1:nky), Frhs(1:rpk, 1:nky) )
  allocate( lap(1:rpk,1:nsuby) )

  ! Get pi.
  pi = 2.0_dp * acos( 0.0_dp )

  ! Reset timestepping coefficient.
  g0 = 1.0_dp

  mm = 8.0_dp * pi
  ll = 8.0_dp * pi
  kk = 8.0_dp * pi
  phi0 = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )

  ! Compute the Laplacian.
  call compute_3D_laplacian( phi0, lap )

  ! Check error on Laplacian.
  dummy1         = -(kk**2.0_dp + ll**2.0_dp + mm**2.0_dp) * sin(kk * cx ) * sin( ll * cy ) * sin( mm * cz )
  linf_laplacian = maxval( abs( dummy1 - lap ) )
  write( caststr, '(A,D17.10)' ) 'Linf Error in Laplacian for 3D Complex GMRES', linf_laplacian
  call notify( caststr )

  ! Generate the right hand side.
  call apply_fft( phi0, Fphi0 )
  call apply_3D_diffusive_fourier( Frhs, Fphi0 )
  call apply_ifft( Frhs, rhs )

  ! Check linear operator error.
  dummy2        = nu_d * dt / g0 * dummy1 - phi0
  linf_operator = maxval( abs( dummy2 - rhs ) )
  write( caststr, '(A,D17.10)' ) 'Linf Error in Linear Operator for 3D Complex GMRES', linf_operator
  call notify( caststr )

  ! Write operator to file.
  open( 65, file=trim( output_prefix ) // '_operator_complex.txt' )
  write( 65, * ) rhs
  close( 65 )

  ! Dirichlet boundary conditions.
!   dummy3 = - ( nu_d * dt / g0 ) * (kk**2.0_dp + ll**2.0_dp + mm**2.0_dp) * sin( kk * cx) * sin( ll * cy) * sin( mm * cz) &
!          - sin( kk * cx ) * sin( ll * cy) * sin( mm * cz)

  ! Neumann boundary conditions.
!  dummy3 = - ( nu_d * dt / g0 ) * (mm**2.0_dp + ll**2.0_dp + kk**2.0_dp) * cos( mm * cx) * cos( ll * cy) * cos( kk * cz) &
!         - cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz )

  ! Solve System via GMRES
  call solve_3D_diffusion_complex( rhs, phi, gmres_tol, gmres_diffusion_iterations, rhs_error )

  ! Write solution to file.
  open( 65, file=trim( output_prefix ) // '_gmres_complex.txt' )
  write( 65, * ) phi0
  write( 65, * ) phi
  write( 65, * ) rhs
  close( 65 )

  linf_diffusion = maxval( abs( phi0 + phi ) )
  write( caststr, '(A,D17.10,A,D17.10,A,D17.10,A,I0)' ) &
          'Linf error of phi (Complex):', linf_diffusion, &
          'RHS Error:', rhs_error, &
          'Convergence:', gmres_tol, &
          'Iterations:', gmres_diffusion_iterations
  call notify( caststr )

  success_flag = ((linf_laplacian < tolerance) .and. &
                  (linf_operator  < tolerance) .and. &
                  (linf_diffusion < tolerance))

end subroutine validate_complex_gmres_3D

subroutine solve_3D_diffusion_complex( rhs, phi, gmres_tol, gmres_diffusion_iterations, rhs_error )
! This routine solves the 3D Helmholtz Eqn, which is used for the diffusion step
! The Helmholtz Equation is stated as:
!
!    alpha * Lap f - f = phi

  use constants, only:         nsg, nsuby, nky, nu_d
  use options, only:           gmres_maxit_viscous, gmres_tol_viscous
  use precision, only:         dp
  use timestepping, only:      dt, g0
  use woodbury_matrices, only: rpk

  implicit none

  ! The input/output variables.
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(in)  :: rhs
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out) :: phi
  real(kind=dp), intent(out)                            :: rhs_error
  real(kind=dp), intent(out)                            :: gmres_tol
  integer, intent(out)                                  :: gmres_diffusion_iterations

  ! Arrays for storing Fourier transform.
  complex(kind=dp), allocatable, dimension(:, :)        :: Frhs, Fphi, Frhs_check
  complex(kind=dp), allocatable, dimension(:)           :: Fbuffer

  ! Other variables.
  integer                                               :: gmres_restart
  real(kind=dp), dimension(1:rpk, 1:nsuby)              :: rhs_check

  ! Declare the function used with GMRES to solve the capacitance problem (i.e. Linear Operator)
  external                                              :: apply_3D_diffusive_fourier

  allocate( Fbuffer(1:rpk) )
  allocate( Frhs(1:rpk, 1:nky), Fphi(1:rpk, 1:nky), Frhs_check(1:rpk, 1:nky) )

  ! Initialize.
  Frhs = 0.0_dp
  Fphi = 0.0_dp

  ! Set GMRES tolerance.
  gmres_tol                  = gmres_tol_viscous
  gmres_diffusion_iterations = min( gmres_maxit_viscous, nsg - 1 )
  gmres_restart              = gmres_diffusion_iterations

  ! Step 1: Fourier transform of RHS.
  call apply_fft( rhs , Frhs )
  Fphi = 0.0_dp * Frhs * (1 - exp( -nu_d * dt )) ! This initialization converges faster and
                                                ! yields a more accurate result

  ! Step 2.1: Solve the linear system via complex GMRES.
  call compute_gmres_householder_complex( Fphi , -Frhs / g0, rpk * nky, &
                                          gmres_tol, gmres_diffusion_iterations, &
                                          gmres_restart, apply_3D_diffusive_fourier )

  ! Step 2.2: Check the error by recomputing RHS.
  call apply_3D_diffusive_fourier(Frhs_check, Fphi)
  call apply_ifft(Frhs_check, rhs_check)
  rhs_error = maxval( abs( rhs_check + rhs / g0 ) )

  ! Step 3: Inverse transform solution buffer.
  call apply_ifft( Fphi, phi )

end subroutine solve_3D_diffusion_complex

subroutine compute_3D_laplacian( u, Lu )
! Computes the Laplacian Lu of the scalar field u.  In computing the
! Laplacian, the gradient operator is implicitly computed, so this subroutine
! also returns the gradient in components (u_x,u_y,u_z).
!
! NOTE: This is a stripped down version of compute_3D_gradient_and_laplacian().

  use constants, only:             nky, nsuby
  use precision, only:             dp
  use woodbury_matrices, only:     rpk
  use transverse, only:            qy

  implicit none

  ! Input variables.
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out) :: Lu
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(in)  :: u

  ! Internal variables.
  integer                                               :: ii
  complex(kind=dp)                                      :: i
  complex(kind=dp), allocatable, dimension(:, :)        :: Fu
  real(kind=dp), allocatable, dimension(:, :)           :: u_yy

  ! Set the complex unit.
  i = cmplx( 0.0_dp, 1.0_dp, kind=dp )

  allocate( Fu(1:rpk, 1:nky) )
  allocate( u_yy(1:rpk, 1:nsuby) )

  ! Compute the 2D gradient and Laplacian plane-by-plane.
  do ii = 1, nsuby
     call compute_laplacian( u(:, ii), Lu(:, ii) )
  enddo

  ! Compute the transverse first and second derivatives.
  call apply_fft( u, Fu )
  do ii = 1, nky
     Fu(:, ii) = i * i * qy(ii) * qy(ii) * Fu(:, ii)
  enddo
  call apply_ifft( Fu, u_yy )

  ! Add the second derivative to the Laplacian.
  Lu = Lu + u_yy

end subroutine compute_3D_laplacian
