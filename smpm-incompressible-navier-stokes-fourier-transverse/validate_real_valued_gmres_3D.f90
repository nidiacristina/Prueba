subroutine validate_real_valued_gmres_3D( output_prefix, tolerance, success_flag )
! Validates the real-valued GMRES routine that supports 3D.

  use constants, only:               nsuby, nky, nu_d
  use geom, only:                    cx, cz
  use parallel_linear_algebra, only: pdot_product, pmaxval, pnorm2, znorm2
  use precision, only:               dp
  use timestepping, only:            g0, dt
  use transverse, only:              cy, qy
  use validation, only:              wave
  use woodbury_matrices, only:       rpk

  implicit none

  character(len=*), intent(in)                   :: output_prefix
  real(kind=dp), intent(in)                      :: tolerance
  logical, intent(out)                           :: success_flag

  ! General use variables.
  real(kind=dp)                                  :: pi, mm, ll, kk

  real(kind=dp)                                  :: linf_operator, linf_diffusion

  ! Variables for gmres.
  integer                                        :: ii

  ! Define variables for testing real/complex GMRES implementation.
  real(kind=dp), allocatable, dimension(:, :)    :: phi0, phi, rhs, dummy1, dummy2, dummy3
  complex(kind=dp), allocatable, dimension(:, :) :: Fphi0, Frhs

  real(kind=dp), allocatable, dimension(:)       :: real_buff, imag_buff

  character(len=256)                             :: caststr

  ! Declare the functions used with GMRES to solve the viscous problems.
  external                                       :: apply_smpm_3D_diffusion
  external                                       :: apply_smpm_3D_diffusion_real_and_complex
  external                                       :: apply_smpm_3D_viscous
  external                                       :: apply_smpm_diffusion
  external                                       :: apply_3D_Diffusive_fourier

  ! Allocate matrices.
  allocate( phi0(1:rpk, 1:nsuby), phi(1:rpk, 1:nsuby) )
  allocate( dummy1(1:rpk, 1:nsuby), dummy2(1:rpk, 1:nsuby), dummy3(1:rpk, 1:nsuby) )

  ! Allocate diffusion statistics arrays.
  allocate( rhs(1:rpk,1:nsuby) )
  allocate( real_buff(1:rpk) )
  allocate( imag_buff(1:rpk) )

  allocate( Fphi0(1:rpk, 1:nky), Frhs(1:rpk, 1:nky) )

  ! Get pi.
  pi = 2.0_dp * acos( 0.0_dp )

  ! Reset the timestepping coefficient.
  g0 = 1.0_dp

  mm = 8.0_dp * pi
  ll = 8.0_dp * pi
  kk = 8.0_dp * pi
  phi0 = sin( mm * cx ) * sin( ll * cy ) * sin( kk * cz )

  ! Build rhs
  call apply_fft( phi0, Fphi0 )
  do ii = 1, nky

     ! Set the wavenumber for the linear operators below.
     wave = qy(ii)

     ! Apply the operator twice, one for real and imaginary.
     call apply_smpm_3D_diffusion_real_valued( real_buff, real( Fphi0(:, ii), kind=dp ) )
     call apply_smpm_3D_diffusion_real_valued( imag_buff, aimag( Fphi0(:, ii) ) )

     ! Store real and imaginary buffers into complex RHS buffer.
     Frhs(:, ii) = cmplx( real_buff, imag_buff, kind=dp )

  enddo
  call apply_ifft( Frhs, rhs )

  ! Check the linear operator's error.
  dummy2 = (-(nu_d * dt / g0) * (mm**2.0_dp + ll**2.0_dp + kk**2.0_dp) * &
            sin( mm * cx ) * sin( ll * cy ) * sin( kk * cz ) - phi0)

  ! Neumann boundary condition.
  !dummy3 = (-(nu_d * dt / g0) * (mm**2.0_dp + ll**2.0_dp + kk**2.0_dp) * &
  !          cos( mm * cx ) * cos( ll * cy ) * cos( kk * cz ) - phi0)

  ! Write the operator to file.
  open( 65, file=trim( output_prefix ) // '_opeator_real_valued.txt' )
  write( 65, * ) rhs
  close( 65 )

  ! Display operator error.
  linf_operator = maxval( abs( dummy2 - rhs ) )
  write( caststr, '(A,D17.10)' ) 'Linf error in Linear Operator 3D (real-valued):', linf_operator
  call notify( caststr )

  ! Apply the boundary condition.
  !call apply_bc( phi0, rhs, bc_flag_diffusion, nu_d * dt / g0, dummy1, dummy2 )

  ! Apply the patching condition.
  !call apply_patching( phi0, rhs, nu_d * dt / g0, facrobin, dummy1, dummy2)

  call solve_3D_diffusion_GAR( rhs, phi )

  ! Write the solution to disk.
  open( 65, file=trim( output_prefix ) // '_gmres_3D_real_valued.txt' )
  write( 65, * ) phi0
  write( 65, * ) phi
  write( 65, * ) rhs
  close( 65 )

  ! Display the error in the solution.
  linf_diffusion = maxval( abs( phi0 - phi ) )
  write( caststr, '(A,D17.10)' ) 'Linf Error of phi (real-valued):', linf_diffusion
  call notify( caststr )

  success_flag = ((linf_operator  < tolerance) .and. &
                  (linf_diffusion < tolerance))

end subroutine validate_real_valued_gmres_3D

subroutine solve_3D_diffusion_GAR( rhs, phi )
! This routine solves the 3D Helmholtz Eqn, which is used for the diffusion step
! The Helmholtz Equation is stated as: alpha * Lap f - f = phi

  use constants, only:         nsg, nky, nsuby
  use options, only:           gmres_maxit_viscous, gmres_tol_viscous
  use precision, only:         dp
  use transverse, only:        qy
  use validation, only:        wave
  use woodbury_matrices, only: rpk

  implicit none

  ! Set the input/output variables
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(in)      :: rhs
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out)     :: phi
  real(kind=dp), dimension(1:nky)                           :: real_tol, imag_tol
  integer, dimension(1:nky)                                 :: real_iters, imag_iters

  ! Arrays for storing fourier transform
  complex(kind=dp), dimension(1:rpk,1:nky)                  :: Frhs, Fphi
  real(kind=dp), dimension(1:rpk)                           :: real_buff, imag_buff
  real(kind=dp), dimension(1:rpk)                           :: rhs_real, rhs_imag

  ! Other Variables
  real(kind=dp)                                             :: real_error, imag_error
  real(kind=dp)                                             :: gmres_tol
  integer                                                   :: gmres_diffusion_iterations
  integer                                                   :: gmres_restart

  ! Looping variable
  integer                                                   :: ii

  character(len=256)                                        :: caststr

  ! Declare the function used with GMRES to solve the capacitance problem (i.e. Linear Operator)
  external                                                  :: apply_smpm_3D_diffusion_real_valued

  ! Initialize arrays.
  Frhs      = 0.0_dp
  Fphi      = 0.0_dp
  imag_buff = 0.0_dp
  real_buff = 0.0_dp

  ! Set GMRES tolerance.
  gmres_tol                  = gmres_tol_viscous
  gmres_diffusion_iterations = min( gmres_maxit_viscous, nsg - 1 )
  gmres_restart              = gmres_diffusion_iterations

  ! Step 1: Take Fourier Transform of RHS
  call apply_fft( rhs, Frhs )

  ! Step 2: Loop over each wavenumber.
  !
  ! NOTE: Due to the choice of wavenumber (i.e. 8pi), only location 5 is active.
  ! XXX: This isn't actually true and needs to be fixed.
  do ii = 5,5

     ! Step 2.1: Set the wavenumber.
     wave = qy(ii)

     ! Step 2.2: Solve the real linear system via GMRES.
     call compute_gmres_householder( real_buff, real( Frhs(:, ii), kind=dp ), rpk, &
                                     gmres_tol, gmres_diffusion_iterations, &
                                     gmres_restart, apply_smpm_3D_diffusion_real_valued )

     ! Step 2.2.1: Store the real global tolerance and iterations
     real_tol(ii)   = gmres_tol
     real_iters(ii) = int(gmres_diffusion_iterations)

     ! Step 2.2.2: Reset the GMRES tolerance for imaginary solve.
     gmres_tol                  = gmres_tol_viscous
     gmres_diffusion_iterations = min( gmres_maxit_viscous, nsg - 1 )
     gmres_restart              = gmres_diffusion_iterations

     ! Step 2.3: Solve the imaginary linear system via GMRES.
     call compute_gmres_householder( imag_buff , aimag( Frhs(:, ii) ), rpk, &
                                     gmres_tol, gmres_diffusion_iterations, &
                                     gmres_restart, apply_smpm_3D_diffusion_real_valued )

     ! Step 2.3.1: Store the imaginary global tolerance and iterations.
     imag_tol(ii)     = gmres_tol
     imag_iters(ii)   = gmres_diffusion_iterations

     ! Step 2.3.2: Reset the GMRES Tolerance for the real solve.
     gmres_tol                  = gmres_tol_viscous
     gmres_diffusion_iterations = min( gmres_maxit_viscous, nsg - 1 )
     gmres_restart              = gmres_diffusion_iterations

     ! Step 2.4: Store the vector of unknowns, f, for this wavenumber.
     Fphi(:, ii) = cmplx( real_buff, imag_buff, kind=dp )

     ! Step 2.5.1: Check the solution error for the real component
     call apply_smpm_3D_diffusion_real_valued( rhs_real, real_buff )
     real_error = maxval( abs( rhs_real - real( Frhs(:, ii), kind=dp ) ) )

     ! Step 2.5.1: Check the solution error for the imaginary component.
     call apply_smpm_3D_diffusion_real_valued( rhs_imag, imag_buff )
     imag_error = maxval( abs( rhs_imag - aimag( Frhs(:, ii) ) ) )

     ! Step 2.6: Clear the buffers for next wavenumber.
     real_buff = 0.0_dp
     imag_buff = 0.0_dp
     rhs_real  = 0.0_dp
     rhs_imag  = 0.0_dp

     ! Step 2.6: Display the statistics.
     write( caststr, '(A,I0,A,I0,A,D17.10,A,D17.10)' ) &
          'Real: Wave#', ii, &
          'GMRES Iters', real_iters(ii), &
          'Ext CVG', real_tol(ii), &
          'Error', real_error
     call notify( caststr )
     write( caststr, '(A,I0,A,I0,A,D17.10,A,D17.10)' ) &
          'Imag: Wave#', ii, &
          'GMRES Iters', imag_iters(ii), &
          'Ext CVG', imag_tol(ii), &
          'Error', imag_error
     call notify( caststr )

  enddo

  ! Step 3: Inverse transform the solution.
  call apply_ifft( Fphi, phi )

end subroutine solve_3D_diffusion_GAR

subroutine apply_smpm_3D_diffusion_real_valued( Ax, x )
! Applies the diffusive operator on x for a single wave number, as specified via
! the validation module's wave variable.
!
! NOTE: This is derived from apply_smpm_3D_diffusion().

  use constants, only:               nu_d
  use options, only:                 bc_flag_diffusion, facrobin
  use parallel_linear_algebra, only: pdot_product, pmaxval, pnorm2, znorm2
  use precision, only:               dp
  use timestepping, only:            dt, g0
  use validation, only:              wave
  use woodbury_matrices, only:       rpk

  implicit none

  real(kind=dp), dimension(1:rpk), intent(out) :: Ax
  real(kind=dp), dimension(1:rpk), intent(in)  :: x

  real(kind=dp)                                :: delta
  real(kind=dp), dimension(1:rpk)              :: rho_x, rho_z

  ! Set the viscosity to the value specified by the input file.
  delta = nu_d * dt / g0

  ! Compute the Laplacian in x and z.
  call compute_gradient_and_laplacian( x, Ax, rho_x, rho_z )

  ! Apply the shifted Laplacian.
  Ax = delta * Ax - (delta * x * real( wave, kind=dp )**2.0_dp) - x

  ! Apply the boundary condition.
  call apply_bc( x, Ax, bc_flag_diffusion, delta, rho_x, rho_z )

  ! Apply the patching condition.
  call apply_patching( x, Ax, delta, facrobin, rho_x, rho_z)

end subroutine apply_smpm_3D_diffusion_real_valued
