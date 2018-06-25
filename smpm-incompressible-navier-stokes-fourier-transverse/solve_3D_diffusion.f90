subroutine solve_3D_diffusion( rho_int )
! This routine solves the 3D Helmholtz Eqn, which is used for the diffusion step.
! The Helmholtz Equation is stated as:
!
!   alpha * Lap f - f = phi
!

  use constants, only:         nsg, nsuby, nky, nu_d
  use errors, only:            gmres_diffusion_iterations_real, gmres_diffusion_iterations_imag, &
                               linf_diffusion_error
  use options, only:           check_numerical_error, gmres_maxit_viscous, gmres_tol_viscous
  use parallel_linear_algebra
  use precision, only:         dp
  use woodbury_matrices, only: rpk
  use timestepping, only:      dt, g0, logflag

  implicit none

  ! Set the input/output variables.
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(inout)   :: rho_int

  ! Buffers for error check
  real(kind=dp), dimension(1:rpk,1:nky)                     :: rho_check_real, rho_check_imag, rho_test

  ! Buffers for RHSs
  real(kind=dp), dimension(1:rpk,1:nky)                     :: rho_int_real, rho_int_imag

  ! Buffers for solution
  real(kind=dp), dimension(1:rpk,1:nky)                     :: rho_real, rho_imag

  ! Buffers for complex space
  complex(kind=dp), dimension(1:rpk, 1:nky)                 :: Frho_int, Frho

  real(kind=dp)                                             :: linf_diffusion_error_real, &
                                                               linf_diffusion_error_imag

  ! Other variables.
  integer                                                   :: gmres_restart
  real(kind=dp)                                             :: gmres_tol
  character(len=32)                                         :: caststr
  integer                                                   :: ii 

  ! Declare the function used with GMRES to solve the capacitance problem
  ! (i.e. Linear Operator).
  external                                                  :: apply_3D_diffusive_fourier

  ! Initialize Iterations
  gmres_diffusion_iterations_real = 0
  gmres_diffusion_iterations_imag = 0

  ! Initialize complex arrays.
  Frho       = cmplx( 0.0_dp, 0.0_dp, kind=dp )
  Frho_int   = cmplx( 0.0_dp, 0.0_dp, kind=dp )

  ! Set GMRES tolerance.
  gmres_tol                       = gmres_tol_viscous
  gmres_diffusion_iterations_real = min( gmres_maxit_viscous, nsg - 1 )
  gmres_restart                   = gmres_diffusion_iterations_real

  ! Apply boundary conditions on the x and z boundary to the density
  ! Note: Eventhough we are not applying bc's to the diffusion equation because
  !       technically they are 0 for our current problems, we have to leave 
  !       the environment ready for future problems. 

  ! Fourier transform the density field
  call apply_fft( rho_int, Frho_int )

  ! Decompose into real and imaginary components
  do ii = 1,nky
     rho_int_real( : , ii ) = real( Frho_int( : , ii ), kind=dp )
     rho_int_imag( : , ii ) = aimag( Frho_int( : , ii ) )
  enddo

  ! Store Solution for checking error
  if ( check_numerical_error) then
     rho_check_real = rho_int_real
     rho_check_imag = rho_int_imag
  endif
 
  ! Initialize the solution
  rho_real = rho_int_real * (1 - exp( -nu_d * dt ))
  rho_imag = 0.0_dp

  ! XXX: this could be folded into the real/imaginary decomposition above if the
  !      check for numerical error isn't affected by it.
  rho_int_real = rho_int_real / (-g0)
  rho_int_imag = rho_int_imag / (-g0)

  ! Solve for the real part.
  call compute_gmres_householder( rho_real, rho_int_real,  rpk * nky, &
                                  gmres_tol, gmres_diffusion_iterations_real, gmres_restart, &
                                  apply_3D_diffusive_fourier )
              
  write( caststr, '(I10) '   ) gmres_diffusion_iterations_real
  call notify_cond( logflag, '         Diffusive solve GMRES iters.   (real): ' // caststr )
  write( caststr, '(D17.10)') gmres_tol
  call notify_cond( logflag, '         Diffusive solve GMRES rel. err (real): ' // caststr )

 ! Check the error in GMRES solve.
  if ( check_numerical_error) then

     ! Reset Buffers
     rho_test = 0.0_dp

     ! Check real part
     call apply_3D_diffusive_fourier( rho_test, rho_real )

     linf_diffusion_error_real = pmaxval( abs(  reshape( rho_test + rho_check_real / g0, &
                                                         (/ nky * rpk/) ) ) )

     write( caststr, '(D17.10)' ) linf_diffusion_error_real
     call notify_cond( logflag, '         Linf error in diffusive solve  (real): ' // caststr )
  else
     linf_diffusion_error_real = 0.0_dp
  endif

  
  call notify_cond( logflag, ' ' )

  ! Solve for the imaginary part.
  ! Check if there is transverse domain.
  if ( nsuby > 1 ) then

     ! Initialize the solution
     rho_imag = rho_int_imag * (1 - exp( -nu_d * dt ))

     ! Reset GMRES tolerance.
     gmres_tol                       = gmres_tol_viscous
     gmres_diffusion_iterations_imag = min( gmres_maxit_viscous, nsg - 1 )
     gmres_restart                   = gmres_diffusion_iterations_imag

     ! Solve for the real part
     call compute_gmres_householder( rho_imag, rho_int_imag,  rpk * nky, &
                                     gmres_tol, gmres_diffusion_iterations_imag, gmres_restart, &
                                     apply_3D_diffusive_fourier )
              
     write( caststr, '(I10) '   ) gmres_diffusion_iterations_imag
     call notify_cond( logflag, '         Diffusive solve GMRES iters.   (imag): ' // caststr )
     write( caststr, '(D17.10)') gmres_tol
     call notify_cond( logflag, '         Diffusive solve GMRES rel. err (imag): ' // caststr )

     ! Check the error in GMRES solve.
     if ( check_numerical_error) then

        ! Reset Buffers
        rho_test = 0.0_dp

        ! Check real part
        call apply_3D_diffusive_fourier( rho_test, rho_imag )

        linf_diffusion_error_imag = pmaxval( abs(  reshape( rho_test + rho_check_imag / g0 , &
                                                            (/ nky * rpk/) ) ) )

        write( caststr, '(D17.10)' )  linf_diffusion_error_imag
        call notify_cond( logflag, '         Linf error in diffusive solve  (imag): ' // caststr )
     else
        linf_diffusion_error_imag = 0.0_dp
     endif

  else ! Skip the imag solve associated with the transverse direction

     gmres_diffusion_iterations_imag = 0
     gmres_tol                       = 0.0_dp
     linf_diffusion_error_imag       = 0.0_dp

  endif

  call notify_cond( logflag, ' ' )

  ! Reassemble solution and error into complex form.
  do ii = 1, nky
     Frho( : , ii) = cmplx( rho_real( :, ii), rho_imag( :, ii) , kind=dp  )
  enddo
  linf_diffusion_error = cmplx( linf_diffusion_error_real, &
                                linf_diffusion_error_imag, &
                                kind=dp )

  ! Inverse Transform Solution into physica space
  call apply_ifft( Frho, rho_int )


end subroutine solve_3D_diffusion
