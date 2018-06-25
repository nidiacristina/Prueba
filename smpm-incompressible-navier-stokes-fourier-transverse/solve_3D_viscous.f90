subroutine solve_3D_viscous( ux_int, uy_int, uz_int )
! Solves the 3D viscous equations with a 2D SMPM discretization and a 1D
! Fourier discretization.

! Notes:
!      - GAR(1/9/16): Added the error computation in the function to avoid
!                     having to recompute in main.f90.

  use constants, only:         nsg, nsuby, nky, nu
  use errors, only:            gmres_viscous_iterations_real, &
                               gmres_viscous_iterations_imag, &
                               linf_viscous_error
  use options, only:           bc_flag_viscous_x, bc_flag_viscous_y, bc_flag_viscous_z, check_numerical_error, &
                               gmres_maxit_viscous, gmres_tol_viscous
  use precision, only:         dp
  use woodbury_matrices, only: rpk
  use parallel_linear_algebra
  use field_variables, only:   ux_b, uy_b, uz_b
  use timestepping, only:      dt, g0, logflag

  implicit none

  ! Set the input/output variables.
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(inout) :: ux_int, uy_int, uz_int

  ! Arrays for storing fourier transformed velocities in complex space
  complex(kind=dp), dimension(1:rpk, 1:nky)               :: Fux, Fuy, Fuz
  complex(kind=dp), dimension(1:3*rpk, 1:nky)             :: Fq, Fq_int

  real(kind=dp), dimension(1:3*rpk, 1:nky)                :: q_real, q_imag
  real(kind=dp), dimension(1:3*rpk, 1:nky)                :: q_int_real, q_int_imag
  real(kind=dp), dimension(1:3*rpk, 1:nky)                :: q_check_real, q_check_imag
  real(kind=dp), dimension(1:3*rpk, 1:nky)                :: q_test

  real(kind=dp)                                             :: linf_viscous_error_real, &
                                                               linf_viscous_error_imag

  ! Other variables.
  real(kind=dp)                                           :: gmres_tol
  integer                                                 :: gmres_restart
  character(len=32)                                       :: caststr
  integer                                                 :: ii

  ! Declare the function used with GMRES to solve the capacitance problem.
  external                                                :: apply_3D_viscous_fourier

  ! Initialize the viscous iteration scalar
  gmres_viscous_iterations_real = 0
  gmres_viscous_iterations_imag = 0

  ! Set some constants.
  gmres_tol                     = gmres_tol_viscous
  gmres_viscous_iterations_real = min( gmres_maxit_viscous, nsg - 1 )
  gmres_restart                 = gmres_viscous_iterations_real

  ! Set the viscous boundary conditions on the x and z boundaries.
  call apply_viscous_bc( ux_int, ux_b, bc_flag_viscous_x, nu * dt / g0 )
  call apply_viscous_bc( uy_int, uy_b, bc_flag_viscous_y, nu * dt / g0 )
  call apply_viscous_bc( uz_int, uz_b, bc_flag_viscous_z, nu * dt / g0 )

  ! Fourier transform the input velocities.
  call apply_fft( ux_int, Fux )
  if (nsuby > 1 ) then
     call apply_fft( uy_int, Fuy )
  else
     Fuy = 0.0_dp
  endif
  call apply_fft( uz_int, Fuz )

  ! Set up an array containing all three velocities.
  Fq_int(1:rpk, :)         = Fux
  Fq_int(rpk+1:2*rpk, :)   = Fuy
  Fq_int(2*rpk+1:3*rpk, :) = Fuz

  ! Decompose into real and imaginary components
  do ii = 1,nky
     q_int_real( :, ii ) = real(  Fq_int( :, ii ), kind=dp )
     q_int_imag( :, ii ) = aimag( Fq_int( :, ii ) )
  enddo

  ! Store velocities for the error calculation after GMRES solve.
  if ( check_numerical_error) then
     q_check_real = q_int_real
     q_check_imag = q_int_imag
  endif

  ! Initialize the solution
  q_real = q_int_real * (1 - exp( -nu * dt ))
  q_imag = 0.0_dp

  ! XXX: this could be folded into the real/imaginary decomposition above if the
  !      check for numerical error isn't affected by it.
  q_int_real = q_int_real / (-g0)
  q_int_imag = q_int_imag / (-g0)

  ! Solve the real part 
  call compute_gmres_householder( q_real, q_int_real, 3 * rpk * nky, &
                                  gmres_tol, gmres_viscous_iterations_real, gmres_restart, &
                                  apply_3D_viscous_fourier )

  write( caststr, '(I10) '   ) gmres_viscous_iterations_real
  call notify_cond( logflag, '         Viscous solve GMRES iterations (real): ' // caststr )
  write( caststr, '(D17.10)') gmres_tol
  call notify_cond( logflag, '         Viscous solve GMRES rel. error (real): ' // caststr )

 
 ! Check the error in GMRES solve.
  if ( check_numerical_error) then

     ! Reset Buffers
     q_test = 0.0_dp

     ! Check real part
     call apply_3D_viscous_fourier( q_test, q_real )

     linf_viscous_error_real = pmaxval( abs(  reshape( q_test + q_check_real / g0 , &
                                                       (/3 * nky * rpk/) ) ) )

     write( caststr, '(D17.10)' ) linf_viscous_error_real
     call notify_cond( logflag, '         Linf error in viscous solve    (real): ' // caststr )

  endif

  ! Solve the imaginary part.
  ! Check if there is a transverse domain first.
  if ( nsuby > 1 ) then

     ! Initialize the solution
     q_imag = q_int_imag * (1 - exp( -nu * dt ))

     ! Reset some constants.
     gmres_tol                     = gmres_tol_viscous
     gmres_viscous_iterations_imag = min( gmres_maxit_viscous, nsg - 1 )
     gmres_restart                 = gmres_viscous_iterations_imag

     ! Solve the imaginary part
     call compute_gmres_householder( q_imag, q_int_imag, 3 * rpk * nky, &
                                     gmres_tol, gmres_viscous_iterations_imag, gmres_restart, &
                                     apply_3D_viscous_fourier )

     call notify_cond( logflag, ' ' )
     write( caststr, '(I10) '   ) gmres_viscous_iterations_imag
     call notify_cond( logflag, '         Viscous solve GMRES iterations (imag): ' // caststr )
     write( caststr, '(D17.10)') gmres_tol
     call notify_cond( logflag, '         Viscous solve GMRES rel. error (imag): ' // caststr )

     ! Check the error in GMRES solve.
     if ( check_numerical_error) then

        ! Reset Buffers
        q_test = 0.0_dp

        call apply_3D_viscous_fourier( q_test, q_imag )

        linf_viscous_error_imag = pmaxval( abs(  reshape( q_test + q_check_imag / g0 , &
                                                          (/ 3 * nky * rpk/) ) ) )

        write( caststr, '(D17.10)' ) linf_viscous_error_imag
        call notify_cond( logflag, '         Linf error in viscous solve    (imag): ' // caststr )

     endif

  else ! Skip the imag solve associated with the transverse direction

     gmres_viscous_iterations_imag = 0
     gmres_tol                     = 0.0_dp
     linf_viscous_error_imag       = 0.0_dp

  endif


  call notify_cond( logflag, ' ' )

  ! Reassemble solution and error into complex form.
  do ii = 1,nky
     Fq( :, ii ) = cmplx(  q_real( :,ii ) , q_imag( :, ii ), kind=dp )
  enddo
  linf_viscous_error = cmplx( linf_viscous_error_real, &
                              linf_viscous_error_imag, &
                              kind=dp )

  ! Unpack the three velocities from the packed array.
  Fux = Fq(1:rpk, :)
  Fuy = Fq(rpk+1:2*rpk, :)
  Fuz = Fq(2*rpk +1:3*rpk, :)

  ! Inverse Fourier transform the solution.
  call apply_ifft( Fux, ux_int )
  if ( nsuby > 1 ) then
     call apply_ifft( Fuy, uy_int )
  else
     uy_int = 0.0_dp
  endif
  call apply_ifft( Fuz, uz_int )


end subroutine solve_3D_viscous
