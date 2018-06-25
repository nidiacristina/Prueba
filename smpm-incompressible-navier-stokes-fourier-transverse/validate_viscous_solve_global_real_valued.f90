subroutine validate_viscous_solve_global_real_valued( output_prefix, tolerance, success_flag )
! Validates the viscous solver using XXX?

  use constants, only:               nsuby, nky, nu
  use field_variables, only:         ux0, uy0, uz0, ux, uy, uz
  use geom, only:                    cx, cz
  use parallel_linear_algebra, only: pdot_product, pmaxval, pnorm2, znorm2
  use precision, only:               dp
  use timestepping, only:            g0, dt
  use transverse, only:              cy, qy
  use woodbury_matrices, only:       rpk
  use validation, only:              wave

  implicit none

  character(len=*), intent(in)                   :: output_prefix
  real(kind=dp), intent(in)                      :: tolerance
  logical, intent(out)                           :: success_flag

  ! General use variables.
  real(kind=dp)                                  :: pi, mm, ll, kk
  integer                                        :: ii

  real(kind=dp)                                  :: linf_ux, linf_uy, linf_uz, &
                                                    linf_ux_single, linf_uy_single, linf_uz_single

  ! Define variables for testing real/complex GMRES implementation.
  real(kind=dp), allocatable, dimension(:, :)    :: rhs_ux, rhs_uy, rhs_uz, dummy1, dummy2, dummy3
  complex(kind=dp), allocatable, dimension(:, :) :: Fq_int, Frhs, Fux_rhs, Fuy_rhs, Fuz_rhs
  complex(kind=dp), allocatable, dimension(:, :) :: Fux, Fuy, Fuz

  real(kind=dp), allocatable, dimension(:)       :: real_buff, imag_buff

  character(len=256)                             :: caststr

  external                                       :: apply_3D_viscous_fourier_real_valued

  allocate( Fux(1:rpk, 1:nky), Fuy(1:rpk, 1:nky), Fuz(1:rpk, 1:nky) )
  allocate( Fux_rhs(1:rpk, 1:nky), Fuy_rhs(1:rpk, 1:nky), Fuz_rhs(1:rpk, 1:nky) )

  allocate( Fq_int(1:3*rpk, 1:nky), Frhs(1:3*rpk, 1:nky) )
  allocate( rhs_ux(1:rpk, 1:nsuby), rhs_uy(1:rpk, 1:nsuby), rhs_uz(1:rpk, 1:nsuby) )

  allocate( real_buff(1:3*rpk), imag_buff(1:3*rpk) )

  ! Dummy Arrays to Store Analytical Solution
  allocate(dummy1(1:rpk, 1:nsuby), dummy2(1:rpk, 1:nsuby), dummy3(1:rpk, 1:nsuby) )

  ! Get pi.
  pi = 2.0_dp * acos( 0.0_dp )

  ! Set the timestepping coefficient.
  g0 = 1.0_dp

  ! Set wavenumber.
  kk = 8.0_dp * pi
  ll = 8.0_dp * pi
  mm = 8.0_dp * pi

  ! Dirichlet
  ux0 = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
  uy0 = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
  uz0 = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )

  ! Initialize arrays
  dummy1 = 0.0_dp
  dummy2 = 0.0_dp
  dummy3 = 0.0_dp

  Fux = 0.0_dp
  Fuy = 0.0_dp
  Fuz = 0.0_dp

  Fq_int = 0.0_dp
  Frhs   = 0.0_dp

  Fux_rhs = 0.0_dp
  Fuy_rhs = 0.0_dp
  Fuz_rhs = 0.0_dp

  rhs_ux = 0.0_dp
  rhs_uz = 0.0_dp
  rhs_uy = 0.0_dp

  call notify( 'Solving via a global vector for all velocities with real-valued GMRES.' )
  ! Construct the RHS

  ! Fourier transform the input velocities.
  call apply_fft( ux0, Fux )
  call apply_fft( uy0, Fuy )
  call apply_fft( uz0, Fuz )

  ! Set up an array containing all three velocities.
  Fq_int(1:rpk, :)         = Fux
  Fq_int(rpk+1:2*rpk, :)   = Fuy
  Fq_int(2*rpk+1:3*rpk, :) = Fuz

  do ii = 1, nky
     ! Set the wave number for the linear operators below.
     wave = qy(ii)

     ! Apply the linear operator to real component.
     call apply_3D_viscous_fourier_real_valued( real_buff, real( Fq_int(:, ii), kind=dp ) )
     call apply_3D_viscous_fourier_real_valued( imag_buff, aimag( Fq_int(:, ii) ) )

     ! Store the complex RHS.
     Frhs(:, ii) = cmplx( real_buff, imag_buff, kind=dp )

  enddo

  ! Set up an array containing all three velocities.
  Fux_rhs = Frhs(1:rpk, :)
  Fuy_rhs = Frhs(rpk+1:2*rpk, :)
  Fuz_rhs = Frhs(2*rpk+1:3*rpk, :)

  call apply_ifft( Fux_rhs, rhs_ux )
  call apply_ifft( Fuy_rhs, rhs_uy )
  call apply_ifft( Fuz_rhs, rhs_uz )

  ! Check the linear operator error.
  dummy1 = -(nu * dt / g0) * (kk**2.0_dp + ll**2.0_dp + mm**2.0_dp) * sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) &
           -(sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) )

  dummy2 = -(nu * dt / g0) * (kk**2.0_dp + ll**2.0_dp + mm**2.0_dp) * sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) &
           -(sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) )

  dummy3 = -(nu * dt / g0) * (kk**2.0_dp + ll**2.0_dp + mm**2.0_dp) * sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) &
           -(sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) )

  ! Report errors.
  linf_ux = maxval( abs( dummy1 - ux ) )
  linf_uy = maxval( abs( dummy2 - uy ) )
  linf_uz = maxval( abs( dummy3 - uz ) )

  write( caststr, '(A,D17.10,A,D17.10,A,D17.10)' ) &
       'Linf error of Linear Operator: ux', linf_ux, &
       'Linf error of Linear Operator: uy', linf_uy, &
       'Linf error of Linear Operator: uy', linf_uz
  call notify( caststr )

  ! Step 3.1: viscous diffusion solve.
  call solve_3D_viscous_real_valued( rhs_ux, rhs_uy, rhs_uz )

  dummy1 = 0.0_dp
  dummy2 = 0.0_dp
  dummy3 = 0.0_dp

  ! Dirichlet boundary conditions.
  dummy1 =  sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
  dummy2 =  sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
  dummy3 =  sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )

  ! Neumann boundary conditions.
  !dummy1 =  cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz )
  !dummy2 =  cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz )
  !dummy3 =  cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz )

  ux = -rhs_ux
  uy = -rhs_uy
  uz = -rhs_uz

  linf_ux_single = maxval( abs( dummy1 - ux ) )
  linf_uy_single = maxval( abs( dummy2 - uy ) )
  linf_uz_single = maxval( abs( dummy3 - uz ) )

  write( caststr, '(A,D17.10,D17.10,D17.10)' ) 'Error in single vector (ux,uy,uz):', &
       linf_ux_single, linf_uy_single, linf_uz_single
  call notify( caststr )

  ! Write the solution to disk.
  open( 65, file=trim( output_prefix ) // '_viscous_solution_single_vector.txt' )
  write( 65, * ) ux
  write( 65, * ) uy
  write( 65, * ) uz
  close( 65 )

   success_flag = ((linf_ux        < tolerance) .and. &
                   (linf_uy        < tolerance) .and. &
                   (linf_uz        < tolerance) .and. &
                   (linf_ux_single < tolerance) .and. &
                   (linf_uy_single < tolerance) .and. &
                   (linf_uz_single < tolerance))

end subroutine validate_viscous_solve_global_real_valued

subroutine solve_3D_viscous_real_valued( ux_int, uy_int, uz_int )
! Solves the 3D viscous equations with a 2D SMPM discretization and a 1D Fourier
! discretization.
!
! NOTE: Derived from solve_3D_viscous().  XXX fill in what has changed

  use constants, only:         nsg, nsuby, nky, nu
  use options, only:           bc_flag_viscous_x, bc_flag_viscous_y, bc_flag_viscous_z, &
                               gmres_maxit_viscous, gmres_tol_viscous
  use woodbury_matrices, only: rpk
  use parallel_linear_algebra
  use field_variables, only:   ux_b, uy_b, uz_b
  use timestepping, only:      dt, g0
  use transverse, only:        qy
  use validation, only:        wave

  implicit none

  ! Set the input/output variables.
  real(kind=dp), dimension( 1:rpk, 1:nsuby ), intent(inout) :: ux_int, uy_int, uz_int
  real(kind=dp), dimension(1:nky)                           :: real_tol, imag_tol
  integer, dimension(1:nky)                                 :: real_iters, imag_iters

  ! Arrays for storing fourier transformed velocities.
  complex(kind=dp), dimension( 1:rpk, 1:nky )               :: Fux, Fuy, Fuz
  complex(kind=dp), dimension( 1:3*rpk, 1:nky )             :: Fq, Fq_int
  real(kind=dp), dimension(1:3*rpk)                         :: real_buff, imag_buff
  real(kind=dp), dimension(1:3*rpk)                         :: rhs_real, rhs_imag
  real(kind=dp)                                             :: imag_error, real_error

  ! Other variables.
  real(kind=dp)                                             :: gmres_tol
  integer                                                   :: gmres_viscous_iterations
  integer                                                   :: gmres_restart

  integer                                                   :: ii

  character(len=256)                                        :: caststr

  ! Declare the function used with GMRES to solve the capacitance problem.
  external                                                  :: apply_3D_viscous_fourier_real_valued

  ! Set some constants.
  gmres_tol                = gmres_tol_viscous
  gmres_viscous_iterations = min( gmres_maxit_viscous, nsg - 1 )
  gmres_restart            = gmres_viscous_iterations

  ! Set the viscous boundary conditions on the x and z boundaries.
  call apply_viscous_bc( ux_int, ux_b, bc_flag_viscous_x, nu * dt / g0 )
  call apply_viscous_bc( uy_int, uy_b, bc_flag_viscous_y, nu * dt / g0 )
  call apply_viscous_bc( uz_int, uz_b, bc_flag_viscous_z, nu * dt / g0 )

  ! Fourier transform the input velocities.
  call apply_fft( ux_int, Fux )
  call apply_fft( uy_int, Fuy )
  call apply_fft( uz_int, Fuz )

  ! Set up an array containing all three velocities.
  Fq_int(1:rpk, :)         = Fux
  Fq_int(rpk+1:2*rpk, :)   = Fuy
  Fq_int(2*rpk+1:3*rpk, :) = Fuz

  ! Apply the time-stepping coefficient.
  Fq_int = -Fq_int / g0

  ! Set the initial guess to zero.
  Fq = Fq_int * (1 - exp( -nu * dt ))

  gmres_tol_viscous = 10**(-6)

  ! Iterate through all of the wavenumbers.
  do ii = 1, nky

     ! Set the wavenumber for the linear solves below.
     wave = qy(ii)

     !Initialize the real buffer.
     real_buff = real( Fq(:, ii), kind=dp )

     ! Step 2.2: Solve the real linear system via GMRES.
     call compute_gmres_householder( real_buff, real( Fq_int(:, ii), kind=dp ), 3 * rpk, &
                                     gmres_tol, gmres_viscous_iterations, gmres_restart, &
                                     apply_3D_viscous_fourier_real_valued )

     ! Step 2.2.1: Store global tolerance and iterations for the real solve.
     real_tol(ii)     = gmres_tol
     real_iters(ii)   = int(gmres_viscous_iterations)

     ! Step 2.2.2: Reset the GMRES tolerance for the imaginary solve.
     gmres_tol                = gmres_tol_viscous
     gmres_viscous_iterations = min( gmres_maxit_viscous, nsg - 1 )
     gmres_restart            = gmres_viscous_iterations

     ! Initialize the imaginary buffer.
     imag_buff = aimag( Fq(:, ii) )

     ! Step 2.3: Solve complex linear system via GMRES.
     call compute_gmres_householder( imag_buff, aimag( Fq_int(:, ii) ), 3 * rpk, &
                                           gmres_tol, gmres_viscous_iterations, gmres_restart, &
                                           apply_3D_viscous_fourier_real_valued )

     ! Step 2.3.2: Reset the GMRES tolerance for real solve.
     gmres_tol                  = gmres_tol_viscous
     gmres_viscous_iterations = min( gmres_maxit_viscous, nsg - 1 )
     gmres_restart              = gmres_viscous_iterations

     ! Step 2.4: Store the vector of unknowns, f, for this wavenumber.
     Fq(:, ii) = cmplx( real_buff, imag_buff, kind=dp )

     ! Step 2.5.1: Check the solution error in the real component.
     call apply_3D_viscous_fourier_real_valued(rhs_real ,real_buff )
     real_error = maxval( abs( rhs_real - real( Fq_int(:, ii), kind=dp ) ) )

     ! Step 2.5.1: Check the solution error in the imaginary component.
     call apply_3D_viscous_fourier_real_valued(rhs_imag, imag_buff )
     imag_error = maxval( abs( rhs_imag - aimag( Fq_int(:, ii) ) ) )

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

  ! Unpack the three velocities from the packed array.
  Fux = Fq(1:rpk, :)
  Fuy = Fq(rpk+1:2*rpk, :)
  Fuz = Fq(2*rpk +1:3*rpk, :)

  ! Inverse Fourier transform the solution.
  call apply_ifft( Fux, ux_int )
  call apply_ifft( Fuy, uy_int )
  call apply_ifft( Fuz, uz_int )

end subroutine solve_3D_viscous_real_valued

subroutine apply_3D_viscous_fourier_real_valued( Aq, q )
! Applies the Poisson-Neumann operator on q.

  use constants,only:          nu
  use options, only:           bc_flag_viscous, facrobin
  use precision, only:         dp
  use timestepping, only:      dt, g0
  use validation, only:        wave
  use woodbury_matrices, only: rpk

  implicit none

  real(kind=dp), dimension(1:3*rpk), intent(out) :: Aq
  real(kind=dp), dimension(1:3*rpk), intent(in)  :: q

  real(kind=dp), dimension(1:2*rpk)              :: qx_and_qz, Aq_ux_and_uz
  real(kind=dp), dimension(1:3*rpk)              :: q_x, q_z
  real(kind=dp), parameter                       :: delta_poisson = 1.0_dp
  integer                                        :: jj, a, z
  real(kind=dp), dimension(1:2*rpk)              :: dqdx, dqdz
  real(kind=dp)                                  :: delta

  ! Set the viscosity to the value specified by the input file.
  delta = nu * dt / g0 ! gamma0 because of KIO time-splitting

  ! Loop over the three components of velocity.
  do jj = 1, 3

     ! Set array bounds for this component.
     a = (jj - 1) * rpk + 1
     z = (jj - 1) * rpk + rpk

     ! Compute the Laplacian in x and z.
     call compute_gradient_and_laplacian( q(a:z) , Aq(a:z), q_x(a:z), q_z(a:z) )

     ! Apply the shifted Laplacian to the real and imaginary parts separately.
     Aq(a:z) = delta * Aq(a:z) - (delta * q(a:z) * real( wave, kind=dp )**2.0_dp ) - q(a:z)

     ! Apply the patching conditions.
     call apply_patching( q(a:z), Aq(a:z), delta, facrobin, q_x(a:z), q_z(a:z) )

  enddo

  ! x derivative of ux and uz, not uy.
  dqdx(1:rpk)       = q_x(1:rpk)
  dqdx(rpk+1:2*rpk) = q_x(2*rpk+1:3*rpk)

  ! z derivative of ux and uz, not uy.
  dqdz(1:rpk)       = q_z(1:rpk)
  dqdz(rpk+1:2*rpk) = q_z(2*rpk+1:3*rpk)

  ! Set up an array that contains ux and uz, not uy.
  qx_and_qz(1:rpk)       = q(1:rpk)
  qx_and_qz(rpk+1:2*rpk) = q(2*rpk+1:3*rpk)

  ! Set up an array to hold the viscous operator applied to the x and z
  ! components.
  Aq_ux_and_uz(1:rpk)       = Aq( 1:rpk )
  Aq_ux_and_uz(rpk+1:2*rpk) = Aq(2*rpk+1:3*rpk)

  ! Apply the vector viscous boundary conditions in (ux,uz).
  call apply_vector_viscous_bc( qx_and_qz, Aq_ux_and_uz, dqdx, dqdz, &
                                bc_flag_viscous, delta )

  ! Store the result in the output vector in the right place.
  Aq(1:rpk)         = Aq_ux_and_uz(1:rpk)
  Aq(2*rpk+1:3*rpk) = Aq_ux_and_uz(rpk+1:2*rpk)

end subroutine apply_3D_viscous_fourier_real_valued
