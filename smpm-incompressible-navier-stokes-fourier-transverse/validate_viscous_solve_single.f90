subroutine validate_viscous_solve_single( output_prefix, tolerance, success_flag )
! Validates the viscous solver using XXX?

  use constants, only:               nsg, nsuby, nky, nu
  use field_variables, only:         ux0, uy0, uz0, ux, uy, uz, ux_b, uy_b, uz_b
  use geom, only:                    cx, cz
  use options, only:                 gmres_maxit_viscous, gmres_tol_viscous, &
                                     bc_flag_viscous_x, bc_flag_viscous_z
  use parallel_linear_algebra, only: pdot_product, pmaxval, pnorm2, znorm2
  use precision, only:               dp
  use timestepping, only:            g0, dt
  use transverse, only:              cy
  use woodbury_matrices, only:       rpk

  implicit none

  character(len=*), intent(in)                  :: output_prefix
  real(kind=dp), intent(in)                     :: tolerance
  logical, intent(out)                          :: success_flag

  ! General use variables.
  real(kind=dp)                                 :: pi, mm, ll, kk

  real(kind=dp)                                 :: linf_ux, linf_uy, linf_uz, &
                                                   linf_ux_multiple, linf_uy_multiple, linf_uz_multiple

  ! Variables for GMRES.
  real(kind=dp)                                 :: gmres_tol
  integer                                       :: gmres_viscous_iterations

  ! Define variables for testing real/complex GMRES implementation.
  real(kind=dp), allocatable, dimension(:,:)    :: rhs_ux, rhs_uy, rhs_uz, dummy1, dummy2, dummy3
  complex(kind=dp), allocatable, dimension(:,:) :: Fq_int, Fqrhs, Fux_rhs, Fuy_rhs, Fuz_rhs
  complex(kind=dp), allocatable, dimension(:,:) :: Fux, Fuy, Fuz

  character(len=256)                            :: caststr

  external                                      :: apply_3D_viscous_fourier_single_velocity

  allocate( Fux(1:rpk, 1:nky), Fuy(1:rpk, 1:nky), Fuz(1:rpk, 1:nky) )
  allocate( Fux_rhs(1:rpk, 1:nky), Fuy_rhs(1:rpk, 1:nky), Fuz_rhs(1:rpk, 1:nky) )

  allocate( Fq_int(1:3*rpk, 1:nky), Fqrhs(1:3*rpk, 1:nky) )
  allocate( rhs_ux(1:rpk, 1:nsuby), rhs_uy(1:rpk, 1:nsuby), rhs_uz(1:rpk, 1:nsuby) )

  ! Dummy arrays to store the analytical solution.
  allocate( dummy1(1:rpk, 1:nsuby), dummy2(1:rpk, 1:nsuby), dummy3(1:rpk, 1:nsuby) )

  ! Get pi.
  pi = 2.0_dp * acos( 0.0_dp )

  ! Set the timestepping coefficient.
  g0 = 1.0_dp

  ! Set the wavenumber.
  kk = 8.0_dp * pi
  ll = 8.0_dp * pi
  mm = 8.0_dp * pi

  ! Dirichlet boundary conditions.
  ux0 = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
  uy0 = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
  uz0 = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )

  ! Initialize arrays.
  dummy1 = 0.0_dp
  dummy2 = 0.0_dp
  dummy3 = 0.0_dp

  Fux = 0.0_dp
  Fuy = 0.0_dp
  Fuz = 0.0_dp

  Fq_int = 0.0_dp
  Fqrhs  = 0.0_dp

  Fux_rhs = 0.0_dp
  Fuy_rhs = 0.0_dp
  Fuz_rhs = 0.0_dp

  rhs_ux = 0.0_dp
  rhs_uz = 0.0_dp
  rhs_uy = 0.0_dp

  call notify( 'Solving via single vector per velocity component.' )

  ! Construct the RHS by Fourier transforming the input velocities.
  call apply_fft( ux0, Fux )
  call apply_fft( uy0, Fuy )
  call apply_fft( uz0, Fuz )

  call apply_3D_viscous_fourier_single_velocity( Fux_rhs, Fux )
  call apply_3D_viscous_fourier_single_velocity( Fuy_rhs, Fuy )
  call apply_3D_viscous_fourier_single_velocity( Fuz_rhs, Fuz )

  call apply_ifft( Fux_rhs, rhs_ux )
  call apply_ifft( Fuy_rhs, rhs_uy )
  call apply_ifft( Fuz_rhs, rhs_uz )

   ! Dirichlet boundary conditions.
   dummy1 = -(nu * dt / g0) * (kk**2.0_dp + ll**2.0_dp + mm**2.0_dp) * sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) &
            -(sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) )

   dummy2 = -(nu * dt / g0) * (kk**2.0_dp + ll**2.0_dp + mm**2.0_dp) * sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) &
            -(sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) )

   dummy3 = -(nu * dt / g0) * (kk**2.0_dp + ll**2.0_dp + mm**2.0_dp) * sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) &
            -(sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) )

   ! Neumann boundary conditions.
!   dummy1 = -(nu * dt / g0) * (kk**2.0_dp + ll**2.0_dp + mm**2.0_dp) * cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz ) &
!            -(cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz ) )

!   dummy2 = -(nu * dt / g0) * (kk**2.0_dp + ll**2.0_dp + mm**2.0_dp) * cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz ) &
!            -(cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz ) )

!   dummy3 = -(nu * dt / g0) * (kk**2.0_dp + ll**2.0_dp + mm**2.0_dp) * cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz ) &
!            -(cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz ) )

   ! Report errors.
   linf_ux = maxval( abs( dummy1 - rhs_ux ) )
   linf_uy = maxval( abs( dummy2 - rhs_uy ) )
   linf_uz = maxval( abs( dummy3 - rhs_uz ) )

   write( caststr, '(A,D17.10,D17.10,D17.10)' ) 'Error in multiple velocity vectors linear operator(ux,uy,uz):', &
        linf_ux, linf_uy, linf_uz
   call notify( caststr )

   ! Step 3.1: Viscous diffusion solve.
   gmres_tol                = gmres_tol_viscous
   gmres_viscous_iterations = min( gmres_maxit_viscous, nsg - 1 )

   call solve_3D_viscous_single_velocity( rhs_ux, ux_b, bc_flag_viscous_x, &
                                          gmres_viscous_iterations, gmres_tol )
   write( caststr, '(A,I0,A,D17.10)' ) &
        'Viscous U: GMRES Iters=', gmres_viscous_iterations, 'GMRES Res=', gmres_tol
   call notify( caststr )

   gmres_tol                = gmres_tol_viscous
   gmres_viscous_iterations = min( gmres_maxit_viscous, nsg - 1 )

   call solve_3D_viscous_single_velocity( rhs_uy, uy_b, bc_flag_viscous_z, &
                                          gmres_viscous_iterations, gmres_tol )
   write( caststr, '(A,I0,A,D17.10)' ) &
        'Viscous V: GMRES Iters=', gmres_viscous_iterations, 'GMRES Res=', gmres_tol
   call notify( caststr )

   gmres_tol                = gmres_tol_viscous
   gmres_viscous_iterations = min( gmres_maxit_viscous, nsg - 1 )
   call solve_3D_viscous_single_velocity( rhs_uz, uz_b, bc_flag_viscous_z, &
                                          gmres_viscous_iterations, gmres_tol )
   write( caststr, '(A,I0,A,D17.10)' ) &
        'Viscous W: GMRES Iters=', gmres_viscous_iterations, 'GMRES Res=', gmres_tol
   call notify( caststr )

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

   ! Set the viscous boundary conditions on the x and z boundaries.
   !call apply_viscous_bc( dummy1, ux_b, bc_flag_viscous_x, nu * dt / g0 )
   !call apply_viscous_bc( dummy2, uy_b, bc_flag_viscous_y, nu * dt / g0 )
   !call apply_viscous_bc( dummy3, uz_b, bc_flag_viscous_z, nu * dt / g0 )

   ux = -rhs_ux
   uy = -rhs_uy
   uz = -rhs_uz

   linf_ux_multiple = maxval( abs( dummy1 - ux ) )
   linf_uy_multiple = maxval( abs( dummy2 - uy ) )
   linf_uz_multiple = maxval( abs( dummy3 - uz ) )

   write( caststr, '(A,D17.10,D17.10,D17.10)' ) 'Error in multiple velocities (ux,uy,uz):', &
        linf_ux_multiple, linf_uy_multiple, linf_uz_multiple
   call notify( caststr )

   ! Write Solution to File
   open( 65, file=trim( output_prefix ) // '_viscous_solution.txt' )
   write( 65, * ) ux
   write( 65, * ) uy
   write( 65, * ) uz
   close( 65 )

   success_flag = ((linf_ux          < tolerance) .and. &
                   (linf_uy          < tolerance) .and. &
                   (linf_uz          < tolerance) .and. &
                   (linf_ux_multiple < tolerance) .and. &
                   (linf_uy_multiple < tolerance) .and. &
                   (linf_uz_multiple < tolerance))

 end subroutine validate_viscous_solve_single


subroutine apply_3D_viscous_fourier_single_velocity( Aq, q )
! Applies the Poisson-Neumann operator on q.

  use options, only:           bc_flag_viscous, facrobin
  use woodbury_matrices, only: rpk
  use constants,only:          nky, nu
  use precision, only:         dp
  use transverse,only:         qy
  use timestepping, only:      dt, g0

  implicit none

  complex(kind=dp), dimension(1:rpk, 1:nky), intent(out) :: Aq
  complex(kind=dp), dimension(1:rpk, 1:nky), intent(in)  :: q

  real(kind=dp), dimension(1:rpk)                        :: real_buff, imag_buff, Rq_x, Rq_z, Iq_x, Iq_z
  integer                                                :: ii
  real(kind=dp)                                          :: delta

  ! Set the viscosity to the value specified by the input file.
  delta = nu * dt / g0 ! gamma0 because of KIO time-splitting.

  ! Loop over each wavenumber.
  do ii = 1, nky

     ! Apply the shifted Laplacian to the real and imaginary parts separately.
     call compute_gradient_and_laplacian(  real( q(:, ii), kind=dp ), real_buff, Rq_x, Rq_z )
     call compute_gradient_and_laplacian( aimag( q(:, ii) ),         imag_buff, Iq_x, Iq_z )
     Aq(:, ii) = (delta * cmplx( real_buff, imag_buff, kind=dp ) - (delta * qy(ii)**2 * q(:, ii)) &
                  - q(:, ii))

     ! Set up an array to hold the diffusive operator applied to the x and z
     ! components.
     real_buff       = real(  Aq(:, ii), kind=dp )
     imag_buff       = aimag( Aq(:, ii) )

     ! Apply the patching conditions.
     call apply_patching( real( q(:, ii), kind=dp ), real_buff, delta, facrobin, Rq_x, Rq_z )
     call apply_patching( aimag( q(:, ii) ),        imag_buff, delta, facrobin, Iq_x, Iq_z )

     ! Apply the diffusive boundary conditions in rho for both real and imaginary components.
     call apply_bc( real( q(:, ii), kind=dp ), real_buff, bc_flag_viscous, delta,  Rq_x, Rq_z )
     call apply_bc( aimag( q(:, ii) ),        imag_buff, bc_flag_viscous, delta,  Iq_x, Iq_z )

     ! Store the result in the output vector in the right place.
     Aq(:, ii) = cmplx( real_buff, imag_buff, kind=dp )

  enddo

end subroutine apply_3D_viscous_fourier_single_velocity

subroutine solve_3D_viscous_single_velocity( u_int, u_b, bc_flag, iterations, residual )
! Solves the 3D viscous equations with a 2D SMPM discretization and a 1D Fourier
! discretization.

  use constants, only:         nsg, nsuby, nky, nu
  use options, only:           gmres_maxit_viscous, gmres_tol_viscous
  use woodbury_matrices, only: rpk
  use parallel_linear_algebra
  use precision, only:         dp
  use timestepping, only:      dt, g0

  implicit none

  ! Set the input/output variables.
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(inout) :: u_int
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(in)    :: u_b
  integer, dimension(1:4), intent(in)                     :: bc_flag
  integer, intent(out)                                    :: iterations
  real(kind=dp), intent(out)                              :: residual

  ! Arrays for storing Fourier transformed velocities.
  complex(kind=dp), dimension(1:rpk, 1:nky)               :: Fq, Fq_int

  ! Other variables.
  real(kind=dp)                                           :: gmres_tol
  integer                                                 :: gmres_viscous_iterations
  integer                                                 :: gmres_restart

  ! Declare the function used with GMRES to solve the capacitance problem.
  external                                                :: apply_3D_viscous_fourier_single_velocity

  ! Set some constants.
  gmres_tol                = gmres_tol_viscous
  gmres_viscous_iterations = min( gmres_maxit_viscous, nsg - 1 )
  gmres_restart            = gmres_viscous_iterations

  ! Set the viscous boundary conditions on the x and z boundaries.
  call apply_viscous_bc( u_int, u_b, bc_flag, nu * dt / g0 )

  ! Fourier transform the input velocities.
  call apply_fft( u_int, Fq_int )

  ! Set the initial guess to zero.
  Fq = 0.0_dp

  ! Solve the decoupled Helmholtz problems in Fourier space using GMRES.
  call compute_gmres_householder_complex( Fq, -Fq_int / g0, rpk * nky, &
                                          gmres_tol, gmres_viscous_iterations, gmres_restart, &
                                          apply_3D_viscous_fourier_single_velocity )

  ! Inverse Fourier transform the solution.
  call apply_ifft( Fq, u_int )

  ! Set the constants related to GMRES for output.
  iterations = gmres_viscous_iterations
  residual   = gmres_tol

end subroutine solve_3D_viscous_single_velocity
