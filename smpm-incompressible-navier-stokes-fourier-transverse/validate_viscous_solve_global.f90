subroutine validate_viscous_solve_global( output_prefix, tolerance, success_flag )
! Validates the viscous solver using a real-valued GMRES by Fourier transforming
! the RHS and storing the real and complex part separately in real datatype arrays.

  use constants, only:               nsuby, nky, nu
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
  real(kind=dp)                                  :: linf_ux_rhs, linf_uy_rhs, linf_uz_rhs, &
                                                    linf_ux_sol, linf_uy_sol, linf_uz_sol

  ! Variables for testing real/complex GMRES implementation.
  real(kind=dp), allocatable, dimension(:, :)    :: ux_ana, uy_ana, uz_ana
  real(kind=dp), allocatable, dimension(:, :)    :: ux_num, uy_num, uz_num

  real(kind=dp), allocatable, dimension(:, :)    :: rhs_ana_ux, rhs_ana_uy, rhs_ana_uz
  real(kind=dp), allocatable, dimension(:, :)    :: rhs_num_ux, rhs_num_uy, rhs_num_uz

  real(kind=dp), allocatable, dimension(:, :)    :: Fqrhs_real, Fqrhs_imag
  real(kind=dp), allocatable, dimension(:, :)    :: Fq_int_real, Fq_int_imag

  complex(kind=dp), allocatable, dimension(:, :) :: Fq_int, Fqrhs, Fux_rhs, Fuy_rhs, Fuz_rhs
  complex(kind=dp), allocatable, dimension(:, :) :: Fux, Fuy, Fuz

  character(len=256)                             :: caststr

  ! Other variables
  integer                                        :: ii
  external                                       :: apply_3D_viscous_fourier


  ! Allocate arrays
  allocate( ux_ana(1:rpk, 1:nsuby), uy_ana(1:rpk, 1:nsuby), uz_ana(1:rpk, 1:nsuby) )
  allocate( ux_num(1:rpk, 1:nsuby), uy_num(1:rpk, 1:nsuby), uz_num(1:rpk, 1:nsuby) )

  allocate( Fux_rhs(1:rpk, 1:nky), Fuy_rhs(1:rpk, 1:nky), Fuz_rhs(1:rpk, 1:nky) )
  allocate( Fux(1:rpk, 1:nky), Fuy(1:rpk, 1:nky), Fuz(1:rpk, 1:nky) )

  allocate( Fq_int(1:3*rpk, 1:nky), Fq_int_real(1:3*rpk, 1:nky), Fq_int_imag(1:3*rpk,1:nky ) )
  allocate( Fqrhs(1:3*rpk, 1:nky), Fqrhs_real(1:3*rpk,1:nky), Fqrhs_imag(1:3*rpk,1:nky) )
  allocate( rhs_ana_ux(1:rpk, 1:nsuby), rhs_ana_uy(1:rpk, 1:nsuby), rhs_ana_uz(1:rpk, 1:nsuby) )
  allocate( rhs_num_ux(1:rpk, 1:nsuby), rhs_num_uy(1:rpk, 1:nsuby), rhs_num_uz(1:rpk, 1:nsuby) )

  ! Get pi.
  pi = 2.0_dp * acos( 0.0_dp )

  ! Set the timestepping coefficient. 
  g0 = 1.0_dp

  ! Set wavenumber
  ! This is an arbitrary parameter. The higher the wavenumber the more resolution we must have
  kk = 2.0_dp * pi
  ll = 2.0_dp * pi
  mm = 2.0_dp * pi

  ! Set the velocity arrays.

  ! Check for transverse domain
  if (nsuby > 1) then

     ! Dirichlet boundary conditions.
     ux_ana = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     uy_ana = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     uz_ana = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )

     ! Neumann boundary conditions.
     !ux_ana = cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz )
     !uy_ana = cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz )
     !uz_ana = cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz )

  else

     ! Dirichlet boundary conditions.
     ux_ana = sin( kk * cx ) * sin( mm * cz )
     uy_ana = 0.0_dp
     uz_ana = sin( kk * cx ) * sin( mm * cz )

     ! Neumann boundary conditions.
     !ux_ana = cos( kk * cx ) * cos( mm * cz )
     !uy_ana = 0.0_dp
     !uz_ana = cos( kk * cx ) * cos( mm * cz )

  endif


  call notify( 'Solving via a global vector for all velocities with complex GMRES.' )

  ! Generate the right hand side.

  ! Fourier transform the input velocities.
  call apply_fft( ux_ana, Fux )
  call apply_fft( uy_ana, Fuy )
  call apply_fft( uz_ana, Fuz )

  ! Set up global array containing all three velocities.
  Fq_int(1:rpk, :)         = Fux
  Fq_int(rpk+1:2*rpk, :)   = Fuy
  Fq_int(2*rpk+1:3*rpk, :) = Fuz

  ! Store the real and complex part separate
  do ii = 1,nky
     Fq_int_real( :, ii ) = real(  Fq_int( :, ii), kind=dp )
     Fq_int_imag( :, ii ) = aimag( Fq_int( :, ii) )
  enddo

  ! Apply the operator to generate RHS for both the real and imaginary part
  call apply_3D_viscous_fourier( Fqrhs_real, Fq_int_real )
  call apply_3D_viscous_fourier( Fqrhs_imag, Fq_int_imag )

  ! Assemble RHS in Complex Buffer
  do ii = 1,nky
     Fqrhs( :, ii ) = cmplx( Fqrhs_real( :, ii ), Fqrhs_imag( :, ii ), kind=dp )
  enddo 

  ! Redistribute the velocity components from the RHS into separate arrays.
  Fux_rhs = Fqrhs(1:rpk, :)
  Fuy_rhs = Fqrhs(rpk+1:2*rpk, :)
  Fuz_rhs = Fqrhs(2*rpk+1:3*rpk, :)

  ! Inverse transform each RHS.
  call apply_ifft( Fux_rhs, rhs_num_ux )
  call apply_ifft( Fuy_rhs, rhs_num_uy )
  call apply_ifft( Fuz_rhs, rhs_num_uz )

  ! Check for transverse
  if (nsuby > 1) then
     ! Dirichlet boundary conditions.
     rhs_ana_ux = -(nu * dt / g0) * (kk**2.0_dp + ll**2.0_dp + mm**2.0_dp) * sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) &
                  -(sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) )

     rhs_ana_uy = -(nu * dt / g0) * (kk**2.0_dp + ll**2.0_dp + mm**2.0_dp) * sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) &
                  -(sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) )

     rhs_ana_uz = -(nu * dt / g0) * (kk**2.0_dp + ll**2.0_dp + mm**2.0_dp) * sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) &
                  -(sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) )

     ! Neumann boundary conditions.
     !   rhs_ana_ux = -(nu * dt / g0) * (kk**2.0_dp + ll**2.0_dp + mm**2.0_dp) * cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz ) &
     !                -(cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz ) )

     !   rhs_ana_uy = -(nu * dt / g0) * (kk**2.0_dp + ll**2.0_dp + mm**2.0_dp) * cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz ) &
     !                -(cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz ) )

     !   rhs_ana_uz = -(nu * dt / g0) * (kk**2.0_dp + ll**2.0_dp + mm**2.0_dp) * cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz ) &
     !                -(cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz ) )

  else
     ! Dirichlet boundary conditions.
     rhs_ana_ux = -(nu * dt / g0) * (kk**2.0_dp + mm**2.0_dp) * sin( kk * cx ) * sin( mm * cz ) - (sin( kk * cx ) * sin( mm * cz ) )

     rhs_ana_uy = 0.0_dp

     rhs_ana_uz = -(nu * dt / g0) * (kk**2.0_dp + mm**2.0_dp) * sin( kk * cx ) * sin( mm * cz ) - (sin( kk * cx ) * sin( mm * cz ) )

     ! Neumann boundary conditions.
     !   rhs_ana_ux = -(nu * dt / g0) * (kk**2.0_dp + mm**2.0_dp) * cos( kk * cx ) * cos( mm * cz ) -(cos( kk * cx ) * cos( mm * cz ) )

     !   rhs_ana_uy = 0.0_dp

     !   rhs_ana_uz = -(nu * dt / g0) * (kk**2.0_dp + mm**2.0_dp) * cos( kk * cx ) * cos( mm * cz ) -(cos( kk * cx ) * cos( mm * cz ) )


   endif

   ! Report errors between analytical RHS and numerical RHS.
   linf_ux_rhs = pmaxval( reshape( abs( rhs_ana_ux - rhs_num_ux ), (/rpk * nsuby/) ) )
   linf_uy_rhs = pmaxval( reshape( abs( rhs_ana_uy - rhs_num_uy ), (/rpk * nsuby/) ) )
   linf_uz_rhs = pmaxval( reshape( abs( rhs_ana_uz - rhs_num_uz ), (/rpk * nsuby/) ) )
   write( caststr, '(A,D17.10,D17.10,D17.10)' ) &
        'Error of operator for global vector linear operator (ux,uy,uz):', linf_ux_rhs, linf_uy_rhs, linf_uz_rhs
   call notify( caststr )

   ! Viscous solve.
   call solve_3D_viscous( rhs_num_ux, rhs_num_uy, rhs_num_uz )

   ! Correct the sign due to the time-stepping integration
   ux_num = -rhs_num_ux
   uy_num = -rhs_num_uy
   uz_num = -rhs_num_uz

   linf_ux_sol = pmaxval( reshape( abs( ux_ana - ux_num ), (/rpk * nsuby/) ) )
   linf_uy_sol = pmaxval( reshape( abs( uy_ana - uy_num ), (/rpk * nsuby/) ) )
   linf_uz_sol = pmaxval( reshape( abs( uz_ana - uz_num ), (/rpk * nsuby/) ) )

   write( caststr, '(A,D17.10,D17.10,D17.10)' ) &
        'Error in single vector (ux,uy,uz):', linf_ux_sol, linf_uy_sol, linf_uz_sol
   call notify( caststr )

   ! Write the solution to file.
   !open( 65, file=trim( output_prefix ) // '_viscous_solution_single_vector.txt' )
   !write( 65, * ) ux_num
   !write( 65, * ) uy_num
   !write( 65, * ) uz_num
   !close( 65 )

   success_flag = ((linf_ux_rhs        < tolerance) .and. &
                   (linf_uy_rhs        < tolerance) .and. &
                   (linf_uz_rhs        < tolerance) .and. &
                   (linf_ux_sol < tolerance) .and. &
                   (linf_uy_sol < tolerance) .and. &
                   (linf_uz_sol < tolerance))

   ! Reset the timestepping coefficient
   g0 = 11.0_dp/6.0_dp

end subroutine validate_viscous_solve_global
