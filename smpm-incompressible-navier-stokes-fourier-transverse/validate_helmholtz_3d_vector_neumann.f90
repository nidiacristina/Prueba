subroutine validate_helmholtz_3d_vector_neumann( tolerance, success_flag )
! Validates the viscous solver using a real-valued GMRES by Fourier transforming
! the RHS and storing the real and complex part separately in real datatype arrays.
! The viscous solver is based on the the Helmholtz Equation, which is given as,
!
!        d * Lap( phi ) - phi = f 
! 
!  where f is the RHS containing the concatenated velocity field, phi is the 
!  unknown velocity at u^{n+1} and d is the coefficient contianing the timestep
!  and viscous information. 
!
!  April 2017
!  Greg Thomsen & Gustavo Rivera

  use constants, only:               nsuby, nky, nu, rank, root
  use geom, only:                    cx, cz
  use mpi, only:                     MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_MAX
  use parallel_linear_algebra, only: pdot_product, pmaxval, pminval, pnorm2, znorm2
  use precision, only:               dp
  use timestepping, only:            g0, dt
  use transverse, only:              cy, Ly
  use woodbury_matrices, only:       rpk

  implicit none

  real(kind=dp), intent(in)                      :: tolerance
  logical, intent(out)                           :: success_flag

  ! General use variables.
  real(kind=dp)                                  :: pi, mm, ll, kk
  real(kind=dp)                                  :: linf_ux_rhs, linf_uy_rhs, linf_uz_rhs, &
                                                    linf_ux_sol, linf_uy_sol, linf_uz_sol
  real(kind=dp)                                  :: linf_ux_rhs_max, linf_uy_rhs_max, linf_uz_rhs_max, &
                                                    linf_ux_sol_max, linf_uy_sol_max, linf_uz_sol_max
  real(kind=dp)                                  :: l2_ux_rhs, l2_uy_rhs, l2_uz_rhs, &
                                                    l2_ux_sol, l2_uy_sol, l2_uz_sol
  real(kind=dp)                                  :: l2_ux_rhs_max, l2_uy_rhs_max, l2_uz_rhs_max, &
                                                    l2_ux_sol_max, l2_uy_sol_max, l2_uz_sol_max


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
  integer                                        :: ierr

  real(kind=dp)                                  :: Lx, Lz

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

  ! Compute domain lengths
  Lx = pmaxval( reshape( cx, (/rpk * nsuby/) ) ) - pminval( reshape( cx, (/rpk * nsuby/) ) )
  Lz = pmaxval( reshape( cz, (/rpk * nsuby/) ) ) - pminval( reshape( cz, (/rpk * nsuby/) ) )

  ! Set the timestepping coefficient. 
  g0 = 1.0_dp

  ! Set wavenumber
  ! This is an arbitrary parameter. The higher the wavenumber the more resolution we must have
  kk = 2.0_dp * pi / Lx
  if (nsuby > 1) then
     ll = 2.0_dp * pi / Ly
  else
     ll = 0.0_dp
  end if
  mm = 2.0_dp * pi / Lz

  ! Initialize error variables
  linf_ux_rhs = 0.0_dp
  linf_uy_rhs = 0.0_dp
  linf_uz_rhs = 0.0_dp
  linf_ux_sol = 0.0_dp
  linf_uy_sol = 0.0_dp
  linf_uz_sol = 0.0_dp

  ! Set the velocity arrays.

  ! Check for transverse domain
  if (nsuby > 1) then

     ! Neumann boundary conditions.
     ux_ana = sin( kk * cx ) * sin( ll * cy ) * cos( mm * cz )
     uy_ana = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     uz_ana = cos( kk * cx ) * sin( ll * cy ) * sin( mm * cz )

  else

     ! Neumann boundary conditions.
     ux_ana = sin( kk * cx ) * cos( mm * cz )
     uy_ana = 0.0_dp
     uz_ana = cos( kk * cx ) * sin( mm * cz )

  endif


  call notify( 'Solving via a global vector for all velocities with complex GMRES.' )

  ! Generate the right hand side.

  ! Fourier transform the input velocities.
  call apply_fft( ux_ana, Fux )
  call apply_fft( uy_ana, Fuy )
  call apply_fft( uz_ana, Fuz )

  ! Set up global array containing all three velocities.
  Fq_int(1:rpk, 1:nky)         = Fux(1:rpk,1:nky)
  Fq_int(rpk+1:2*rpk, 1:nky)   = Fuy(1:rpk,1:nky)
  Fq_int(2*rpk+1:3*rpk, 1:nky) = Fuz(1:rpk,1:nky)

  ! Store the real and complex part separate
  do ii = 1,nky
     Fq_int_real( :, ii ) = real(  Fq_int( :, ii) , kind=dp )
     Fq_int_imag( :, ii ) = aimag( Fq_int( :, ii) )
  enddo

  ! Apply the operator to generate RHS for both the real and imaginary part
  call apply_3D_viscous_fourier( Fqrhs_real, Fq_int_real )
  call apply_3D_viscous_fourier( Fqrhs_imag, Fq_int_imag )

  ! Assemble RHS in Complex Buffer
  do ii = 1,nky
     Fqrhs( :, ii ) = cmplx( Fqrhs_real( :, ii ), Fqrhs_imag( :, ii ) , kind=dp )
  enddo 

  ! Redistribute the velocity components from the RHS into separate arrays.
  Fux_rhs(1:rpk,1:nky) = Fqrhs(1:rpk,1:nky)
  Fuy_rhs(1:rpk,1:nky) = Fqrhs(rpk+1:2*rpk,1:nky)
  Fuz_rhs(1:rpk,1:nky) = Fqrhs(2*rpk+1:3*rpk,1:nky)

  ! Inverse transform each RHS.
  call apply_ifft( Fux_rhs, rhs_num_ux )
  call apply_ifft( Fuy_rhs, rhs_num_uy )
  call apply_ifft( Fuz_rhs, rhs_num_uz )

  ! Check for transverse
  if (nsuby > 1) then

     ! Neumann boundary conditions.
     rhs_ana_ux = -(nu * dt / g0) * (kk**2.0 + ll**2.0 + mm**2.0) * sin( kk * cx ) * sin( ll * cy ) * cos( mm * cz ) &
                  -(sin( kk * cx ) * sin( ll * cy ) * cos( mm * cz ) )

     rhs_ana_uy = -(nu * dt / g0) * (kk**2.0 + ll**2.0 + mm**2.0) * sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) &
                  -(sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) )

     rhs_ana_uz = -(nu * dt / g0) * (kk**2.0 + ll**2.0 + mm**2.0) * cos( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) &
                  -(cos( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) )

  else

     ! Neumann boundary conditions.
     rhs_ana_ux = -(nu * dt / g0) * (kk**2.0 + mm**2.0) * sin( kk * cx ) * cos( mm * cz ) -(sin( kk * cx ) * cos( mm * cz ) )

     rhs_ana_uy = 0.0_dp

     rhs_ana_uz = -(nu * dt / g0) * (kk**2.0 + mm**2.0) * cos( kk * cx ) * sin( mm * cz ) -(cos( kk * cx ) * sin( mm * cz ) )

   endif

   ! Report errors between analytical RHS and numerical RHS.
   linf_ux_rhs = pmaxval( reshape( abs( rhs_ana_ux - rhs_num_ux ), (/rpk * nsuby/) ) )
   linf_uy_rhs = pmaxval( reshape( abs( rhs_ana_uy - rhs_num_uy ), (/rpk * nsuby/) ) )
   linf_uz_rhs = pmaxval( reshape( abs( rhs_ana_uz - rhs_num_uz ), (/rpk * nsuby/) ) )

   l2_ux_rhs = pnorm2( reshape( abs( rhs_ana_ux - rhs_num_ux ), (/rpk * nsuby/) ) )
   if (nsuby .gt. 1) then
      l2_uy_rhs = pnorm2( reshape( abs( rhs_ana_uy - rhs_num_uy ), (/rpk * nsuby/) ) )/ &
                  pnorm2( reshape( abs( rhs_ana_uy ), (/rpk * nsuby/) ) )
   else 
      l2_uy_rhs = 0.0_dp
   endif
   l2_uz_rhs = pnorm2( reshape( abs( rhs_ana_uz - rhs_num_uz ), (/rpk * nsuby/) ) )


   ! Compute the maximum error over all ranks and distribute across all ranks 
   call MPI_ALLREDUCE( linf_ux_rhs, linf_ux_rhs_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
   call MPI_ALLREDUCE( linf_uy_rhs, linf_uy_rhs_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
   call MPI_ALLREDUCE( linf_uz_rhs, linf_uz_rhs_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
   call MPI_ALLREDUCE( l2_ux_rhs, l2_ux_rhs_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
   call MPI_ALLREDUCE( l2_uy_rhs, l2_uy_rhs_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
   call MPI_ALLREDUCE( l2_uz_rhs, l2_uz_rhs_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

   write( caststr, '(A,D17.10,D17.10,D17.10)' ) &
        'Error of operator for global vector linear operator (ux,uy,uz):', linf_ux_rhs_max, linf_uy_rhs_max, linf_uz_rhs_max
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

   l2_ux_sol = pnorm2( reshape( abs( ux_ana - ux_num ), (/rpk * nsuby/) ) ) /&
               pnorm2( reshape( abs( ux_ana ), (/rpk * nsuby/) ) ) 

   if (nsuby .gt. 1 ) then
      l2_uy_sol = pnorm2( reshape( abs( uy_ana - uy_num ), (/rpk * nsuby/) ) ) /&
                  pnorm2( reshape( abs( uy_ana ), (/rpk * nsuby/) ) ) 
   else 
      l2_uy_sol = 0.0_dp
   endif

   l2_uz_sol = pnorm2( reshape( abs( uz_ana - uz_num ), (/rpk * nsuby/) ) ) /&
               pnorm2( reshape( abs( uz_ana ), (/rpk * nsuby/) ) ) 

   ! Compute the maximum error over all ranks and distribute across all ranks 
   call MPI_ALLREDUCE( linf_ux_sol, linf_ux_sol_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
   call MPI_ALLREDUCE( linf_uy_sol, linf_uy_sol_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
   call MPI_ALLREDUCE( linf_uz_sol, linf_uz_sol_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

   call MPI_ALLREDUCE( l2_ux_sol, l2_ux_sol_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
   call MPI_ALLREDUCE( l2_uy_sol, l2_uy_sol_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
   call MPI_ALLREDUCE( l2_uz_sol, l2_uz_sol_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

   write( caststr, '(A,D17.10,D17.10,D17.10)' ) &
        'Error in single vector (ux,uy,uz):', linf_ux_sol_max, linf_uy_sol_max, linf_uz_sol_max
   call notify( caststr )

   success_flag = ((linf_ux_sol_max < tolerance) .and. &
                   (linf_uy_sol_max < tolerance) .and. &
                   (linf_uz_sol_max < tolerance))

   ! Reset the timestepping coefficient
   g0 = 11.0_dp/6.0_dp

   if (rank .eq. root ) then
      open(unit=65,file='HVN_error.txt')
      write(65,*) linf_ux_sol_max
      write(65,*) linf_uy_sol_max
      write(65,*) linf_uz_sol_max
      write(65,*) linf_ux_rhs_max
      write(65,*) linf_uy_rhs_max
      write(65,*) linf_uz_rhs_max

      write(65,*) l2_ux_sol_max
      write(65,*) l2_uy_sol_max
      write(65,*) l2_uz_sol_max
      write(65,*) l2_ux_rhs_max
      write(65,*) l2_uy_rhs_max
      write(65,*) l2_uz_rhs_max
      close(65)
   endif

end subroutine validate_helmholtz_3d_vector_neumann
