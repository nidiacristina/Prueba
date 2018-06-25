subroutine validate_helmholtz_3d_scalar_neumann( tolerance, success_flag )
! Validates the scalar solver with neumann bcs using a real-valued GMRES 
! by Fourier transforming the RHS and storing the real and complex part separately 
! in real datatype arrays. The diffusive solver is based on the the Helmholtz Equation, 
! which is given as,
!
!        d * Lap( phi ) - phi = f 
! 
!  where f is the RHS containing the density field, phi is the 
!  unknown density at rho^{n+1} and d is the coefficient contianing the timestep
!  and diffusion information. 
!
!  April 2017
!  Greg Thomsen & Gustavo Rivera

  use constants, only:               nsuby, nky, nu_d, rank, root
  use geom, only:                    cx, cz
  use mpi, only:                     MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_MAX
  use parallel_linear_algebra, only: pmaxval, pminval, pnorm2
  use precision, only:               dp
  use timestepping, only:            g0, dt
  use transverse, only:              cy, Ly
  use woodbury_matrices, only:       rpk

  implicit none

  real(kind=dp), intent(in)                     :: tolerance
  logical, intent(out)                          :: success_flag

  ! General use variables.
  real(kind=dp)                                 :: pi, mm, ll, kk
  real(kind=dp)                                 :: linf_rhs, linf_sol
  real(kind=dp)                                 :: linf_rhs_max, linf_sol_max
  real(kind=dp)                                 :: l2_rhs, l2_sol
  real(kind=dp)                                 :: l2_rhs_max, l2_sol_max

  ! Variables for testing real/complex GMRES implementation.
  real(kind=dp), allocatable, dimension(:, :)   :: rho_ana
  real(kind=dp), allocatable, dimension(:, :)   :: rho_num

  complex(kind=dp), allocatable, dimension(:,:) :: Frho_ana
  real(kind=dp), allocatable, dimension(:,:)    :: Frho_ana_real
  real(kind=dp), allocatable, dimension(:,:)    :: Frho_ana_imag

  complex(kind=dp), allocatable, dimension(:,:) :: Frho_num
  real(kind=dp), allocatable, dimension(:,:)    :: Frho_num_real
  real(kind=dp), allocatable, dimension(:,:)    :: Frho_num_imag

  complex(kind=dp), allocatable, dimension(:,:) :: Frhs_num
  real(kind=dp), allocatable, dimension(:, :)   :: Frhs_num_real
  real(kind=dp), allocatable, dimension(:, :)   :: Frhs_num_imag

  real(kind=dp), allocatable, dimension(:, :)   :: rhs_ana
  real(kind=dp), allocatable, dimension(:, :)   :: rhs_num

  character(len=256)                            :: caststr
  integer                                       :: ierr

  real(kind=dp)                                 :: Lx, Lz

  ! Other variables
  integer                                       :: ii


  ! Allocate arrays
  allocate( rho_ana(1:rpk, 1:nsuby) )
  allocate( rho_num(1:rpk, 1:nsuby) )
 
  allocate( Frho_ana(1:rpk,1:nky ) )
  allocate( Frho_ana_real(1:rpk,1:nky ) )
  allocate( Frho_ana_imag(1:rpk,1:nky ) )

  allocate( Frho_num(1:rpk,1:nky ) )
  allocate( Frho_num_real(1:rpk,1:nky ) )
  allocate( Frho_num_imag(1:rpk,1:nky ) )

  allocate( rhs_ana(1:rpk, 1:nsuby) )
  allocate( rhs_num(1:rpk, 1:nsuby) )

  allocate( Frhs_num(1:rpk,1:nky ) )
  allocate( Frhs_num_real(1:rpk,1:nky ) )
  allocate( Frhs_num_imag(1:rpk,1:nky ) )

  ! Get pi.
  pi = 2.0_dp * acos( 0.0_dp )

  ! Compute domain lengths.
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
  linf_rhs = 0.0_dp
  linf_sol = 0.0_dp

  ! Set the velocity arrays.

  ! Check for transverse domain
  if (nsuby > 1) then

     ! Neumann boundary conditions.
     rho_ana = cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz )

  else

     ! Neumann boundary conditions.
     rho_ana = cos( kk * cx ) * cos( mm * cz )

  endif

  ! Fourier transform the input velocities.
  call apply_fft( rho_ana, Frho_ana )

  ! Store the real and complex part separate
  do ii = 1,nky
     Frho_ana_real( :, ii ) = real(  Frho_ana( :, ii) , kind=dp )
     Frho_ana_imag( :, ii ) = aimag( Frho_ana( :, ii) )
  enddo

  ! Apply the operator to generate RHS for both the real and imaginary part
  call apply_3D_diffusive_fourier( Frhs_num_real, Frho_ana_real )
  call apply_3D_diffusive_fourier( Frhs_num_imag, Frho_ana_imag )

  ! Assemble RHS in Complex Buffer
  do ii = 1,nky
     Frhs_num( :, ii ) = cmplx( Frhs_num_real( :, ii ), Frhs_num_imag( :, ii ) , kind=dp )
  enddo 

  ! Inverse transform each RHS.
  call apply_ifft( Frhs_num, rhs_num )

  ! Check for transverse
  if (nsuby > 1) then

     ! Neumann boundary conditions.
     rhs_ana = -(nu_d * dt / g0) * (kk**2.0 + ll**2.0 + mm**2.0) * cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz ) &
               -(cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz ) )

  else

     ! Neumann boundary conditions.
     rhs_ana = -(nu_d * dt / g0) * (kk**2.0 + mm**2.0) * cos( kk * cx ) * cos( mm * cz ) -(cos( kk * cx ) * cos( mm * cz ) )

   endif

   ! Report errors between analytical RHS and numerical RHS.
   linf_rhs = pmaxval( reshape( abs( rhs_ana - rhs_num ), (/rpk * nsuby/) ) )
   l2_rhs   = pnorm2( reshape( abs( rhs_ana - rhs_num ), (/rpk * nsuby/) ) ) /&
              pnorm2( reshape( rhs_ana, (/rpk * nsuby/) ) )

  ! Compute the maximum error over all ranks and distribute across all ranks 
  call MPI_ALLREDUCE( linf_rhs, linf_rhs_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( l2_rhs, l2_rhs_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

   write( caststr, '(A,D17.10,D17.10,D17.10)' ) &
        'Error of operator in diffusive solve:', linf_rhs_max
   call notify( caststr )

   ! Diffusive solve.
   call solve_3D_diffusion( rhs_num )

   ! Correct the sign due to the time-stepping integration
   rho_num = -rhs_num

   linf_sol = pmaxval( reshape( abs( rho_ana - rho_num ), (/rpk * nsuby/) ) )
   l2_sol   = pnorm2( reshape( abs( rho_ana - rho_num ), (/rpk * nsuby/) ) ) /&
              pnorm2( reshape( rho_ana , (/rpk * nsuby/) ) )

   ! Reset the timestepping coefficient
   g0 = 11.0_dp/6.0_dp

  ! Compute the maximum error over all ranks and distribute across all ranks
  call MPI_ALLREDUCE( linf_sol, linf_sol_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( l2_sol, l2_sol_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

   success_flag = (linf_sol_max        < tolerance)

   write( caststr, '(A,D17.10,D17.10,D17.10)' ) &
        'Error in diffusive solve:', linf_sol_max
   call notify( caststr )



   if (rank .eq. root) then
      open(unit=65,file='HSN_error.txt')
      write(65,*) linf_sol_max
      write(65,*) linf_rhs_max
      write(65,*) l2_sol_max
      write(65,*) l2_rhs_max
      close(65)
   endif

end subroutine validate_helmholtz_3d_scalar_neumann
