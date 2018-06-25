subroutine validate_poisson_3d_iteratively( tolerance, success_flag )
! Validates the poisson equation in 3d with deflation using a real-valued GMRES 
! by Fourier transforming the RHS and storing the real and complex part separately 
! in real datatype arrays. The poisson equation is given as, 
!
!        Lap( phi ) = f 
! 
!  where f is the RHS (divergence of the velocity) and phi os the unknown (the pressure) 
!
!  April 2017
!  Greg Thomsen & Gustavo Rivera


  use constants, only:         nsg, nsuby, nky, rank, root
  use geom, only:              cx, cz
  use mpi, only:               MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_MAX
  use transverse, only:        cy, Ly
  use parallel_linear_algebra
  use precision, only:         dp
  use woodbury_matrices, only: rpk


  implicit none

  real(kind=dp), intent(in)                                  :: tolerance
  logical, intent(out)                                       :: success_flag

  ! Arrays for store the analytical and numerical phi
  real(kind=dp), allocatable, dimension(:,:)                 :: phi_ana
  real(kind=dp), allocatable, dimension(:,:)                 :: phi_num

  ! Arrays for store the analytical and numerical rhs
  real(kind=dp), allocatable, dimension(:,:)                 :: rhs_ana
  real(kind=dp), allocatable, dimension(:,:)                 :: rhs_num

  ! Arrays to store the real and imag right hand side
  real(kind=dp), allocatable, dimension(:,:)                 :: Frhs_num_real
  real(kind=dp), allocatable, dimension(:,:)                 :: Frhs_num_imag

  ! Array to store complex phi
  complex(kind=dp), allocatable, dimension(:,:)              :: Fphi_ana
  complex(kind=dp), allocatable, dimension(:,:)              :: Fphi_num

  ! Array to store complex RHS
  complex(kind=dp), allocatable, dimension(:,:)              :: Frhs_num

  ! Arraty to store the real and imag part of phu
  real(kind=dp), allocatable, dimension(:,:)                 :: Fphi_ana_real
  real(kind=dp), allocatable, dimension(:,:)                 :: Fphi_ana_imag

  ! Buffers with whole domain
  real(kind=dp), allocatable, dimension(:,:)                 :: phi_ana_full, phi_num_full
  real(kind=dp), allocatable, dimension(:,:)                 :: rhs_ana_full, rhs_num_full


  ! The wave numbers for the analytical solution
  real(kind=dp)                                              :: kk, mm, ll

  ! Constant pi
  real(kind=dp)                                              :: pi

  ! Error Variable
  real(kind=dp)                                              :: linf_rhs, l2_rhs
  real(kind=dp)                                              :: linf_rhs_max, l2_rhs_max
  real(kind=dp)                                              :: linf_solution, l2_solution
  real(kind=dp)                                              :: linf_solution_max, l2_solution_max

  ! Looping Variable
  integer                                                    :: ii

  integer                                                    :: ierr
  character(len=64)                                          :: caststr

  ! Domain Lengths
  real(kind=dp)                                              :: Lx, Lz

  ! GMRES Stats
  integer                                                    :: gmres_poisson_iterations_real
  integer                                                    :: gmres_poisson_iterations_imag
  real(kind=dp)                                              :: l2_poisson_schur_error_real
  real(kind=dp)                                              :: l2_poisson_schur_error_imag

  ! Allocate Arrays
  allocate( phi_ana( 1:rpk, 1:nsuby ) )
  allocate( phi_num( 1:rpk, 1:nsuby ) )

  allocate( rhs_ana( 1:rpk, 1:nsuby ) )
  allocate( rhs_num( 1:rpk, 1:nsuby ) )

  allocate( Fphi_ana(1:rpk, 1:nky ) )
  allocate( Fphi_ana_real(1:rpk,1:nky ) )
  allocate( Fphi_ana_imag(1:rpk,1:nky ) )

  allocate( Frhs_num(1:rpk, 1:nky ) )
  allocate( Frhs_num_real( 1:rpk, 1:nky ) )
  allocate( Frhs_num_imag( 1:rpk, 1:nky ) )

  allocate( Fphi_num( 1:rpk,1:nky ) )

  allocate( rhs_ana_full(1:nsg,1:nsuby))
  allocate( rhs_num_full(1:nsg,1:nsuby))

  allocate( phi_ana_full(1:nsg,1:nsuby))
  allocate( phi_num_full(1:nsg,1:nsuby))


  ! Compute the domain lengths
  Lx = pmaxval( reshape( cx, (/rpk * nsuby/) ) ) - pminval( reshape( cx, (/rpk * nsuby/) ) )
  Lz = pmaxval( reshape( cz, (/rpk * nsuby/) ) ) - pminval( reshape( cz, (/rpk * nsuby/) ) )

  ! Get pi.
  pi = 2.0_dp * acos( 0.0_dp )

  ! Set the wavenumber. This in an arbitrary parameter, the higher the wavenumber the more resolution
  ! we need
  kk = 2.0_dp * pi / Lx ! x
  if (nsuby > 1) then
     ll = 2.0_dp * pi / Ly
  else
     ll = 0.0_dp
  end if
  mm = 2.0_dp * pi / Lz ! z

  ! Initialize error variables
  linf_rhs = 0.0_dp
  linf_solution = 0.0_dp

  ! Generate phi
  if (nsuby > 1) then
     phi_ana  = cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz )
  else
     phi_ana  = cos( kk * cx ) * cos( mm * cz )
  endif

  ! Transform into complex space
  call apply_fft( phi_ana, Fphi_ana)
  
  ! Separate into real and imaginary components
  do ii = 1, nky
     Fphi_ana_real( :, ii ) = real(  Fphi_ana( :, ii ), kind=dp )
     Fphi_ana_imag( :, ii ) = aimag( Fphi_ana( :, ii ) )
  enddo

  ! Compute the numerical RHS
  call apply_3D_poisson(Frhs_num_real, Fphi_ana_real)
  call apply_3D_poisson(Frhs_num_imag, Fphi_ana_imag)

  ! Compute analytical RHS
  if (nsuby > 1) then
     rhs_ana = -( kk**2 + ll**2 + mm**2 ) * cos( kk * cx ) * cos( ll * cy ) * cos( mm * cz ) 
  else
     rhs_ana = -( kk**2 + mm**2 ) * cos( kk * cx ) * cos( mm * cz ) 
  endif

  ! Assemble the rhs in complex space for error check
  do ii = 1,nky
     Frhs_num( :, ii) = cmplx( Frhs_num_real( :, ii ), Frhs_num_imag( :, ii), kind=dp )
  enddo

  ! Inverse Fourier Transform the rhs
  call apply_ifft(Frhs_num, rhs_num)


  ! Display operator error.
  linf_rhs = pmaxval( reshape( abs( rhs_ana - rhs_num ), (/rpk * nsuby/) ) )
  l2_rhs   = pnorm2( reshape( abs( rhs_ana - rhs_num ), (/rpk * nsuby/) ) )/ pnorm2( reshape( rhs_ana, (/rpk * nsuby/) ) )


  ! Compute the maximum error over all ranks and distribute across all ranks
  call MPI_ALLREDUCE( linf_rhs, linf_rhs_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  ! Compute the maximum error over all ranks and distribute across all ranks
  call MPI_ALLREDUCE( l2_rhs, l2_rhs_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  write( caststr, '(A,D17.10)' ) 'L2 error in Laplacian Operator 3D:', l2_rhs_max
  call notify( caststr )

  write( caststr, '(A,D17.10)' ) 'Linf error in Laplacian Operator 3D:', linf_rhs_max
  call notify( caststr )

   ! Solve Poisson equation for the real coefficients
  call solve_3D_poisson_iteratively( Frhs_num_real, gmres_poisson_iterations_real, l2_poisson_schur_error_real, 'R' )

  ! Solve Poisson equation for the imag coefficients
  call solve_3D_poisson_iteratively( Frhs_num_imag, gmres_poisson_iterations_imag, l2_poisson_schur_error_imag, 'I' )

  ! Assemble the solution, phi, in complex space
  do ii = 1,nky
     Fphi_num( :, ii) = cmplx( Frhs_num_real( :, ii ), Frhs_num_imag( :, ii), kind=dp )
  enddo

  ! Inverse Fourier Transform the Solution
  call apply_ifft(Fphi_num, phi_num)

  ! Compute Error
  l2_solution = pnorm2( reshape( abs(phi_num - phi_ana), (/nsuby * rpk/) ) ) / pnorm2( reshape(  phi_ana , (/nsuby * rpk/) ) )
  linf_solution = pmaxval( reshape( abs(phi_num - phi_ana), (/nsuby * rpk/) ) )

  ! Compute the maximum error over all ranks and distribute across all ranks
  call MPI_ALLREDUCE( linf_solution, linf_solution_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  ! Compute the maximum error over all ranks and distribute across all ranks
  call MPI_ALLREDUCE( l2_solution, l2_solution_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  ! Notify Error
  write( caststr, '(A,D17.10)' ) 'Linf error in solution:', linf_solution_max
  call notify( caststr )

  write( caststr, '(A,D17.10)' ) 'L2 error in solution:', l2_solution_max
  call notify( caststr )

  write( caststr, '(A,D17.10)' ) 'L2 error in real solve:', l2_poisson_schur_error_real
  call notify( caststr )

  write( caststr, '(A,D17.10)' ) 'L2 error in imag solve:', l2_poisson_schur_error_imag
  call notify( caststr )

  write( caststr, '(A,I5)' ) 'Iterations in real solve:', gmres_poisson_iterations_real
  call notify( caststr )

  write( caststr, '(A,I5)' ) 'Iterations in imag solve:', gmres_poisson_iterations_imag
  call notify( caststr )

  success_flag = (l2_solution_max  < tolerance)
              
  ! Write to File
  if (rank .eq. root) then
     open( 65, file='pressure_error.txt' )
     write(65, *) l2_solution_max
     write(65, *) linf_solution_max
     write(65, *) l2_rhs_max
     write(65, *) linf_rhs_max
     write(65, *) gmres_poisson_iterations_real
     write(65, *) gmres_poisson_iterations_imag
     write(65, *) l2_poisson_schur_error_real
     write(65, *) l2_poisson_schur_error_imag
     close(65)
  endif

  ! Write Output to file
  call gather_3D_array( phi_num, phi_num_full)
  call gather_3D_array( phi_ana, phi_ana_full)
  call gather_3D_array( rhs_num, rhs_num_full)
  call gather_3D_array( rhs_ana, rhs_ana_full)

  if (rank .eq. root) then
     open(65, file='poisson_solution.txt')
     write(65,*) phi_num_full
     write(65,*) phi_ana_full
     write(65,*) rhs_num_full
     write(65,*) rhs_ana_full
   close(65)
  endif


end subroutine validate_poisson_3d_iteratively

