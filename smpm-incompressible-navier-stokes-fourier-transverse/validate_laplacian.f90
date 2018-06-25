subroutine validate_laplacian( tolerance, success_flag )
! Validates the laplacian function in 2D and 3D using a sinusoidal 
! with a user-specified wavenumber. 

  use constants, only:               nsuby, rank, root
  use geom, only:                    cx, cz
  use mpi, only:                     MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_MAX
  use parallel_linear_algebra
  use precision, only:               dp
  use transverse, only:              cy, Ly
  use woodbury_matrices, only:       rpk

  implicit none

  real(kind=dp), intent(in)                   :: tolerance
  logical, intent(out)                        :: success_flag

  integer                                     :: ierr

  real(kind=dp)                               :: pi, kk, ll, mm
  real(kind=dp)                               :: linf_error, l2_error
  real(kind=dp)                               :: linf_error_max, l2_error_max

  real(kind=dp), allocatable, dimension(:, :) :: phi, lapphi_ana, lapphi_num

  character(len=256)                          :: caststr

  real(kind=dp)                               :: Lx, Lz

  allocate( phi(1:rpk, 1:nsuby), lapphi_ana(1:rpk, 1:nsuby), lapphi_num(1:rpk, 1:nsuby) )

  ! Get pi.
  pi = 2.0_dp * acos( 0.0_dp )

  ! Compute domain lengths.
  Lx = pmaxval( reshape( cx, (/rpk * nsuby/) ) ) - pminval( reshape( cx, (/rpk * nsuby/) ) )
  Lz = pmaxval( reshape( cz, (/rpk * nsuby/) ) ) - pminval( reshape( cz, (/rpk * nsuby/) ) )

  ! Set wavenumber
  kk = 2.0_dp * pi / Lx
  if (nsuby > 1) then
     ll = 2.0_dp * pi / Ly
  else
     ll = 0.0_dp
  end if
  mm = 2.0_dp * pi / Lz

  ! Set the function values and check for transverse direction
  if (nsuby > 1) then
     phi = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
  else
     phi = sin( kk * cx ) * sin( mm * cz )
  endif

  ! Check for transverse domain and define analytical laplacian
  if (nsuby > 1) then
     lapphi_ana = - (kk**2  + ll**2  + mm**2) * sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) 
  else
     lapphi_ana = - (kk**2.0 + mm**2.0) * sin( kk * cx ) * sin( mm * cz )
  endif

  ! Compute numerical laplacian
  call compute_3D_laplacian( phi, lapphi_num )

  ! Compute errors.
  linf_error = pmaxval( reshape( abs(lapphi_ana - lapphi_num) , (/rpk * nsuby/) )  )
  l2_error = pnorm2( reshape( abs( lapphi_ana - lapphi_num ) , (/rpk * nsuby/) )  ) / &
             pnorm2( reshape( abs( lapphi_ana ) , (/rpk * nsuby/) )  )

  ! Compute the maximum error over all ranks and distribute across all ranks
  call MPI_ALLREDUCE( linf_error, linf_error_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  ! Compute the maximum error over all ranks and distribute across all ranks
  call MPI_ALLREDUCE( l2_error, l2_error_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  write( caststr, '(A,D17.10)' ) 'Error in laplacian:', linf_error_max
  call notify( caststr )

  success_flag = ((linf_error_max < tolerance))

  if (rank .eq. root) then
     open( 65, file='laplacian_error.txt' )
     write(65, *) linf_error_max
     write(65, *) l2_error_max
     close(65)
  endif


end subroutine validate_laplacian
