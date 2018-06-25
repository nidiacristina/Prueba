subroutine validate_divergence( tolerance, success_flag )
! Validates the divergence function in 2D and 3D using a sinusoidal 
! with a user-specified wavenumber. 

  use constants, only:               nsuby, rank, root
  use geom, only:                    cx, cz
  use mpi, only:                     MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_MAX
  use field_variables, only:         ux, uy, uz
  use parallel_linear_algebra
  use precision, only:               dp
  use transverse, only:              cy, Ly
  use woodbury_matrices, only:       rpk

  implicit none

  real(kind=dp), intent(in)                   :: tolerance
  logical, intent(out)                        :: success_flag

  real(kind=dp)                               :: pi, kk, ll, mm
  real(kind=dp)                               :: linf_error, linf_error_max
  real(kind=dp)                               :: l2_error, l2_error_max

  integer                                     :: ierr

  real(kind=dp), allocatable, dimension(:, :) :: divu_ana, divu_num
  real(kind=dp), allocatable, dimension(:, :) :: dudx_ana, dvdy_ana, dwdz_ana

  character(len=256)                          :: caststr

  real(kind=dp)                               :: Lx, Lz

  allocate( divu_ana(1:rpk, 1:nsuby), divu_num(1:rpk,1:nsuby) )
  allocate( dudx_ana(1:rpk, 1:nsuby), dvdy_ana(1:rpk,1:nsuby), dwdz_ana(1:rpk,1:nsuby) )

  ! Get pi.
  pi = 2.0_dp * acos( 0.0_dp )

  ! Get the domain lengths
  Lx = pmaxval( reshape( cx , (/rpk * nsuby/) ) ) - pminval( reshape( cx, (/rpk * nsuby/) ) )
  Lz = pmaxval( reshape( cz , (/rpk * nsuby/) ) ) - pminval( reshape( cz, (/rpk * nsuby/) ) )

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
     ux = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     uy = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     uz = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
  else
     ux = sin( kk * cx ) * sin( mm * cz )
     uy = sin( kk * cx ) * sin( mm * cz )
     uz = sin( kk * cx ) * sin( mm * cz )
  endif

  ! Compute analytical gradients
  if (nsuby > 1) then
     dudx_ana = kk * cos( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     dvdy_ana = ll * sin( kk * cx ) * cos( ll * cy ) * sin( mm * cz )
     dwdz_ana = mm * sin( kk * cx ) * sin( ll * cy ) * cos( mm * cz )
  else
     dudx_ana = kk * cos( kk * cx ) * sin( mm * cz )
     dvdy_ana = 0.0_dp
     dwdz_ana = mm * sin( kk * cx ) * cos( mm * cz )
  endif

  ! Compute the analytical divergence
  divu_ana = dudx_ana + dvdy_ana + dwdz_ana

  ! Compute the numerical divergence
  call compute_3D_divergence( ux, uy, uz, divu_num )

  ! Compute errors.
  linf_error = pmaxval( reshape( abs( divu_ana - divu_num ), (/rpk * nsuby/) )  )
  l2_error = pnorm2( reshape( abs( divu_ana - divu_num ), (/rpk * nsuby/) )  ) /&
             pnorm2( reshape( divu_ana , (/rpk * nsuby/) )  )


  ! Compute the maximum error over all ranks and distribute across all ranks
  call MPI_ALLREDUCE( linf_error, linf_error_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( l2_error, l2_error_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )


  write( caststr, '(A,D17.10)' ) 'Error in divergence:', linf_error_max
  call notify( caststr )

  success_flag = (linf_error_max < tolerance)

  if (rank .eq. root) then
     open(unit=65,file='divergence_error.txt')
     write(65,*) linf_error_max
     write(65,*) l2_error_max
     close(65)
  endif

end subroutine validate_divergence
