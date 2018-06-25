subroutine validate_transport_term( tolerance, success_flag )
! Validates the advective terms in X, Y, and Z.

  use constants, only:               nsuby, rank, root
  use field_variables, only:         Nrho0, rho, ux, uy, uz
  use geom, only:                    cx, cz
  use mpi, only:                     MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_MAX
  use parallel_linear_algebra
  use precision, only:               dp
  use transverse, only:              Ly, cy
  use woodbury_matrices, only:       rpk

  implicit none

  real(kind=dp), intent(in)                   :: tolerance
  logical, intent(out)                        :: success_flag

  real(kind=dp)                               :: pi, kk, ll, mm
  real(kind=dp)                               :: linf
  real(kind=dp)                               :: linf_max
  real(kind=dp)                               :: l2
  real(kind=dp)                               :: l2_max

  real(kind=dp), allocatable, dimension(:, :) :: analytic
  real(kind=dp), allocatable, dimension(:, :) :: drhodx, drhody, drhodz

  real(kind=dp)                               :: Lx, Lz

  character(len=256)                          :: caststr
  integer                                     :: ierr

  allocate( drhodx(1:rpk, 1:nsuby), drhody(1:rpk, 1:nsuby), drhodz(1:rpk, 1:nsuby) )
  allocate( analytic(1:rpk, 1:nsuby) )

  ! Get pi.
  pi = 2.0_dp * acos( 0.0_dp )

  ! Get the domain lengths
  Lx = pmaxval( reshape( cx, (/rpk * nsuby/) ) ) - pminval( reshape( cx, (/rpk * nsuby/) ) )
  Lz = pmaxval( reshape( cz, (/rpk * nsuby/) ) ) - pminval( reshape( cz, (/rpk * nsuby/) ) )

  ! Set coefficients.
  kk = 2.0_dp * pi / Lx
  if (nsuby > 1) then
     ll = 2.0_dp * pi / Ly
  else
     ll = 0.0_dp
  end if
  mm = 2.0_dp * pi / Lz

  ! Set the velocity values and check for transverse direction
  if (nsuby > 1) then
     ux = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     uy = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     uz = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
  else
     ux = sin( kk * cx ) * sin( mm * cz )
     uy = 0.0_dp
     uz = sin( kk * cx ) * sin( mm * cz )
  endif

  ! Set the density values and check for transverse direction
  if (nsuby > 1) then
     rho = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
  else
     rho = sin( kk * cx ) * sin( mm * cz )
  endif

  ! Compute the transport term
  call apply_smpm_transport( Nrho0, rho, ux, uy, uz )

  ! Compute Analytical Gradient of Scalar Field
  if ( nsuby > 1 ) then
     drhodx = kk * cos( kk * cx) * sin( ll * cy ) * sin( mm * cz )
     drhody = ll * sin( kk * cx) * cos( ll * cy ) * sin( mm * cz )
     drhodz = mm * sin( kk * cx) * sin( ll * cy ) * cos( mm * cz )
  else
     drhodx = kk * cos( kk * cx) * sin( mm * cz )
     drhody = 0.0_dp
     drhodz = mm * sin( kk * cx) * cos( mm * cz )
  endif

  ! Analytical transport.
  analytic = -1.0_dp * (ux * drhodx + uy * drhody + uz * drhodz )

  ! Compute errors.
  linf = pmaxval( reshape( abs( analytic - Nrho0 ), (/rpk * nsuby/) )  )

  l2 = pnorm2( reshape( abs( analytic - Nrho0 ), (/rpk * nsuby/) )  ) /&
       pnorm2( reshape( abs( analytic ), (/rpk * nsuby/) )  )

  ! Compute the maximum error over all ranks and distribute across all ranks
  call MPI_ALLREDUCE( linf, linf_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  call MPI_ALLREDUCE( l2, l2_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  write( caststr, '(A,D17.10,D17.10,D17.10)' ) 'Linf Error in transport term', linf_max
  call notify( caststr )
  write( caststr, '(A,D17.10,D17.10,D17.10)' ) 'L2 Error in transport term', l2_max
  call notify( caststr )

  success_flag = ( linf_max < tolerance )

  if (rank .eq. root) then
     open( unit=65, file='transport_error.txt')
     write(65,*) linf_max
     write(65,*) l2_max
     close(65)
  endif

end subroutine validate_transport_term
