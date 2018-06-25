subroutine validate_advection_terms( tolerance, success_flag )
! Validates the advective terms in X, Y, and Z.

  use constants, only:               nsuby, rank, root
  use field_variables, only:         Nux0, Nuy0, Nuz0, ux, ux0, uy, uy0, uz, uz0, &
                                     dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
  use geom, only:                    cx, cz
  use mpi, only:                     MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_MAX
  use parallel_linear_algebra
  use precision, only:               dp
  use transverse, only:              cy, Ly
  use woodbury_matrices, only:       rpk

  implicit none

  real(kind=dp), intent(in)                   :: tolerance
  logical, intent(out)                        :: success_flag

  real(kind=dp)                               :: pi, kk, ll, mm
  real(kind=dp)                               :: linf_x, linf_y, linf_z
  real(kind=dp)                               :: linf_x_max, linf_y_max, linf_z_max
  real(kind=dp)                               :: l2_x, l2_y, l2_z
  real(kind=dp)                               :: l2_x_max, l2_y_max, l2_z_max

  real(kind=dp), allocatable, dimension(:, :) :: analytic_x, analytic_y, analytic_z

  real(kind=dp), allocatable, dimension(:, :) :: duudx, duvdy, duwdz
  real(kind=dp), allocatable, dimension(:, :) :: dvudx, dvvdy, dvwdz
  real(kind=dp), allocatable, dimension(:, :) :: dwudx, dwvdy, dwwdz

  real(kind=dp)                               :: Lx, Lz

  character(len=256)                          :: caststr
  integer                                     :: ierr

  allocate( analytic_x(1:rpk, 1:nsuby), analytic_y(1:rpk, 1:nsuby), analytic_z(1:rpk, 1:nsuby) )

  allocate( duudx(1:rpk, 1:nsuby), duvdy(1:rpk,1:nsuby), duwdz(1:rpk,1:nsuby) )
  allocate( dvudx(1:rpk, 1:nsuby), dvvdy(1:rpk,1:nsuby), dvwdz(1:rpk,1:nsuby) )
  allocate( dwudx(1:rpk, 1:nsuby), dwvdy(1:rpk,1:nsuby), dwwdz(1:rpk,1:nsuby) )

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

  ! Set the function values.
  if (nsuby > 1) then
     ux0 = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     uy0 = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     uz0 = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
  else
     ux0 = sin( kk * cx ) * sin( mm * cz )
     uy0 = 0.0_dp
     uz0 = sin( kk * cx ) * sin( mm * cz )
  endif
  ! Save the analytic values so we can compare against them later.
  ux = ux0
  uy = uy0
  uz = uz0

  ! Compute the advective terms.
  if (nsuby > 1) then
     call apply_smpm_advection( Nux0, ux0, ux0, uy0, uz0 )
     call apply_smpm_advection( Nuy0, uy0, ux0, uy0, uz0 )
     call apply_smpm_advection( Nuz0, uz0, ux0, uy0, uz0 )
  else
     call apply_smpm_advection( Nux0, ux0, ux0, uy0, uz0 )
     Nuy0 = 0.0_dp
     call apply_smpm_advection( Nuz0, uz0, ux0, uy0, uz0 )
  endif

  ! Analytical gradients in x-momentum.
  if (nsuby > 1) then
     analytic_x = kk * cos( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     analytic_y = ll * sin( kk * cx ) * cos( ll * cy ) * sin( mm * cz )
     analytic_z = mm * sin( kk * cx ) * sin( ll * cy ) * cos( mm * cz )
  else
     analytic_x = kk * cos( kk * cx ) * sin( mm * cz )
     analytic_y = 0.0_dp
     analytic_z = mm * sin( kk * cx ) * cos( mm * cz )
  endif

  ! Check error on derivatives in x-direction.
  call compute_3D_gradient( ux0, dudx, dudy, dudz )

  ! Analytical gradients in y-momentum.
  if (nsuby > 1) then
     analytic_x = kk * cos( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     analytic_y = ll * sin( kk * cx ) * cos( ll * cy ) * sin( mm * cz )
     analytic_z = mm * sin( kk * cx ) * sin( ll * cy ) * cos( mm * cz )
  else
     analytic_x = 0.0_dp
     analytic_y = 0.0_dp
     analytic_z = 0.0_dp
  endif

  if (nsuby > 1) then
     ! Check error on derivatives in y-direction.
     call compute_3D_gradient( uy0, dvdx, dvdy, dvdz )
  else
     dvdx = 0.0_dp
     dvdy = 0.0_dp
     dvdz = 0.0_dp
  endif

  ! Analytical gradients in x-momentum.
  if (nsuby > 1) then
     analytic_x = kk * cos( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     analytic_y = ll * sin( kk * cx ) * cos( ll * cy ) * sin( mm * cz )
     analytic_z = mm * sin( kk * cx ) * sin( ll * cy ) * cos( mm * cz )
  else
     analytic_x = kk * cos( kk * cx ) * sin( mm * cz )
     analytic_y = 0.0_dp
     analytic_z = mm * sin( kk * cx ) * cos( mm * cz )
  endif


  ! Analytical skew-symmetric term.
  if (nsuby > 1) then
     duudx = kk * sin( 2.0_dp * kk * cx ) * sin( ll * cy )**2.0_dp  * sin( mm * cz )**2.0_dp
     duvdy = ll * sin( kk * cx )**2.0_dp  * sin( 2.0_dp * ll * cy ) * sin( mm * cz )**2.0_dp
     duwdz = mm * sin( kk * cx )**2.0_dp  * sin( ll * cy )**2.0_dp  * sin( 2.0_dp * mm * cz )

     dvudx = kk * sin( 2.0_dp * kk * cx ) * sin( ll * cy )**2.0_dp  * sin( mm * cz )**2.0_dp
     dvvdy = ll * sin( kk * cx )**2.0_dp  * sin( 2.0_dp * ll * cy ) * sin( mm * cz )**2.0_dp
     dvwdz = mm * sin( kk * cx )**2.0_dp  * sin( ll * cy )**2.0_dp  * sin( 2.0_dp * mm * cz )

     dwudx = kk * sin( 2.0_dp * kk * cx ) * sin( ll * cy )**2.0_dp  * sin( mm * cz )**2.0_dp
     dwvdy = ll * sin( kk * cx )**2.0_dp  * sin( 2.0_dp * ll * cy ) * sin( mm * cz )**2.0_dp
     dwwdz = mm * sin( kk * cx )**2.0_dp  * sin( ll * cy )**2.0_dp  * sin( 2.0_dp * mm * cz )

  else
     duudx = kk * sin( 2.0_dp * kk * cx ) * sin( mm * cz )**2.0_dp
     duvdy = 0.0_dp
     duwdz = mm * sin( kk * cx )**2.0_dp  * sin( 2.0_dp * mm * cz )

     dvudx = 0.0_dp 
     dvvdy = 0.0_dp 
     dvwdz = 0.0_dp 

     dwudx = kk * sin( 2.0_dp * kk * cx ) * sin( mm * cz )**2.0_dp
     dwvdy = 0.0_dp 
     dwwdz = mm * sin( kk * cx )**2.0_dp  * sin( 2.0_dp * mm * cz )

  endif

  ! Check error on derivatives in z-direction.
  call compute_3D_gradient( uz0, dwdx, dwdy, dwdz )


  ! Analytical advection in X, Y, and Z.
  if (nsuby > 1) then
     analytic_x = -0.5_dp * (ux * dudx + uy * dudy + uz * dudz + (duudx + duvdy + duwdz))
     analytic_y = -0.5_dp * (ux * dvdx + uy * dvdy + uz * dvdz + (dvudx + dvvdy + dvwdz))
     analytic_z = -0.5_dp * (ux * dwdx + uy * dwdy + uz * dwdz + (dwudx + dwvdy + dwwdz))
  else
     analytic_x = -0.5_dp * (ux * dudx + uz * dudz + (duudx + duwdz))
     analytic_y =  0.0_dp
     analytic_z = -0.5_dp * (ux * dwdx + uz * dwdz + (dwudx + dwwdz))
  endif


  ! Compute errors.
  linf_x = pmaxval( reshape( abs( analytic_x - Nux0 ), (/rpk * nsuby/) )  )
  linf_y = pmaxval( reshape( abs( analytic_y - Nuy0 ), (/rpk * nsuby/) )  )
  linf_z = pmaxval( reshape( abs( analytic_z - Nuz0 ), (/rpk * nsuby/) )  )

  l2_x = pnorm2( reshape( abs( analytic_x - Nux0 ), (/rpk * nsuby/) )  ) /&
         pnorm2( reshape( abs(analytic_x), (/rpk * nsuby/) )  )

  if (nsuby .gt. 1) then
        l2_y = pnorm2( reshape( abs( analytic_y - Nuy0 ), (/rpk * nsuby/) )  ) /&
               pnorm2( reshape( abs( analytic_y) , (/rpk * nsuby/) )  )
  else
        l2_y = 0.0_dp
  endif
  l2_z = pnorm2( reshape( abs( analytic_z - Nuz0 ), (/rpk * nsuby/) )  ) /&
         pnorm2( reshape( abs(analytic_z), (/rpk * nsuby/) )  )

  ! Compute the maximum error over all ranks and distribute across all ranks
  call MPI_ALLREDUCE( linf_x, linf_x_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( linf_y, linf_y_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( linf_z, linf_z_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  call MPI_ALLREDUCE( l2_x, l2_x_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( l2_y, l2_y_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( l2_z, l2_z_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  write( caststr, '(A,D17.10,D17.10,D17.10)' ) 'Error in advective terms (x,y,z):', &
                                               linf_x_max, linf_y_max, linf_z_max
  call notify( caststr )

  success_flag = ((linf_x_max < tolerance) .and. &
                  (linf_y_max < tolerance) .and. &
                  (linf_z_max < tolerance))

  if (rank .eq. root) then
     open( unit=65, file='skew_symmetric_errors.txt')
     write(65,*) linf_x_max
     write(65,*) linf_y_max
     write(65,*) linf_z_max
     write(65,*) l2_x_max
     write(65,*) l2_y_max
     write(65,*) l2_z_max
     close(65)
  endif

end subroutine validate_advection_terms
