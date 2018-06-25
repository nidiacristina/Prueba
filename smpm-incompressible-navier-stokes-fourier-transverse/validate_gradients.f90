subroutine validate_gradients( tolerance, success_flag )
! Validates the gradient function in 2D and 3D using a sinusoidal 
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

  real(kind=dp)                               :: pi, kk, ll, mm
  real(kind=dp)                               :: linf_x, linf_y, linf_z
  real(kind=dp)                               :: l2_x, l2_y, l2_z
  real(kind=dp)                               :: linf_x_max, linf_y_max, linf_z_max
  real(kind=dp)                               :: l2_x_max, l2_y_max, l2_z_max


  integer                                     :: ierr

  real(kind=dp), allocatable, dimension(:, :) :: phi

  real(kind=dp), allocatable, dimension(:, :) :: dphidx_num, dphidy_num, dphidz_num
  real(kind=dp), allocatable, dimension(:, :) :: dphidx_ana, dphidy_ana, dphidz_ana

  character(len=256)                          :: caststr

  real(kind=dp)                               :: Lx, Lz

  allocate( phi(1:rpk, 1:nsuby) )

  allocate( dphidx_num(1:rpk, 1:nsuby), dphidy_num(1:rpk,1:nsuby), dphidz_num(1:rpk,1:nsuby) )
  allocate( dphidx_ana(1:rpk, 1:nsuby), dphidy_ana(1:rpk,1:nsuby), dphidz_ana(1:rpk,1:nsuby) )

  ! Get pi.
  pi = 2.0_dp * acos( 0.0_dp )

  ! Get domain lengths.
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

  ! Compute Analytical Gradients
  if (nsuby > 1) then
     dphidx_ana = kk * cos( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     dphidy_ana = ll * sin( kk * cx ) * cos( ll * cy ) * sin( mm * cz )
     dphidz_ana = mm * sin( kk * cx ) * sin( ll * cy ) * cos( mm * cz )
  else
     dphidx_ana = kk * cos( kk * cx ) * sin( mm * cz )
     dphidy_ana = 0.0_dp
     dphidz_ana = mm * sin( kk * cx ) * cos( mm * cz )
  endif

  ! Check Numerical Gradients
  call compute_3D_gradient( phi, dphidx_num, dphidy_num, dphidz_num )

  ! Compute errors.
  linf_x = pmaxval( reshape( abs( dphidx_ana - dphidx_num ), (/rpk * nsuby/) )  )
  linf_y = pmaxval( reshape( abs( dphidy_ana - dphidy_num ), (/rpk * nsuby/) )  )
  linf_z = pmaxval( reshape( abs( dphidz_ana - dphidz_num ), (/rpk * nsuby/) )  )

  l2_x = pnorm2( reshape( abs(dphidx_ana - dphidx_num) , (/rpk * nsuby/) )  ) /&
         pnorm2( reshape( abs(dphidx_ana), (/rpk * nsuby/) )  )

  if (nsuby .gt. 1) then
     l2_y = pnorm2( reshape( abs(dphidy_ana - dphidy_num) , (/rpk * nsuby/) )  ) /&
            pnorm2( reshape( abs(dphidy_ana), (/rpk * nsuby/) )  )
  else
     l2_y = 0.0_dp
  endif

  l2_z = pnorm2( reshape( abs(dphidz_ana - dphidz_num) , (/rpk * nsuby/) )  ) /&
         pnorm2( reshape( abs(dphidz_ana), (/rpk * nsuby/) )  )

  ! Compute the maximum error over all ranks and distribute across all ranks 
  call MPI_ALLREDUCE( linf_x, linf_x_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  call MPI_ALLREDUCE( linf_y, linf_y_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  call MPI_ALLREDUCE( linf_z, linf_z_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  call MPI_ALLREDUCE( l2_x, l2_x_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  call MPI_ALLREDUCE( l2_y, l2_y_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  call MPI_ALLREDUCE( l2_z, l2_z_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  write( caststr, '(A,D17.10,D17.10,D17.10)' ) 'Error in gradients (x,y,z):', linf_x_max, linf_y_max, linf_z_max
  call notify( caststr )

  success_flag = ((linf_x_max < tolerance) .and. &
                  (linf_y_max < tolerance) .and. &
                  (linf_z_max < tolerance))

  if (rank .eq. root) then
     open(unit=65,file='gradients_error.txt')
     write(65,*) linf_x_max
     write(65,*) linf_y_max
     write(65,*) linf_z_max
     write(65,*) l2_x_max
     write(65,*) l2_y_max
     write(65,*) l2_z_max    
     close(65)
  endif

end subroutine validate_gradients
