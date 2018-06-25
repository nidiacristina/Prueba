subroutine validate_double_curl( tolerance, success_flag )
! Validates the double curl function in 2D and 3D using a sinusoidal 
! with a user-specified wavenumber. 

  use constants, only:               nsuby, rank, root
  use geom, only:                    cx, cz
  use field_variables, only:         ux, uy, uz
  use mpi, only:                     MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_MAX
  use parallel_linear_algebra
  use precision, only:               dp
  use transverse, only:              cy, Ly
  use woodbury_matrices, only:       rpk

  implicit none

  real(kind=dp), intent(in)                   :: tolerance
  logical, intent(out)                        :: success_flag

  real(kind=dp)                               :: pi, kk, ll, mm

  real(kind=dp)                               :: linf_error_x, linf_error_y, linf_error_z
  real(kind=dp)                               :: linf_error_xx, linf_error_yy, linf_error_zz

  real(kind=dp)                               :: l2_error_x, l2_error_y, l2_error_z
  real(kind=dp)                               :: l2_error_xx, l2_error_yy, l2_error_zz

  real(kind=dp)                               :: linf_error_x_max, linf_error_y_max, linf_error_z_max
  real(kind=dp)                               :: linf_error_xx_max, linf_error_yy_max, linf_error_zz_max

  real(kind=dp)                               :: l2_error_x_max, l2_error_y_max, l2_error_z_max
  real(kind=dp)                               :: l2_error_xx_max, l2_error_yy_max, l2_error_zz_max

  integer                                     :: ierr

  real(kind=dp), allocatable, dimension(:, :) :: curli_ana, curlj_ana, curlk_ana
  real(kind=dp), allocatable, dimension(:, :) :: curlii_ana, curljj_ana, curlkk_ana

  real(kind=dp), allocatable, dimension(:, :) :: curli_num, curlj_num, curlk_num
  real(kind=dp), allocatable, dimension(:, :) :: curlii_num, curljj_num, curlkk_num

  real(kind=dp), allocatable, dimension(:, :) :: grad_of_div_x, grad_of_div_y, grad_of_div_z

  real(kind=dp), allocatable, dimension(:, :) :: lap_ux, lap_uy, lap_uz

  real(kind=dp), allocatable, dimension(:, :) :: ux_x_ana, ux_y_ana, ux_z_ana
  real(kind=dp), allocatable, dimension(:, :) :: uy_x_ana, uy_y_ana, uy_z_ana
  real(kind=dp), allocatable, dimension(:, :) :: uz_x_ana, uz_y_ana, uz_z_ana

  character(len=256)                          :: caststr

  real(kind=dp)                               :: Lx, Lz

  allocate( curli_ana(1:rpk, 1:nsuby), curlj_ana(1:rpk,1:nsuby), curlk_ana(1:rpk,1:nsuby) )
  allocate( curlii_ana(1:rpk, 1:nsuby), curljj_ana(1:rpk,1:nsuby), curlkk_ana(1:rpk,1:nsuby) )

  allocate( curli_num(1:rpk, 1:nsuby), curlj_num(1:rpk,1:nsuby), curlk_num(1:rpk,1:nsuby) )
  allocate( curlii_num(1:rpk, 1:nsuby), curljj_num(1:rpk,1:nsuby), curlkk_num(1:rpk,1:nsuby) )

  allocate( grad_of_div_x(1:rpk, 1:nsuby), grad_of_div_y(1:rpk,1:nsuby), grad_of_div_z(1:rpk,1:nsuby) )

  allocate( lap_ux(1:rpk, 1:nsuby), lap_uy(1:rpk,1:nsuby), lap_uz(1:rpk,1:nsuby) )

  allocate( ux_x_ana(1:rpk, 1:nsuby), ux_y_ana(1:rpk,1:nsuby), ux_z_ana(1:rpk,1:nsuby) )
  allocate( uy_x_ana(1:rpk, 1:nsuby), uy_y_ana(1:rpk,1:nsuby), uy_z_ana(1:rpk,1:nsuby) )
  allocate( uz_x_ana(1:rpk, 1:nsuby), uz_y_ana(1:rpk,1:nsuby), uz_z_ana(1:rpk,1:nsuby) )

  ! Get pi.
  pi = 2.0_dp * acos( 0.0_dp )

  ! Compute domain lengths
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
     ux = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     uy = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     uz = sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
  else
     ux = sin( kk * cx ) * sin( mm * cz )
     uy = 0.0_dp
     uz = sin( kk * cx ) * sin( mm * cz )
  endif

  ! Compute analytical gradients
  if (nsuby > 1) then
     ux_x_ana = kk * cos( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     ux_y_ana = ll * sin( kk * cx ) * cos( ll * cy ) * sin( mm * cz )
     ux_z_ana = mm * sin( kk * cx ) * sin( ll * cy ) * cos( mm * cz )

     uy_x_ana = kk * cos( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     uy_y_ana = ll * sin( kk * cx ) * cos( ll * cy ) * sin( mm * cz )
     uy_z_ana = mm * sin( kk * cx ) * sin( ll * cy ) * cos( mm * cz )

     uz_x_ana = kk * cos( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     uz_y_ana = ll * sin( kk * cx ) * cos( ll * cy ) * sin( mm * cz )
     uz_z_ana = mm * sin( kk * cx ) * sin( ll * cy ) * cos( mm * cz )

  else

     ux_x_ana = kk * cos( kk * cx ) * sin( mm * cz )
     ux_y_ana = 0.0_dp
     ux_z_ana = mm * sin( kk * cx ) * cos( mm * cz )

     uy_x_ana = 0.0_dp
     uy_y_ana = 0.0_dp
     uy_z_ana = 0.0_dp

     uz_x_ana = kk * cos( kk * cx ) * sin( mm * cz )
     uz_y_ana = 0.0_dp
     uz_z_ana = mm * sin( kk * cx ) * cos( mm * cz )

  endif

  ! Compute the analytical curl
  if (nsuby > 1) then
     curli_ana =  ( uz_y_ana - uy_z_ana )
     curlj_ana =  ( ux_z_ana - uz_x_ana )
     curlk_ana =  ( uy_x_ana - ux_y_ana )
  else
     curli_ana = 0.0_dp
     curlj_ana = - ( uz_x_ana - ux_z_ana )
     curlk_ana = 0.0_dp
  endif

  ! Compute numerical curl
  call compute_3D_curl( ux, uy, uz, curli_num, curlj_num, curlk_num )

  ! Compute Error of curl
  linf_error_x = pmaxval( reshape( abs( curli_ana - curli_num ), (/rpk * nsuby/) )  )
  linf_error_y = pmaxval( reshape( abs( curlj_ana - curlj_num ), (/rpk * nsuby/) )  )
  linf_error_z = pmaxval( reshape( abs( curlk_ana - curlk_num ), (/rpk * nsuby/) )  )

  ! Compute Error of curl
  if ( nsuby > 1) then
     l2_error_x = pnorm2( reshape( abs(curli_ana - curli_num) , (/rpk * nsuby/) )  ) /&
                  pnorm2( reshape( abs(curli_ana) , (/rpk * nsuby/) )  )
     l2_error_z = pnorm2( reshape( abs(curlk_ana - curlk_num) , (/rpk * nsuby/) )  ) /&
                  pnorm2( reshape( abs(curlk_ana) , (/rpk * nsuby/) )  )
  else
     l2_error_x = 0.0_dp
     l2_error_z = 0.0_dp
  endif

  l2_error_y = pnorm2( reshape( curlj_ana - curlj_num , (/rpk * nsuby/) )  ) /&
               pnorm2( reshape( curlj_ana , (/rpk * nsuby/) )  )



  ! Compute the maximum error over all ranks and distribute across all ranks 
  call MPI_ALLREDUCE( linf_error_x, linf_error_x_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  call MPI_ALLREDUCE( linf_error_y, linf_error_y_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  call MPI_ALLREDUCE( linf_error_z, linf_error_z_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  call MPI_ALLREDUCE( l2_error_x, l2_error_x_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  call MPI_ALLREDUCE( l2_error_y, l2_error_y_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  call MPI_ALLREDUCE( l2_error_z, l2_error_z_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )


  write( caststr, '(A,D17.10,D17.10,D17.10)' ) 'Error in curl:', linf_error_x_max, linf_error_y_max, linf_error_z_max
  call notify( caststr )
 
  ! Compute double curl via the identity
  !
  !     curl( curl( u ) ) = gradient( divergence( u ) ) - laplacian( u )
  !
  ! where u is a vector field with components (ux,uy,uz).


  if (nsuby > 1) then
     ! Compute the double curl in x:
     grad_of_div_x = -kk*kk * sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) &
                   +  kk*ll * cos( kk * cx ) * cos( ll * cy ) * sin( mm * cz ) &
                   +  kk*mm * cos( kk * cx ) * sin( ll * cy ) * cos( mm * cz ) 
                             
     lap_ux        =  -(kk**2  + ll**2  + mm**2) * sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     curlii_ana    = grad_of_div_x - lap_ux

     ! Compute the double curl in y:
     grad_of_div_y =  ll*kk * cos( kk * cx ) * cos( ll * cy ) * sin( mm * cz ) &
                   -  ll*ll * sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) &
                   +  ll*mm * sin( kk * cx ) * cos( ll * cy ) * cos( mm * cz ) 
                             
     lap_uy        =  -(kk**2  + ll**2  + mm**2) * sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     curljj_ana    = grad_of_div_y - lap_uy

     ! Compute the double curl in z:
     grad_of_div_z =  mm*kk * cos( kk * cx ) * sin( ll * cy ) * cos( mm * cz ) &
                   +  mm*ll * sin( kk * cx ) * cos( ll * cy ) * cos( mm * cz ) &
                   -  mm*mm * sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz ) 
                             
     lap_uz        =  -(kk**2  + ll**2  + mm**2) * sin( kk * cx ) * sin( ll * cy ) * sin( mm * cz )
     curlkk_ana    = grad_of_div_z - lap_uz

  else

     ! Compute the double curl in x:
     grad_of_div_x = -kk*kk * sin( kk * cx ) * sin( mm * cz ) &
                   +  kk*mm * cos( kk * cx ) * cos( mm * cz ) 
                             
     lap_ux        =  -(kk**2 + mm**2) * sin( kk * cx ) * sin( mm * cz )
     curlii_ana    = grad_of_div_x - lap_ux

     ! Compute the double curl in y:
     grad_of_div_y =  0.0_dp 
                             
     lap_uy        =  0.0_dp
     curljj_ana    = grad_of_div_y - lap_uy

     ! Compute the double curl in z:
     grad_of_div_z =  mm*kk * cos( kk * cx ) * cos( mm * cz ) &
                   -  mm*mm * sin( kk * cx ) * sin( mm * cz ) 
                             
     lap_uz        =  -(kk**2 + mm**2) * sin( kk * cx ) * sin( mm * cz )
     curlkk_ana    = grad_of_div_z - lap_uz

  endif


  ! Compute the double curl numerically
  call compute_double_curl( ux, uy, uz, curlii_num, curljj_num, curlkk_num)

  ! Compute Error of curl
  linf_error_xx = pmaxval( reshape( abs( curlii_ana - curlii_num ), (/rpk * nsuby/) )  )
  linf_error_yy = pmaxval( reshape( abs( curljj_ana - curljj_num ), (/rpk * nsuby/) )  )
  linf_error_zz = pmaxval( reshape( abs( curlkk_ana - curlkk_num ), (/rpk * nsuby/) )  )

  l2_error_xx = pnorm2( reshape( abs( curlii_ana - curlii_num ), (/rpk * nsuby/) )  ) / &
                pnorm2( reshape( curlii_ana , (/rpk * nsuby/) )  )

  if (nsuby .gt. 1) then
     l2_error_yy = pnorm2( reshape( abs( curljj_ana - curljj_num ), (/rpk * nsuby/) )  ) / &
                  pnorm2( reshape( curljj_ana, (/rpk * nsuby/) )  )
  else
     l2_error_yy = 0.0_dp
  endif
  l2_error_zz = pnorm2( reshape( abs( curlkk_ana - curlkk_num ), (/rpk * nsuby/) )  ) / &
                pnorm2( reshape( curlkk_ana , (/rpk * nsuby/) )  )  

  ! Compute the maximum error over all ranks and distribute across all ranks 
  call MPI_ALLREDUCE( linf_error_xx, linf_error_xx_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  call MPI_ALLREDUCE( linf_error_yy, linf_error_yy_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  call MPI_ALLREDUCE( linf_error_zz, linf_error_zz_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  call MPI_ALLREDUCE( l2_error_xx, l2_error_xx_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  call MPI_ALLREDUCE( l2_error_yy, l2_error_yy_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  call MPI_ALLREDUCE( l2_error_zz, l2_error_zz_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )


  write( caststr, '(A,D17.10,D17.10,D17.10)' ) 'Error in double curl:', linf_error_xx_max, linf_error_yy_max, linf_error_zz_max
  call notify( caststr )

   success_flag = ((linf_error_x_max  < tolerance) .and. &
                   (linf_error_y_max  < tolerance) .and. &
                   (linf_error_z_max  < tolerance) .and. &
                   (linf_error_xx_max < tolerance) .and. &
                   (linf_error_yy_max < tolerance) .and. &
                   (linf_error_zz_max < tolerance))
  
  if (rank .eq. root) then
     open(unit=65,file='curl_error.txt')
     write(65,*) linf_error_x_max
     write(65,*) linf_error_y_max
     write(65,*) linf_error_z_max
     write(65,*) linf_error_xx_max
     write(65,*) linf_error_yy_max
     write(65,*) linf_error_zz_max
     write(65,*) l2_error_x_max
     write(65,*) l2_error_y_max
     write(65,*) l2_error_z_max
     write(65,*) l2_error_xx_max
     write(65,*) l2_error_yy_max
     write(65,*) l2_error_zz_max 
     close(65)
  endif
end subroutine validate_double_curl
