! This file contains a bunch of matvecs related to the deflation
! preconditioner.

subroutine apply_Z( q, Zq, Nrhs )
!
! Applies Z to q and returns the product.  Z is the deflation matrix.

  use constants,         only:     nsubx, rank
  use precision,         only:     dp
  use woodbury_matrices, only:     aCk, arankC, dimCk, numA_per_rank, s, zCk

  implicit none

  integer, intent(in)                                      :: Nrhs
  real(kind=dp), dimension(1:nsubx-1, 1:Nrhs), intent(in)  ::  q
  real(kind=dp), dimension(1:dimCk, 1:Nrhs),   intent(out) :: Zq

  ! Internal variables
  integer                                                  :: ii, jj, a, z, block_number
  integer                                                  :: interface_L, interface_R

  ! Set the output to zero.
  Zq = 0.0_dp

  ! Figure out the interfaces this rank has ownership of.
  do ii = 1, numA_per_rank

     block_number = numA_per_rank * rank + ii
     a            = aCk( block_number ) - arankC + 1
     z            = zCk( block_number ) - arankC + 1! These are indices in [1,k].

     ! Which column(s)/inteface in C are we trying to extract?
     interface_L = block_number - 1  !
     interface_R = block_number      ! "left" and "right" refer to the left/right sides of this block.

     ! Grab the pieces that are needed for this block.
     if ( block_number == 1 ) then
        do jj = 1, Nrhs
           Zq( 1:s, jj ) = q( interface_R, jj )
        enddo
     elseif ( block_number == nsubx ) then
        do jj = 1, Nrhs
           Zq( dimCk - s + 1: dimCk, jj ) = q( interface_L, jj )
        enddo
     else
        do jj = 1, Nrhs
           Zq( a:a+s-1, jj ) = q( interface_L, jj )
           Zq( a + s:z, jj ) = q( interface_R, jj )
        enddo
     endif
  enddo


end subroutine apply_Z



subroutine apply_Z_transpose( q, ZTq, Nrhs )
!
! Applies Z to q and returns the product.  Z is the deflation matrix.

  use constants,         only:     nsubx, rank
  use mpi,               only:     MPI_COMM_WORLD, MPI_SUM, MPI_DOUBLE_PRECISION, MPI_IN_PLACE
  use precision,         only:     dp
  use woodbury_matrices, only:     aCk, arankC, dimCk, numA_per_rank, s, zCk

  implicit none

  integer, intent(in)                                      :: Nrhs
  real(kind=dp), dimension(1:dimCk, 1:Nrhs),   intent(in)  ::  q
  real(kind=dp), dimension(1:nsubx-1, 1:Nrhs), intent(out) :: ZTq

  ! Internal variables
  integer                                                  :: ii, a, z, block_number
  integer                                                  :: interface_L, interface_R
  integer                                                  :: ierr

  ! Set the output to zero.
  ZTq = 0.0_dp

  ! Figure out the interfaces this rank has ownership of.
  do ii = 1, numA_per_rank

     block_number = numA_per_rank * rank + ii
     a            = aCk( block_number ) - arankC + 1
     z            = zCk( block_number ) - arankC + 1! These are indices in [1,k].

     ! Which column(s)/inteface in C are we trying to extract?
     interface_L = block_number - 1  !
     interface_R = block_number      ! "left" and "right" refer to the left/right sides of this block.

     ! Sum over the interfaces.
     if ( block_number == 1 ) then
        ZTq( interface_R, 1:Nrhs ) = ZTq( interface_R, 1:Nrhs ) + sum( q(1:s, 1:Nrhs), 1 )
     elseif ( block_number == nsubx ) then
        ZTq( interface_L, 1:Nrhs ) = ZTq( interface_L, 1:Nrhs ) + sum( q( dimCk - s + 1: dimCk, 1:Nrhs ), 1 )
     else
        ZTq( interface_L, 1:Nrhs ) = ZTq( interface_L, 1:Nrhs ) + sum( q( a:a + s - 1, 1:Nrhs ), 1 )
        ZTq( interface_R, 1:Nrhs ) = ZTq( interface_R, 1:Nrhs ) + sum( q( a + s:z, 1:Nrhs ), 1 )
     endif

  enddo

  ! Sum up over all ranks.
  call MPI_ALLREDUCE( MPI_IN_PLACE, ZTq, (nsubx - 1)*Nrhs, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )

end subroutine apply_Z_transpose


subroutine apply_3D_P( Pq, q )
!
! Applies P to q and returns the product.  P is the deflation left projection matrix, one for
! each transverse wavenumber.

  use constants,         only:     nky, nsubx 
  use precision,         only:     dp
  use woodbury_matrices, only:     S_coarse, dimCk, uC_coarse
 
  implicit none

  real(kind=dp), dimension(1:dimCk, 1:nky), intent(in)  :: q
  real(kind=dp), dimension(1:dimCk, 1:nky), intent(out) :: Pq

  ! Internal variables
  integer                                                  :: ii
  real(kind=dp), dimension( 1:nsubx - 1, 1:nky )        :: ZTq, coarse_rhs
  real(kind=dp), dimension( 1:dimCk, 1:nky )            :: Zq

  real(kind=dp), external                               :: DDOT

  ! Apply the projection onto the coarse grid.
  call apply_Z_transpose( q, ZTq, nky )

  ! Regularize the zero wavenumber.
  ZTq(:, 1) = ZTq( :, 1 ) - uC_coarse * DDOT( size(uC_coarse, 1), uC_coarse, 1, ZTq( :, 1), 1 )

  ! Solve the coarse systems, one for each wavenumber.
  coarse_rhs = ZTq
  do ii = 1, nky
     call solve_tridiagonal_system( nsubx - 1, S_coarse(:, :, ii), ZTq(:, ii ) )
  enddo

  ! Multiply by the inflation.
  call apply_Z( ZTq, Zq, nky )

  ! Apply the Schur matrix.
  call apply_3D_poisson_schur_sparse( Pq, Zq )

  ! Add the identity.
  Pq = q - Pq

end subroutine apply_3D_P


subroutine apply_3D_Q( Qq, q )
!
! Applies Q to q and returns the product.  Q is the deflation right projection matrix.

  use constants,         only:     nky, nsubx
  use precision,         only:     dp
  use woodbury_matrices, only:     S_coarse, dimCk, uC_coarse

  implicit none

  real(kind=dp), dimension(1:dimCk, 1:nky), intent(in)  :: q
  real(kind=dp), dimension(1:dimCk, 1:nky), intent(out) :: Qq

  ! Internal variables
  integer                                               :: ii
  real(kind=dp), dimension( 1:nsubx - 1, 1:nky )        :: ZTq
  real(kind=dp), dimension( 1:dimCk, 1:nky )            :: Sq

  real(kind=dp), external                               :: DDOT

  ! Apply the Schur matrix.
  Sq = 0.0_dp
  call apply_3D_poisson_schur_sparse( Sq, q )

  ! Apply the projection onto the coarse grid.
  call apply_Z_transpose( Sq, ZTq, nky )

  ! Regularize the zero wavenumber.
  ZTq(:, 1) = ZTq( :, 1 ) - uC_coarse * DDOT( size( uC_coarse, 1 ), uC_coarse, 1, ZTq( :, 1), 1 )

  ! Solve the coarse systems, one for each wavenumber.
  do ii = 1, nky
     call solve_tridiagonal_system( nsubx - 1, S_coarse( :, :, ii ), ZTq( :, ii ) )
  enddo

  ! Multiply by the inflation.
  call apply_Z( ZTq, Qq, nky )

  ! Add the identity.
  Qq = q - Qq

end subroutine apply_3D_Q




