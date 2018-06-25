subroutine apply_smpm_poisson_transpose_parallel( Tx, x )
! Applies the transpose of the Spectral Multidomain Penalty method operator
! for the pressure Poisson equation to a vector x and returns the
! matrix-vector product Ax.

  use constants, only:               nsg
  use parallel_linear_algebra, only: pdot_product
  use precision, only:               dp
  use woodbury_matrices, only:       arank, rpk, zrank

  implicit none

  real(kind=dp), dimension(1:rpk), intent(out) :: Tx
  real(kind=dp), dimension(1:rpk), intent(in)  :: x

  real(kind=dp), dimension(1:rpk)              :: e_local
  real(kind=dp), dimension(1:nsg)              :: e_global
  integer                                      :: ii, local_ndx
  real(kind=dp), dimension(1:rpk)              :: iiTx
  real(kind=dp)                                :: LTx

  ! NOTE: You could do this faster by using a global basis vector that looks
  !       like e_i = [i:2*dimA:nsg], since the block-diagonal structure of the
  !       matrix allows two-columns to be computed at once so long as the
  !       columns are in non-adjacent blocks.

  do ii = 1, nsg
     ! Set the global Euclidean basis vector.
     e_global     = 0.0_dp
     e_global(ii) = 1.0_dp

     ! Get the local basis vector.
     e_local = e_global(arank:zrank)

     ! Apply this column.
     iiTx = 0.0_dp
     call apply_smpm_poisson( iiTx, e_local )

     ! Compute the dot product.
     LTx = pdot_product( iiTx, x )

     ! Store the dot product in the correct row.
     if ( ii .ge. arank .and. ii .le. zrank ) then
        local_ndx     = ii - arank + 1
        Tx(local_ndx) = LTx
     endif
  enddo

end subroutine apply_smpm_poisson_transpose_parallel
