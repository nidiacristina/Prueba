subroutine apply_extension( q, Bq, Nrhs )
! Compute the extension operator that takes a vector from the capacitance grid
! to the whole grid by an inclusion.
!
! q  - the input vector.
! Bq - the output vector.

  use constants, only:         nsubx, rank
  use precision, only:         dp
  use woodbury_matrices, only: aCk, dimA, dimCk, numA_per_rank, rpk, s, zCk

  implicit none

  real(kind=dp), dimension(dimCk, Nrhs), intent(in)   :: q
  real(kind=dp), dimension(rpk,  Nrhs), intent(inout) :: Bq
  integer, intent(in)                                 :: Nrhs
  integer                                             :: ii
  integer                                             :: a, z, block_number, aC, zC, first_block

  ! Set the output vector to zero (there is no good reason it should ever not
  ! be preinitialized to zero).
  Bq = 0.0_dp

  ! Loop over the interfaces and put the array in its place in the full grid.
  first_block = rank * numA_per_rank + 1
  do ii = 1, numA_per_rank

     ! Get the ownership range of this block in the full grid on this rank.
     a = (ii - 1) * dimA + 1
     z = (ii - 1) * dimA + dimA

     ! Get the block number.
     block_number = rank * numA_per_rank  + ii

     ! Get the local extent of this block.
     aC = aCk(block_number) - aCk(first_block) + 1
     zC = aC + zCk(block_number) - aCk(block_number)

     ! Apply extension.
     if ( block_number > 1 ) then
       Bq(a:a + s - 1, :) = q(aC:aC + s - 1, :)
     endif
     if ( block_number < nsubx ) then
        Bq(z - s + 1:z, :) = q(zC - s + 1:zC, :)
     endif

  enddo

end subroutine apply_extension
