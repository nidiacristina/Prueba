subroutine set_ownership_ranges()
! Computes ownership ranges (indices) as well as MPI displacement and receive
! counts for each rank.

  use constants,         only:       n, nsubx, nsubz, nprocs, rank
  use woodbury_matrices, only:       aCk, arank, arankC, &
                                     dimA, dimC, dimCk, &
                                     displacements, k, numA_per_rank, &
                                     recv_counts, s, zCk, zrank, zrankC

  implicit none

  integer :: ii

  ! Set the ownership range of this rank within the full grid.
  arank = dimA * numA_per_rank * rank + 1
  zrank = dimA * numA_per_rank * rank + dimA * numA_per_rank

  ! Set the dimension of the local capacitance grid this rank is responsible for.
  dimC = 2 * n * nsubz * (nsubx - 1)
  allocate( displacements(1:nprocs), recv_counts(1:nprocs) )

  if ( nprocs > 1 ) then
     if ( rank == 0 .or. rank == nprocs - 1 ) then
        ! Give this block one interface worth of the k-grid plus any other
        ! remaining blocks it has ownership of.
        dimCk  = s + 2 * s * (numA_per_rank - 1)
     else
        ! Give this block two interfaces worth of the k-grid for every block
        ! it owns.
        dimCk = 2 * s * numA_per_rank
     endif
  else
     ! If we're only running on one rank, this rank gets ownership of the
     ! whole grid.
     dimCk = k
  endif

  ! Set up the arrays for MPI_GATHERV()/MPI_SCATTERV() for operations on the
  ! split k-grid.
  displacements(1) = 0

  if ( nprocs > 1 ) then
     ! The first and last ranks hold one fewer interface each (left and right,
     ! respectively, are missing), while everything in the middle holds two
     ! interfaces per domain.
     recv_counts         = 2 * s * numA_per_rank
     recv_counts(1)      = s + 2 * s * (numA_per_rank - 1)
     recv_counts(nprocs) = s + 2 * s * (numA_per_rank - 1)

     do ii = 2, nprocs
        displacements(ii) = displacements(ii - 1) + recv_counts(ii - 1)
     enddo
  else
     ! A single rank holds the entire k-grid.
     recv_counts(1) = k
  end if

  ! Set the ownership range of the capacitance grid for this rank.
  arankC = displacements(rank + 1) + 1
  zrankC = arankC + dimCk - 1

  if ( nprocs == 1 ) then
     arankC = 1
     zrankC = k
  endif

  ! Set the ownership range of each of the capacitance blocks for the whole
  ! grid.
  allocate( aCk(1:nsubx), zCk(1:nsubx) )
  aCk(1) = 1
  zCk(1) = s

  do ii = 2, nsubx - 1
     aCk(ii) = zCk(ii - 1) + 1
     zCk(ii) = zCk(ii - 1) + 2 * s
  enddo

  aCk(nsubx) = k - s + 1
  zCk(nsubx) = k

end subroutine set_ownership_ranges
