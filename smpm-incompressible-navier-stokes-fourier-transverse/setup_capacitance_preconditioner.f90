subroutine setup_capacitance_preconditioner
! This implements a block-Jacobi preconditioner for the capacitance problem.
! This is the setup routine for the preconditioner.

  use constants, only:         nsubx, rank
  use precision, only:         dp
  use woodbury_matrices, only: C_poisson_block, numA_per_rank, preC, preC_pivot, s

  implicit none

  ! Internal variables.
  integer                            :: ii, jj
  integer                            :: info
  integer                            :: block_number

  ! Flux communication across ranks.
  real(kind=dp), dimension(1:s, 1:s) :: to_left, to_right, from_left, from_right

  ! Pass information between ranks that we'll need to assemble the factors.
  to_left  = C_poisson_block(1:s, 1:s, 1)
  to_right = C_poisson_block(4*s - s + 1:4*s, 2*s - s + 1:2*s, numA_per_rank)
  call sync_ranks( to_left, to_right, s**2 )

  ! Loop over this processes' capacitance blocks, factoring and storing the
  ! factors.
  do ii = 1, numA_per_rank

     ! Get this block's # in the global block sequence.
     block_number = numA_per_rank * rank + ii

     ! Set up the left/right interblock comms.
     if ( ii == 1 ) then
        from_left = to_left
     else
        from_left = C_poisson_block(3*s + 1:4*s, s+1:2*s, ii - 1)
     endif
     if ( ii == numA_per_rank ) then
        from_right = to_right
     else
        from_right = C_poisson_block(1:s, 1:s, ii + 1)
     endif

     ! Get the indices for building the factors.
     if ( block_number == 1 ) then
        ! Get the portion that's in this block.
        preC(2*s + 1:4*s, 2*s + 1:3*s, ii) = C_poisson_block(2*s + 1:4*s, s + 1:2*s, ii)

        ! Apply addition from the right.
        preC(2*s + 1:3*s, 3*s + 1:4*s, ii) = from_right

        ! Add the identity.
        do jj = 2 * s + 1, 4 * s
           preC(jj, jj, ii) = preC(jj, jj, ii) + 1.0_dp
        enddo

        ! Factor this block.
        call DGETRF( 2*s, 2*s, preC(2*s + 1, 2*s + 1, ii), 4*s, &
                     preC_pivot(2*s + 1, ii), info )
     elseif ( block_number == nsubx ) then
        ! Get the portion that's in this block.
        preC(1:2*s, s+1:2*s, ii) = C_poisson_block(1:2*s, 1:s, ii)

        ! Apply the addition from the left.
        preC(s+1:2*s, 1:s, ii)   = from_left

        ! Add the identity.
        do jj = 1, 2 * s
           preC(jj, jj, ii) = preC(jj, jj, ii) + 1.0_dp
        enddo

        ! Factor this block.
        call DGETRF( 2*s, 2*s, preC(1, 1, ii), 4*s, &
                     preC_pivot(1, ii), info )
     else
        ! Get the porition that's in this block.
        preC(1:4*s, s+1:3*s, ii) = C_poisson_block(:, :, ii)

        ! Apply the addition from the left.
        preC(s+1:2*s, 1:s, ii)   = from_left

        ! Apply the addition from the right.
        preC(2*s + 1: 3*s, 3*s + 1:4*s, ii) = from_right

        ! Add the identity.
        do jj = 1, 4 * s
           preC(jj, jj, ii) = preC(jj, jj, ii) + 1.0_dp
        enddo

        ! Factor this block.
        call DGETRF( 4*s, 4*s, preC(1, 1, ii), 4*s, &
                     preC_pivot(1, ii), info )
      endif

  enddo

end subroutine setup_capacitance_preconditioner
