subroutine setup_3D_schur_preconditioner
! This implements a block-Jacobi preconditioner for all of the schur problems.
! This is the setup routine for the preconditioner.

  use constants, only:             nky, rank, nsubx
  use precision, only:             dp
  use woodbury_matrices, only:     BJS, BJS_pivot, numA_per_rank, numA_per_rank, &
                                   s, S_poisson

  implicit none

  ! Internal variables.
  integer                                 :: ii, jj, kk
  integer                                 :: info
  integer                                 :: block_number

  ! Flux communication across ranks.
  real(kind=dp), dimension(1:s, 1:s, nky) :: to_left, to_right, from_left, from_right

  ! Initialize Matrices
  BJS = 0.0_dp
  BJS_pivot = 0

  ! Pass information between ranks that we'll need to assemble the factors.
  to_left  = S_poisson(1:s, 1:s, 1, 1:nky)
  to_right = S_poisson(4*s - s + 1:4*s, 2*s - s + 1:2*s, numA_per_rank, 1:nky)
  call sync_ranks( to_left, to_right, s**2 * nky )

  ! Loop over the wavenumbers, computing the preconditioner for each.
  do kk = 1, nky

     ! Loop over this processes' capacitance blocks, factoring and storing the
     ! factors.
     do ii = 1, numA_per_rank

        ! Get this block's # in the global block sequence.
        block_number = numA_per_rank * rank + ii

        ! Set up the left/right interblock comms.
        if ( ii == 1 ) then
           from_left = to_left
        else
           from_left = S_poisson(3*s + 1:4*s, s+1:2*s , ii - 1, 1:nky)
        endif
        if ( ii == numA_per_rank ) then
           from_right = to_right
        else
           from_right = S_poisson(1:s, 1:s, ii + 1, 1:nky)
        endif

        ! Get the indices for building the factors.
        if ( block_number == 1 ) then

           ! Get the portion that's in this block.
           BJS(2*s + 1:4*s, 2*s + 1:3*s, ii, kk) = S_poisson(2*s + 1:4*s, s + 1:2*s, ii, kk)

           ! Apply addition from the right.
           BJS(2*s + 1:3*s, 3*s + 1:4*s, ii, kk) = from_right(:, :, kk)

           ! Add the identity.
           do jj = 2 * s + 1, 4 * s
              BJS(jj, jj, ii, kk) = BJS(jj, jj, ii, kk) + 1.0_dp
           enddo

           ! Factor this block.
           call DGETRF( 2*s, 2*s, BJS(2*s + 1, 2*s + 1, ii, kk), 4*s, &
                                  BJS_pivot(2*s + 1, ii, kk), info )

        elseif ( block_number == nsubx ) then

           ! Get the portion that's in this block.
           BJS(1:2*s, s+1:2*s, ii, kk) = S_poisson(1:2*s, 1:s, ii, kk)

           ! Apply the addition from the left.
           BJS(s+1:2*s, 1:s, ii, kk) = from_left(:, :, kk)

           ! Add the identity.
           do jj = 1, 2 * s
              BJS(jj, jj, ii, kk) = BJS(jj, jj, ii, kk) + 1.0_dp
           enddo

           ! Factor this block.
           call DGETRF(2*s, 2*s, BJS(1, 1, ii, kk), 4*s, BJS_pivot(1, ii, kk), info )

        else

           ! Get the porition that's in this block.
           BJS(1:4*s, s+1:3*s, ii, kk) = S_poisson(:, :, ii, kk)

           ! Apply the addition from the left.
           BJS(s+1:2*s, 1:s, ii, kk) = from_left(:, :, kk)

           ! Apply the addition from the right.
           BJS(2*s + 1:3*s, 3*s + 1:4*s, ii, kk) = from_right(:, :, kk)

           ! Add the identity.
           do jj = 1, 4 * s
              BJS(jj, jj, ii, kk) = BJS(jj, jj, ii, kk) + 1.0_dp
           enddo

           ! Factor this block.
           call DGETRF(4*s, 4*s, BJS(1, 1, ii, kk), 4*s, BJS_pivot(1, ii, kk), info )

         endif

     enddo
  enddo

end subroutine setup_3D_schur_preconditioner
