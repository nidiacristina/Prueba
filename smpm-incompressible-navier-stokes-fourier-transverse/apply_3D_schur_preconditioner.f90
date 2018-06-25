subroutine apply_3D_schur_preconditioner( Cx, x )
!
! Returns Cx = x + B(A\(Ex))
!
! where Ex is the extension operator, A the block diagonal part of the
! Poisson problem, and B the off-diagonal boundary block of the Poisson
! problem.
!
! This version stores the capacitance matrix in blocks, and applies the
! matvec Cx distributed across all MPI processes.

  use constants, only:         nky, nprocs, rank, nsubx
  use precision, only:         dp
  use woodbury_matrices, only: BJS, BJS_pivot, s, numA_per_rank, dimCk, &
                               aCk, zCk
  implicit none

  real(kind=dp), dimension(dimCk, 1:nky), intent(out) :: Cx
  real(kind=dp), dimension(dimCk, 1:nky), intent(in)  :: x
  real(kind=dp), dimension( 1:s, 1:nky )              :: Cx_left, Cx_right
  real(kind=dp), dimension( 1:s, 1:nky )              :: to_left, to_right, from_left, from_right
  real(kind=dp), dimension( 1:4*s )                   :: iix
  integer                                             :: ii, kk, info, first_block
  integer                                             :: aC, zC
  integer                                             :: block_number, first_block_number, last_block_number

  ! Do inter-rank comms.
  to_left  = x( 1:s, 1:nky )
  to_right = x( dimCk-s+1:dimCk, 1:nky )
  call sync_ranks( to_left, to_right, s * nky )

  ! Initialize output array to zero.
  ! This is overkill, only some of Cx are untouched; investigate this
  ! replace the below line with a targeted intialization.
  Cx = 0.0_dp

  ! Loop over wavenumbers.
  do kk = 1, nky

     ! Loop over blocks, dividing when required.
     first_block = rank * numA_per_rank + 1
     do ii = 1, numA_per_rank

        ! Get the block number.
        block_number = numA_per_rank * rank + ii

        ! Get the local extent of this block.
        aC = aCk( block_number ) - aCk( first_block ) + 1
        zC = aC + zCk( block_number ) - aCk( block_number )

        ! Set up the left/right interblock comms.
        if ( ii == 1 ) then
           from_left( :, kk ) = to_left( :, kk )
        else
           from_left( :, kk ) = x( aC - s: aC - 1, kk )
        endif
        if ( ii == numA_per_rank ) then
           from_right( :, kk ) = to_right( :, kk )
        else
           from_right( :, kk ) = x( zC + 1: zC + s, kk )
        endif

        ! If this blocks rank index is even, do a block LU solve.
        if ( mod( block_number, 2 ) == 0 ) then

           ! If there's an odd number of subdomains you'll never see the last block here.
           ! If there's an even number of subdomains that isn't divisible by four then the last block will be half the size.

           ! If we're at the last block and its supposed to be smaller, take care of this, otherwise proceed as usual.
           if ( block_number == nsubx ) then

              ! Grab the target array for this block.
              iix( s+1:2*s ) = x( aC:zC, kk )
              iix( 1:s )     = from_left( :, kk )

              ! Perform the division.
              call DGETRS( 'N', 2*s, 1, BJS( :, :, ii, kk ), 4*s, BJS_pivot( 1:2*s, ii, kk ), &
                           iix(1:2*s), 2*s, info )

              ! Store this rank's ownership in output array.
              Cx( aC:zC, kk ) = iix( s+1:2*s )

              ! Do the communication to send to the previous block, either internally or to the previous rank.
              if ( ii > 1 ) then

                 ! Intra-rank pass.
                 Cx( aC - s: aC - 1, kk ) = iix( 1:s )

              else

                 ! Inter-rank comm.
                 Cx_left( :, kk ) = iix( 1:s )

              endif

           else

              ! Grab the target subarray for this block.
              iix( s+1:3*s )   = x( aC:zC, kk )
              iix( 1:s )       = from_left( :, kk )
              iix( 3*s+1:4*s ) = from_right( :, kk )

              ! Perform the division.
              call DGETRS( 'N', 4*s, 1, BJS( :, :, ii, kk ), 4*s, BJS_pivot( 1:4*s, ii, kk ), &
                           iix(1:4*s), 4*s, info )

              ! Store this rank's ownership in output array.
              Cx( aC:zC, kk ) = iix( s+1:3*s )

              ! Do the communication to send to the previous and next block, either interally or to the previous/next rank.
              if ( ii > 1 ) then

                 ! Intra-rank pass.
                 Cx( aC - s: aC - 1, kk ) = iix( 1:s )

              else

                 ! Inter-rank comms.
                 Cx_left( :, kk ) = iix( 1:s )

              endif

              if ( ii < numA_per_rank ) then

                 ! Intra-rank pass.
                 Cx( zC + 1: zC + s, kk ) = iix( 3*s + 1: 4*s )

              else

                 ! Inter-rank comms.
                 Cx_right( :, kk ) = iix( 3*s + 1: 4*s )

              endif
           endif
        endif
     enddo

  enddo

  ! Do the inter-rank comms.
  call sync_ranks( Cx_left, Cx_right, s * nky )

  ! Put the communicated fluxes in the right place.
  !
  ! XXX: I need some logic here (maybe) about whether the first and last block_number on this rank are
  ! even or odd in the global block numbering scheme.
  !
  first_block_number = rank * numA_per_rank + 1
  last_block_number  = rank * numA_per_rank + numA_per_rank
  if ( rank > 0 .and. mod( first_block_number, 2 ) == 1 ) then
     Cx( 1:s, 1:nky ) = Cx_left
  endif
  if ( rank < nprocs - 1 .and. mod( last_block_number, 2 ) == 1 ) then
     Cx( dimCk - s + 1: dimCk, 1:nky ) = Cx_right
  endif

endsubroutine apply_3D_schur_preconditioner



