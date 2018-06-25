subroutine apply_3D_poisson_schur_sparse( Sx, x )
!
! Returns Sx = x + B(A\(Ex))
!
! where Ex is the extension operator, A the block diagonal part of the
! Poisson problem, and B the off-diagonal boundary block of the Poisson
! problem.
!
! This version stores the capacitance matrix in blocks, and applies the
! matvec Cx distributed across all MPI processes.

  use constants, only:         nky, nprocs, rank, nsubx
  use precision, only:         dp
  use woodbury_matrices, only: k, S_poisson, s, numA_per_rank, dimCk

  implicit none

  real(kind=dp), dimension(1:dimCk, 1:nky), intent(out) :: Sx
  real(kind=dp), dimension(1:dimCk, 1:nky), intent(in)  :: x
  real(kind=dp), dimension( 4 * s )                     :: iiCx ! XXX: This is incorrect.  
                                                                ! It should be dimCk + 2 * s.
  real(kind=dp), dimension(1:s, 1:nky)                  :: flux_L, flux_R
  integer                                               :: ii, kk
  integer                                               :: a, z
  integer                                               :: CrowA, CrowZ, CcolA, CcolZ
  integer                                               :: Nrows, Ncols
  integer                                               :: block_number

  ! Initialize.
  iiCx = 0.0_dp
  Sx   = 0.0_dp

  ! Loop over all the wavenumbers.
  do kk = 1, nky

     ! Loop over all the capacitance blocks this rank is responsible for.
     do ii = 1, numA_per_rank

        ! Set some constants for the range of this block.
        block_number = rank * numA_per_rank + ii

        ! Get the range of the input vector for this block.
        if ( rank == 0 ) then
           if ( block_number == 1 ) then
              a = 1
              z = s
           elseif ( block_number == nsubx ) then
              a = k - s + 1
              z = k
           else
              a = s + (ii - 2) * 2 * s + 1
              z = s + (ii - 2) * 2 * s + 2 * s
           endif
        elseif ( rank == nprocs - 1 .and. ii == numA_per_rank ) then
           a = dimCk - s + 1
           z = dimCk
        else
           a = (ii - 1) * 2 * s + 1
           z = (ii - 1) * 2 * s + 2 * s
        endif

        ! Get the submatrix in this Cblock we're multiplying.
        if ( block_number == 1 ) then
           CrowA = 2 * s + 1
           CrowZ = 4 * s
           CcolA = s + 1
           CcolZ = 2 * s
        elseif ( block_number == nsubx ) then
           CrowA = 1
           CrowZ = 2 * s
           CcolA = 1
           CcolZ = s
        else
           CrowA = 1
           CrowZ = 4 * s
           CcolA = 1
           CcolZ = 2 * s
        endif

        ! Carry out this block matvec.
        nrows = CrowZ - CrowA + 1
        ncols = CcolZ - CcolA + 1
        iiCx  = 0.0_dp

        ! Carry out the block matvec, skipping the middle parts since they're zero.
        if ( block_number > 1 ) then
           call DGEMV( 'N', s, ncols, 1.0_dp, S_poisson( :, :, ii, kk ), &
                            4*s, x(a:z, kk), 1, 0.0_dp, iiCx( 1:s ), 1 )
        endif
        if ( block_number < nsubx ) then
           if ( block_number > 1 ) then
              call DGEMV( 'N', s, ncols, 1.0_dp, S_poisson( 3 * s + 1, CcolA, ii, kk ), &
                               4 * s, x(a:z, kk), 1, 0.0_dp, iiCx( 3 * s + 1: 4 * s ), 1 )
           else
              call DGEMV( 'N', s, ncols, 1.0_dp, S_poisson( 3 * s + 1, CcolA, ii, kk ), &
                               4 * s, x(a:z, kk), 1, 0.0_dp, iiCx( s + 1: 2 * s ), 1 )
           endif
        endif

       ! If there's something to the left, pass the top quarter block somewhere.
        if ( block_number .ne. 1 ) then

           ! If this is an internal pass, then do it right now.  If not, set up an inter-rank flux.
           if ( ii > 1 ) then
              Sx( a-s:a-1, kk ) = Sx( a-s:a-1, kk ) + iiCx( 1:s )
           else
              flux_L(:, kk) = iiCx( 1:s )
           endif

        endif

        ! If there's something to the right, pass the bottom quarter block somewhere.
        if ( block_number .ne. nsubx ) then

           ! If this is an internal pass, then do it right now. If not, set up an inter-rank flux.
           if ( ii < numA_per_rank ) then
              if ( rank > 0 .or. ii > 1 ) then
                 Sx( z+1:z+s, kk ) = Sx( z+1:z+s, kk ) + iiCx( 3*s+1:4*s )
              else
                 Sx( z+1:z+s, kk ) = Sx( z+1:z+s, kk ) + iiCx( s+1:2*s )
              endif
           else

              ! If this is rank 0 it should send s+2:2*s
              if ( block_number > 1 ) then
                 flux_R(:, kk) = iiCx( 3*s+1:4*s )
              else
                 flux_R(:, kk) = iiCx( s+1:2*s   )
              endif

           endif
        endif

     enddo

  enddo

  ! Do some comms.  These are particularly tricky since different ranks have differing amounts of overlap with others.
  call sync_ranks( flux_L, flux_R, s * nky )

  ! Loop over wavenumbers.
  do kk = 1, nky

     ! Apply the comms to each rank.
     if ( rank > 0 ) then
        Sx( 1:s, kk ) = Sx( 1:s, kk ) + flux_L( :, kk )
     endif
     if ( rank < nprocs - 1 ) then
        Sx( dimCk - s + 1: dimCk, kk ) = Sx( dimCk - s + 1: dimCk, kk ) + flux_R( :, kk )
     endif

  enddo

  ! Add the identity.
  Sx = Sx + x

end subroutine apply_3D_poisson_schur_sparse



