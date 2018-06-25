subroutine apply_poisson_capacitance_sparse( Cx, x )
! Returns Cx = x + B(A\(Ex))
!
! where Ex is the extension operator, A the block diagonal part of the Poisson
! problem, and B the off-diagonal boundary block of the Poisson problem.
!
! This version stores the capacitance matrix in blocks, and applies the matvec
! Cx distributed across all MPI processes.

  use constants, only:         nprocs, nsubx, rank
  use precision, only:         dp
  use woodbury_matrices, only: C_poisson_block, dimCk, numA_per_rank, k, s

  implicit none

  real(kind=dp), dimension(dimCk), intent(out) :: Cx
  real(kind=dp), dimension(dimCk), intent(in)  :: x
  real(kind=dp), dimension(4 * s)              :: iiCx ! XXX: This is incorrect.  It
                                                          !      should be dimCk + 2 * s.
  real(kind=dp), dimension(1:s)                :: flux_L, flux_R
  integer                                      :: ii
  integer                                      :: a, z
  integer                                      :: CrowA, CrowZ, CcolA, CcolZ
  integer                                      :: Nrows, Ncols
  integer                                      :: block_number

  ! Initialize.
  iiCx = 0.0_dp
  Cx   = 0.0_dp

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
     call DGEMV( 'N', nrows, ncols, 1.0_dp, C_poisson_block(CrowA, CcolA, ii), 4*s, &
                      x(a:z), 1, 0.0_dp, iiCx, 1 )

     !! Carry out the block matvec, skipping the middle parts since they're zero.
     !if ( block_number > 1 ) then
     !   call DGEMV( 'N', s, ncols, 1.0_dp, C_poisson_block( 1:s, CcolA:CcolZ, ii ), &
     !                    s, x(a:z), 1, 0.0_dp, iiCx( 1:s ), 1 )
     !endif
     !if ( block_number < nsubx ) then
     !   if ( block_number > 1 ) then
     !      call DGEMV( 'N', s, ncols, 1.0_dp, C_poisson_block( 3 * s + 1:4 * s, CcolA:CcolZ, ii ), &
     !                       s, x(a:z), 1, 0.0_dp, iiCx( 3 * s + 1: 4 * s ), 1 )
     !   else
     !      call DGEMV( 'N', s, ncols, 1.0_dp, C_poisson_block( 3 * s + 1:4 * s, CcolA:CcolZ, ii ), &
     !                       s, x(a:z), 1, 0.0_dp, iiCx( s + 1: 2 * s ), 1 )
     !   endif
     !endif


     ! If there's something to the left, pass the top quarter block somewhere.
     if ( block_number .ne. 1 ) then

        ! If this is an internal pass, then do it right now.  If not, set up
        ! an inter-rank flux.
        if ( ii > 1 ) then
           Cx(a-s:a-1) = Cx(a-s:a-1) + iiCx(1:s)
        else
           flux_L      = iiCx(1:s)
        endif

     endif

     ! If there's something to the right, pass the bottom quarter block
     ! somewhere.
     if ( block_number .ne. nsubx ) then

        ! If this is an internal pass, then do it right now. If not, set up an
        ! inter-rank flux.
        if ( ii < numA_per_rank ) then
           if ( rank > 0 .or. ii > 1 ) then
              Cx(z+1:z+s) = Cx(z+1:z+s) + iiCx(3*s+1:4*s)
           else
              Cx(z+1:z+s) = Cx(z+1:z+s) + iiCx(s+1:2*s)
           endif
        else

           ! If this is rank 0 it should send s+2:2*s
           if ( block_number > 1 ) then
              flux_R = iiCx(3*s+1:4*s)
           else
              flux_R = iiCx(s+1:2*s)
           endif

        endif
     endif

  enddo

  ! Do some comms.  These are particularly tricky since different ranks have
  ! differing amounts of overlap with others.
  call sync_ranks( flux_L, flux_R, s )

  ! Apply the comms to each rank.
  if ( rank > 0 ) then
     Cx(1:s) = Cx(1:s) + flux_L
  endif
  if ( rank < nprocs - 1 ) then
     Cx(dimCk - s + 1:dimCk) = Cx(dimCk - s + 1:dimCk) + flux_R
  endif

  ! Add the identity.
  Cx = Cx + x

end subroutine apply_poisson_capacitance_sparse
