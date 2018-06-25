subroutine apply_interfacial_averaging( q )
! Applies strong, adaptive interfacial averaging.  See Diamessis (2005)
! 298-322 JCP.

  use constants, only:         n, nprocs, nsuby, nsubz, rank
  use precision, only:         dp
  use woodbury_matrices, only: dimA, numA_per_rank, rpk, s

  implicit none

  real(kind=dp), dimension(rpk, 1:nsuby), intent(inout) :: q

  real(kind=dp), dimension(1:s, 1:nsuby)                :: q_left, q_right
  integer                                               :: left, right, top, bottom, ii, jj, kk

  real(kind=dp), parameter                              :: Cave = 0.0_dp

  ! Do some comms to get the left/right fluxes.
  q_left  = q(1:s, 1:nsuby)
  q_right = q(rpk-s+1:rpk, 1:nsuby)
  call sync_ranks( q_left, q_right, s * nsuby )

  do kk = 1, nsuby

     ! Do the interfaces along eta (the horizontal interfaces) first.
     do ii = 1, nsubz - 1
        do jj = 1, n * numA_per_rank

           ! Get the locations of the/bottom side of this interface.
           bottom = (jj - 1) * n * nsubz + ii * n
           top    = bottom + 1

           ! Compute the normalized average value across this interface.
           if ( abs( q(bottom, kk) + q(top, kk) ) > 0.0_dp ) then
              if ( abs( q(bottom, kk) - q(top, kk) ) / abs( q(bottom, kk) + q(top, kk) ) > Cave ) then
                 q(bottom, kk) = 0.5_dp * (q(bottom, kk) + q(top, kk))
                 q(top, kk)    = 0.5_dp * (q(bottom, kk) + q(top, kk))
              endif
           endif
        enddo
     enddo

     ! Do the internal interfaces along xi (the vertical interfaces) first.
     do ii = 1, numA_per_rank - 1
        do jj = 1, s

           ! Get the locations of the left/right sides of this interface.
           left  = ii * dimA     - jj + 1
           right = ii * dimA + s - jj + 1

           ! Compute the normalized average value across this interface.
           if ( abs( q(left, kk) + q(right, kk) ) > 0.0_dp ) then
              if ( abs( q(left, kk) - q(right, kk) ) / abs( q(left, kk) + q(right, kk) ) > Cave ) then
                 q(left, kk)  = 0.5_dp * (q(left, kk) + q(right, kk))
                 q(right, kk) = 0.5_dp * (q(left, kk) + q(right, kk))
              endif
           endif
        enddo
     enddo

     ! Do the inter-rank interfaces.
     if ( rank > 0 ) then
        do jj = 1, s
           ! Compute the normalized average values across this interface.
           if ( abs( q_left(jj, kk) + q(jj, kk) ) > 0.0_dp ) then
              if ( abs( q_left(jj, kk) - q(jj, kk) ) / abs( q_left(jj, kk) + q(jj, kk) ) > Cave ) then
                 q(jj, kk) = 0.5_dp * (q_left(jj, kk) + q(jj, kk))
              endif
           endif
        enddo
     endif

     if ( rank < nprocs - 1 ) then
        do jj = 1, s

           ! Get the index of the right interface.
           right = rpk - s + jj

           ! Compute the normalized average values across this interface.
           if ( abs( q_right(jj, kk) + q(right, kk) ) > 0.0_dp ) then
              if ( abs( q_right(jj, kk) - q(right, kk) ) / abs( q_right(jj, kk) + q(right, kk) ) > Cave ) then
                 q(right, kk) = 0.5_dp * (q_right(jj, kk) + q(right, kk))
              endif
           endif
        enddo
     endif

  enddo

end subroutine apply_interfacial_averaging
