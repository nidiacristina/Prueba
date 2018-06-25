subroutine apply_pB_poisson_transpose( q, Bq )
! Compute the matrix-vector multiply B'q, where B is the boundary matrix in
! the SMW decomposition of the Poisson operator.
!
! (q_x,q_z) - the gradient of q.
! delta     - the viscosity.
! fac       - the multiplicative factor on the penalty coefficient.

  use constants, only:               nsubx, rank
  use precision, only:               dp
  use woodbury_matrices, only:       aCk, arankC, dimA, dimCk, numA_per_rank, rpk, s, zCk

  implicit none

  real(kind=dp), dimension(dimCk), intent(in)  :: q
  real(kind=dp), dimension(rpk), intent(inout) :: Bq
  integer                                      :: ii, ndx
  real(kind=dp), dimension(1:dimCk)            :: Ai
  real(kind=dp), dimension(1:rpk)              :: ei, ei_x, ei_z
  integer                                      :: a, z, block_number, iiblock
  real(kind=dp), dimension(1:s)                :: q_left, q_right, Ai_left, Ai_right
  real(kind=dp), dimension(1:1)                :: dot_left, dot_right

  real(kind=dp), external                      :: DDOT

  ! Loop over the vertical columns of B, computing dot products.
  Bq = 0.0_dp
  do ii = 1, dimA

     ! Build the basis vector corresponding to the ii-th column in each dimA grid.
     ei = 0.0_dp
     do iiblock = 1, numA_per_rank
        ndx     = (iiblock - 1) * dimA + ii
        ei(ndx) = 1.0_dp
     enddo

     ! Compute the gradient for all of the basis vectors at once.
     call compute_gradient( ei, ei_x, ei_z )

     ! Compute the matvec.
     Ai = 0.0_dp
     call apply_pB_poisson( ei, Ai, ei_x, ei_z, 1 )

     ! Pass some information to the left and the right.
     Ai_left   = Ai(1:s)
     Ai_right  = Ai(dimCk-s+1:dimCk)
     q_left    = q(1:s)
     q_right   = q(dimCk-s+1:dimCk)
     dot_left  = DDOT( s, q_left, 1, Ai_left, 1 )
     dot_right = DDOT( s, q_right, 1, Ai_right, 1 )
     call sync_ranks( dot_left, dot_right, 1 )

     ! Compute each of the dot products and store.
     do iiblock = 1, numA_per_rank

       ! Get the block number, and the ownership range of the indices in the capacitance grid.
       block_number = numA_per_rank * rank + iiblock
       a            = aCk(block_number) - arankC + 1
       z            = zCk(block_number) - arankC + 1
       ndx          = (iiblock - 1) * dimA + ii

       ! Compute the top half-block.
       if ( block_number > 1 ) then
          if ( iiblock > 1 ) then
             Bq(ndx) = Bq(ndx) + DDOT( s, q(a-s:a-1), 1, Ai(a-s:a-1), 1 )
          else
             Bq(ndx) = Bq(ndx) + dot_left(1)
          endif
       endif

       ! Compute the bottom half-block.
       if ( block_number < nsubx ) then
          if ( iiblock < numA_per_rank ) then
             Bq(ndx) = Bq(ndx) + DDOT( s, q(z + 1:z + s), 1, Ai(z + 1:z + s), 1 )
          else
             Bq(ndx) = Bq(ndx) + dot_right(1)
          endif
       endif

     enddo

  enddo

! XXX: the implementation below is workable and much faster but needs some
!      fixing because of comms.  I will look at this later, since this is
!      low-priority ( it is only used a few times during setup ).

!  ! NOTE: This routine basically works by getting each column of B and taking
!  !       a dot product with it.  This is really inefficient and can be
!  !       re-coded, but since it is only used in the setup of the solver, I
!  !       have not bothered to do this.
!  !
!  !       B has blocks such that Bi and Bj are disjoint if i, j are in
!  !       seperate blocks.  I use this property to get np columns of B with
!  !       one matvec by computing B*( ei + ej + ek + ... ) for a bunch of
!  !       basis vectors ek and then seperate them into the blocks.  This is
!  !       done canonically by the apply_pB_poisson matvec.
!
!  ! Loop over vertical columns of subdomains owned by this rank.
!  do ii = 1, numA_per_rank
!
!     ! Loop over this vertical strip.
!     do jj = 1, n * n * nsubz
!
!        ! Compute the index in the global grid.
!        ndx = ( ii - 1 ) * dimA + jj
!
!        ! Get the ownership range of this strip.
!        jjstart = ( ii - 1 ) * dimA + 1
!        jjend   = ( ii - 1 ) * dimA + dimA
!
!        ! Build the basis vector.
!        ei      = 0.0_dp
!        ei(ndx) = 1.0_dp
!
!        ! Build the gradient of basis vector.
!        call getcol_IkronA( n * nsubz, D, n, jj, iiD_xi  )
!        call getcol_AkronI( n * nsubz, D, n, jj, iiD_eta )
!        ei_x = 0.0_dp
!        ei_z = 0.0_dp
!        ei_x(jjstart:jjend) = ( z_xi(jjstart:jjend) * iiD_eta - z_eta(jjstart:jjend) * iiD_xi ) / detJ(jjstart:jjend)
!        ei_z(jjstart:jjend) = ( x_eta(jjstart:jjend) * iiD_xi - x_xi(jjstart:jjend) * iiD_eta ) / detJ(jjstart:jjend)
!
!        ! Grab this column of the B matrix.
!        call apply_pB_poisson( ei, Ai, ei_x, ei_z )
!
!        ! Compute the matvec and add it to the transpose matvec.
!        Bq(ndx) = pdot_product( Ai , q )
!
!     enddo
!
!  enddo

end subroutine apply_pB_poisson_transpose
