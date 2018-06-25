subroutine apply_preconditioned_capacitance_shift( CTx, x_in )
! Returns CTx = C'x - sigma * x
!
! where Ex is the extension operator, A the block diagonal part of the Poisson
! problem, and B the off-diagonal boundary block of the Poisson problem.

  use constants, only:         nsubx, rank, sigma
  use precision, only:         dp
  use woodbury_matrices, only: C_poisson_block, dimCk, numA_per_rank, s

  implicit none

  real(kind=dp), dimension(dimCk), intent(out) :: CTx
  real(kind=dp), dimension(dimCk), intent(in)  :: x_in
  real(kind=dp), dimension(dimCk)              :: x
  integer                                      :: ii, counter, kk
  integer                                      :: rowA, rowZ, a, z
  integer                                      :: CrowA, CrowZ
  integer                                      :: offset
  integer                                      :: block_number
  real(kind=dp), dimension(1:s)                :: x_left, x_right, iix_left, iix_right
  real(kind=dp), dimension(1:4*s)              :: iix

  ! Initialize.
  CTx = 0.0_dp

  ! Right precondition.
  call apply_capacitance_preconditioner_transpose( CTx, x_in )
  x = CTx

  ! Before stepping into the computation loop, we need to get pieces of the
  ! target vector from neighboring ranks.
  x_left  = x(1:s)
  x_right = x(dimCk - s + 1:dimCk)
  call sync_ranks( x_left, x_right, s )

  ! Loop over the capacitance blocks owned by this rank.
  do kk = 1, numA_per_rank

     ! Get some constants.
     block_number = rank * numA_per_rank + kk

     ! Get the range of columns we'll loop over, indexed within this rank.
     if ( rank > 0 ) then
        a = 1 + (kk - 1) * 2 * s
     else
        if ( kk == 1 ) then
           a = 1
        else
           a = 1 + s + (kk - 2) * 2 * s
        endif
     endif
     if ( block_number < nsubx .and. block_number > 1 ) then
        z = a + 2 * s - 1
     else
        z = a + s - 1
     endif

     ! Set some bounds for this block.  Identify the range of rows owned by
     ! this rank and the column offset we start at.
     if ( block_number == 1 ) then
        CrowA  = 2 * s + 1
        CrowZ  = 4 * s
        offset = s
     else if ( block_number == nsubx ) then
        CrowA  = 1
        CrowZ  = 2 * s
        offset = 0
     else
        CrowA  = 1
        CrowZ  = 4 * s
        offset = 0
     endif

     ! Grab the left/right pieces.

     ! Get the left piece.
     if ( block_number > 1 ) then
        if ( kk > 1 ) then
           iix_left = x(a - s:a- 1)
        else
           iix_left = x_left
        endif
     endif

     ! Get the right piece.
     if ( block_number < nsubx ) then
        if ( kk < numA_per_rank ) then
           iix_right = x(z + 1:z + s)
        else
           iix_right = x_right
        endif
     endif

     ! Assemble the target vector.
     iix = 0.0_dp
     if ( block_number == 1 ) then

        iix(1:s)     = x(1:s)
        iix(s+1:2*s) = iix_right
        rowA         = 1
        rowZ         = 2 * s

     elseif ( block_number == nsubx ) then

        iix(2*s+1:3*s) = iix_left
        iix(3*s+1:4*s) = x(dimCk - s + 1:dimCk)
        rowA           = 2 * s + 1
        rowZ           = 4 * s

     else
        iix(s + 1:3 * s) = x(a:z)

        iix(1:s)         = iix_left
        iix(3*s + 1:4*s) = iix_right

        rowA = 1
        rowZ = 4 * s

     endif

     ! Compute the dot-products.
     counter = 1
     do ii = a, z

        CTx(ii) = dot_product( iix(rowA:rowZ), C_poisson_block(CrowA:CrowZ, offset + counter, kk) )
        counter  = counter + 1

     enddo

  enddo

  ! Add the identity.
  CTx = CTx + (1.0_dp - sigma) * x

end subroutine apply_preconditioned_capacitance_shift
