subroutine setup_2D_schur( this_schur, diag_shift )
!
! This assembles the 2D Schur problem for the desired diagonal shift.
! Assume that the A matrix has been built and factored.

  use constants, only:             nprocs, rank, nsubx
  use precision, only:             dp
  use woodbury_matrices, only:     k, &
                                   s, dimA, &
                                   numA_per_rank, numA_per_rank, dimA, T_poisson, &
                                   rpk, U_poisson
  use mesh_deformation_maps, only: eta_x, eta_z, x_xi_n, z_xi_n
  use mpi, only:                   MPI_Wtime

  implicit none

  ! Dummy variables.
  real(kind=dp), dimension(1:4*s, 1:2*s, numA_per_rank), intent(out) :: this_schur
  complex(kind=dp), intent(in)                                       :: diag_shift

  ! Internal variables.
  integer                                                            :: ii, jj, kk
  integer                                                            :: info
  integer                                                            :: rowA, rowZ
  integer                                                            :: CrowA, CrowZ, a, z, aL, zL, aR, zR
  integer                                                            :: block_number
  complex(kind=dp), allocatable, dimension(:, :)                     :: QTEblock, Eblock
  real(kind=dp), allocatable, dimension(:, :)                        :: Eblock_real

  ! Flux communication across ranks.
  real(kind=dp), dimension(1:s)                                      :: eta_x_left, eta_x_right, eta_z_left, eta_z_right
  real(kind=dp), dimension(1:s)                                      :: x_xi_left, x_xi_right, z_xi_left, z_xi_right
  real(kind=dp), dimension(1:2*s)                                    :: tau_metric_nbr, x_xi_nbr, z_xi_nbr

  ! Initialize some vectors.
  !
  ! NOTE: This isn't used within this routine, or its children, though sets
  !       the initial guess for solving the capacitance problem.


  ! Do some comms before stepping into the block assembly to get metric terms that we'll need.
  x_xi_left   = x_xi_n(1:s)
  x_xi_right  = x_xi_n(rpk - s + 1:rpk)
  z_xi_left   = z_xi_n(1:s)
  z_xi_right  = z_xi_n(rpk - s + 1:rpk)
  eta_x_left  = eta_x(1:s)
  eta_x_right = eta_x(rpk - s + 1:rpk)
  eta_z_left  = eta_z(1:s)
  eta_z_right = eta_z(rpk - s + 1:rpk)

  call sync_ranks( x_xi_left, x_xi_right, s )
  call sync_ranks( z_xi_left, z_xi_right, s )
  call sync_ranks( eta_x_left, eta_x_right, s )
  call sync_ranks( eta_z_left, eta_z_right, s )

  ! Build the capacitance matrix blocks for this rank.
  allocate( Eblock( dimA, 2*s ) )
  allocate( Eblock_real( dimA, 2*s ) )
  allocate( QTEblock( dimA, 2*s ) )
  do kk = 1, numA_per_rank

     ! Build the block extension operator.
     Eblock = 0.0_dp
     do ii = 1,s
        Eblock(ii, ii)                        = cmplx( 1.0_dp, 0.0_dp, kind=dp )
        Eblock(dimA - ii + 1, 2 * s - ii + 1) = cmplx( 1.0_dp, 0.0_dp, kind=dp )
     enddo

     ! Solve this Ak into the block.

     ! Multiply transpose(U) into the RHS.
     call ZGEMM( 'C', 'N', dimA, 2 * s, dimA, cmplx( 1.0_dp, 0.0_dp, kind=dp ), U_poisson( :, :, kk ), dimA, &
                 Eblock, dimA, cmplx( 0.0_dp, 0.0_dp, kind=dp ), QTEblock, dimA )

     ! Shift the triangular system.
     do jj = 1, dimA
         T_poisson( jj, jj, kk ) = T_poisson( jj, jj, kk ) - diag_shift
     enddo

     ! Solve the upper triangular system.
     call ZTRTRS( 'U', 'N', 'N', dimA, 2 * s, T_poisson( :, :, kk ), dimA, &
                  QTEblock, dimA, info )

     ! Un-shift the triangular system.
     do jj = 1, dimA
         T_poisson( jj, jj, kk ) = T_poisson( jj, jj, kk ) + diag_shift
     enddo

     ! Multiply U into this solution.
     call ZGEMM( 'N', 'N', dimA, 2 * s, dimA, cmplx( 1.0_dp, 0.0_dp, kind=dp ), U_poisson( :, :, kk ), dimA, &
                 QTEblock, dimA, cmplx( 0.0_dp, 0.0_dp, kind=dp ), Eblock, dimA )

     ! Get the ownership range of rows in the Ak\Ek matrix.
     a = ( kk - 1 ) * dimA + 1
     z = ( kk - 1 ) * dimA + dimA

     ! Get the ownership range of rows in the capacitance matrix for this capacitance block.
     block_number = rank * numA_per_rank + kk

     ! Get the ownership range of rows in the B * Ak \ Ek matrix for this capacitance block.
     if ( block_number == 1 ) then
          rowA = 1
          rowZ = 2 * s
     else if ( block_number == nsubx ) then
          rowA = k - 2 * s + 1
          rowZ = k
     else
          rowA = ( block_number - 2 ) * 2 * s + 1
          rowZ = ( block_number - 2 ) * 2 * s + 4 * s
     endif

     ! Get the row range in the capacitance block.
     if ( rank == 0 .and. kk == 1 ) then
          CrowA = 4 * s - 2 * s + 1
          CrowZ = 4 * s
     else if ( ( rank == nprocs - 1 ) .and. ( kk == numA_per_rank ) ) then
          CrowA = 1
          CrowZ = 2 * s
     else
          CrowA = 1
          CrowZ = 4 * s
     endif

     ! Get the metric terms from this block's neighbors.
     tau_metric_nbr = 0.0_dp
     x_xi_nbr       = 0.0_dp
     z_xi_nbr       = 0.0_dp
     if ( kk == 1 ) then
        x_xi_nbr(1:s)       = x_xi_left
        z_xi_nbr(1:s)       = z_xi_left
        tau_metric_nbr(1:s) = hypot( eta_x_left, eta_z_left )
     else
        aL                  = (kk - 1) * dimA - (s - 1)
        zL                  = (kk - 1) * dimA

        x_xi_nbr(1:s)       = x_xi_n(aL:zL)
        z_xi_nbr(1:s)       = z_xi_n(aL:zL)
        tau_metric_nbr(1:s) = hypot( eta_x(aL:zL), eta_z(aL:zL) )
     endif

     if ( kk == numA_per_rank ) then
        x_xi_nbr(s+1:2*s)       = x_xi_right
        z_xi_nbr(s+1:2*s)       = z_xi_right
        tau_metric_nbr(s+1:2*s) = hypot( eta_x_right, eta_z_right )
     else
        aR                      = kk * dimA + 1
        zR                      = kk * dimA + s
        x_xi_nbr(s+1:2*s)       = x_xi_n(aR:zR)
        z_xi_nbr(s+1:2*s)       = z_xi_n(aR:zR)
        tau_metric_nbr(s+1:2*s) = hypot( eta_x(aR:zR), eta_z(aR:zR) )
     endif

     ! Assemble the Poisson capacitance block.
     Eblock_real = real( Eblock, kind=dp )
     call assemble_capacitance_block( this_schur( :, :, kk ), kk, Eblock_real, tau_metric_nbr, x_xi_nbr, z_xi_nbr )

   enddo
   deallocate( Eblock, Eblock_real, QTEblock )

end subroutine setup_2D_schur
