subroutine setup_poisson_solver
! This is the assembly step for the Poisson problem.  The A, B, C and E
! matrices are built, the A matrix is factored, and the C matrix is
! distributed.

  use constants, only:             nprocs, nsubx, rank
  use mesh_deformation_maps, only: eta_x, eta_z, x_xi_n, z_xi_n
  use precision, only:             dp
  use woodbury_matrices, only:     A_poisson, A_poisson_pivot, &
                                   C_poisson_block, dimA, dimCk, Eblock, k, &
                                   numA_per_rank, rpk, s

  implicit none

  ! Internal variables.
  integer                                  :: ii, kk
  integer                                  :: rankndx
  integer                                  :: info
  integer                                  :: rowA, rowZ
  integer                                  :: CrowA, CrowZ, a, z, aL, zL, aR, zR
  real(kind=dp), allocatable, dimension(:) :: jjE
  real(kind=dp), allocatable, dimension(:) :: Ej, Ej_x, Ej_z, BEj
  character(len=32)                        :: caststr
  real(kind=dp)                            :: percent
  integer                                  :: block_number

  ! Flux communication across ranks.
  real(kind=dp), dimension(1:s)            :: eta_x_left, eta_x_right, eta_z_left, eta_z_right
  real(kind=dp), dimension(1:s)            :: x_xi_left, x_xi_right, z_xi_left, z_xi_right
  real(kind=dp), dimension(1:2*s)          :: tau_metric_nbr, x_xi_nbr, z_xi_nbr

  ! Allocate the Woodbury matrices.
  allocate( Ej(1:rpk), Ej_x(1:rpk), EJ_z(1:rpk), BEj(1:dimCk) )
  allocate( jjE(1:dimA) )

  ! Set this processes' rankndx.
  rankndx = rank + 1
  call notify( '   Computing LU-factorizations of the A-matrices. ' )

  ! Factor each of the Ak blocks this rank is responsible for.
  do ii = 1, numA_per_rank

     ! Notify the user of progress.
     percent = 100.0_dp * real( ii, kind=dp ) / real( numA_per_rank, kind=dp )
     write( caststr, '(f6.1)' ) percent
     call notify('      Completed ' // trim(caststr) // ' percent of factorization.' )

     ! Factor the block diagonal component.
     call DGETRF( dimA, dimA, A_poisson(:, :, ii), dimA, &
                  A_poisson_pivot(:, ii), info )

  enddo

  ! Poisson operator capacitance matrix construction.
  call notify( '   Building capacitance matrix for the Poisson operator. ')

  ! Do some comms before stepping into the block assembly to get metric terms
  ! that we'll need.
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
  allocate( Eblock(dimA, 2*s) )
  do kk = 1, numA_per_rank

     ! Build the block extension operator.
     Eblock = 0.0_dp
     do ii = 1,s
        Eblock(ii, ii)                        = 1
        Eblock(dimA - ii + 1, 2 * s - ii + 1) = 1
     enddo

     ! Solve this Ak into the block.
     call DGETRS( 'n', dimA, 2 * s, A_poisson(:, :, kk), dimA, &
                  A_poisson_pivot(:, kk), Eblock, dimA, info )

     ! Get the ownership range of rows in the Ak\Ek matrix.
     a = (kk - 1) * dimA + 1
     z = (kk - 1) * dimA + dimA

     ! Get the ownership range of rows in the capacitance matrix for this
     ! capacitance block.
     block_number = rank * numA_per_rank + kk

     ! Get the ownership range of rows in the B * Ak \ Ek matrix for this
     ! capacitance block.
     if ( block_number == 1 ) then
          rowA = 1
          rowZ = 2 * s
     else if ( block_number == nsubx ) then
          rowA = k - 2 * s + 1
          rowZ = k
     else
          rowA = (block_number - 2) * 2 * s + 1
          rowZ = (block_number - 2) * 2 * s + 4 * s
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
     call assemble_capacitance_block( C_poisson_block(:, :, kk), kk, Eblock, &
                                      tau_metric_nbr, x_xi_nbr, z_xi_nbr )

   enddo
   deallocate( Eblock )

  ! Notify.
  call notify( 'Poisson solver setup complete.' )

end subroutine setup_poisson_solver
