subroutine assemble_laplacian
! Assembles the block diagonal A and off-diagonal EB components of the SMW
! decomposition of the Poisson matrix.

  use constants, only:             n, nsubz
  use legendre, only:              D, D2
  use mesh_deformation_maps, only: detJ, d_eta_to_laplacian, d_etaeta_to_laplacian, &
                                   d_xi_to_laplacian, d_xieta_to_laplacian, &
                                   d_xixi_to_laplacian, x_eta, x_xi, z_eta, z_xi
  use options, only:               bc_flag_lhsgmres, facrobin_PPE
  use precision, only:             dp
  use woodbury_matrices, only:     A_poisson, dimA, dimblock, numA_per_rank, rpk

  implicit none

  ! Internal variables.
  integer                                     :: ii, jj, iistart, iiend, ndx, jjstart, jjend
  real(kind=dp), allocatable, dimension(:)    :: ek, ek_x, ek_z, Aek
  real(kind=dp), allocatable, dimension(:)    :: iiek, iiek_x, iiek_z, iiAek
  real(kind=dp), dimension(1:dimA)            :: iiD_eta, iiD_xi, iiD_etaeta, iiD_xixi, iiD_xieta
  real(kind=dp)                               :: percent = 0.0_dp
  character(len=32)                           :: caststr

  real(kind=dp), parameter                    :: delta_poisson = 1.0_dp

  ! Grab the ownership range of this rank.
  iistart = 1
  iiend   = dimblock

  ! Allocate some variables.
  allocate( iiek(1:dimA), iiek_x(1:dimA), iiek_z(1:dimA), iiAek(1:dimA) )

  ! Loop over the internal blocks, building the Laplacian matrix for each.
  do ii = 1, n * n * nsubz

     ! Notify the user at 10% intervals.
     if ( 100.0_dp * ii / real( n * n * nsubz, kind=dp ) > percent ) then
        write( caststr, '(f6.1)' ) percent
        call notify('      Completed ' // trim(caststr) // ' percent of Laplacian assembly.' )
        percent = percent + 10.0_dp
     endif

     ! Get the row of the differentiation matrices we need.
     call getrow_IkronA( n * nsubz, D,  n, ii, iiD_xi )
     call getrow_IkronA( n * nsubz, D2, n, ii, iiD_xixi )
     call getrow_AkronI( n * nsubz, D,  n, ii, iiD_eta )
     call getrow_AkronI( n * nsubz, D2, n, ii, iiD_etaeta )
     call getrow_cross ( n * nsubz, D,  n, ii, iiD_xieta )

     do jj = 1, numA_per_rank

        ! Get the index of the metric term for this row.
        ndx = iistart + dimA * (jj - 1) + ii - 1

        ! Put them together to form the column of the deformed Laplacian.
        A_poisson(ii, :, jj) = d_eta_to_laplacian(ndx)    * iiD_eta + &
                               d_xi_to_laplacian(ndx)     * iiD_xi + &
                               d_etaeta_to_laplacian(ndx) * iiD_etaeta + &
                               d_xixi_to_laplacian(ndx)   * iiD_xixi + &
                               d_xieta_to_laplacian(ndx)  * iiD_xieta

     enddo

  enddo

  call notify('      Assembling the internal patching/boundary conditions.' )

  ! Build the patching/boundary conditions.
  !
  ! First, compute d/dx and d/dz of e_k by applying D_eta/D_xi matrices to
  ! e_k.  Second, send these derivatives of the canonical basis to
  ! apply_patching and apply_bc.  Third, split up the result into A and B
  ! matrices and send to their respective arrays.
  allocate( ek(1:rpk), ek_x(1:rpk), ek_z(1:rpk), Aek(1:rpk) )
  ndx = iistart - 1

  do jj = 1, dimA

     ! Get the necessary columns of the differentiation matrices.
     call getcol_IkronA( n * nsubz, D, n, jj, iiD_xi )
     call getcol_AkronI( n * nsubz, D, n, jj, iiD_eta )

     do ii = 1, numA_per_rank

        ! Get the ownership range of this sub-block.
        jjstart = (ii - 1) * dimA + 1
        jjend   = (ii - 1) * dimA + dimA

        ! Get the current grid point and build the canonical basis vector.
        ndx     = dimA * (ii - 1) + jj
        ek      = 0.0_dp
        ek(ndx) = 1.0_dp

        ! Construct the derivatives by grabbing the right columns of the
        ! differentiation matrices.
        ek_x                = 0.0_dp
        ek_z                = 0.0_dp
        ek_x(jjstart:jjend) = (z_xi(jjstart:jjend) * iiD_eta - z_eta(jjstart:jjend) * iiD_xi) / detJ(jjstart:jjend)
        ek_z(jjstart:jjend) = (x_eta(jjstart:jjend) * iiD_xi - x_xi(jjstart:jjend) * iiD_eta) / detJ(jjstart:jjend)

        ! Compute the patching and boundary conditions.
        Aek    = 0.0_dp
        iiAek  = 0.0_dp
        iiek   = ek(jjstart:jjend)
        iiek_x = ek_x(jjstart:jjend)
        iiek_z = ek_z(jjstart:jjend)

        call assemble_internal_patching( iiek, iiAek, delta_poisson, &
                                         facrobin_PPE, iiek_x, iiek_z, ii )
        Aek(jjstart:jjend) = iiAek
        call apply_bc( ek, Aek, bc_flag_lhsgmres, delta_poisson, ek_x, ek_z )

        ! Grab the block-diagonal component (A).
        A_poisson(:, jj, ii) = A_poisson(:, jj, ii) + Aek(jjstart:jjend)

     enddo
  enddo

end subroutine assemble_laplacian

! Some support routines for getting portions of the large sub-block
! differentiation matrices.

! Get a row of a matrix that is (I \kron A).
subroutine getrow_IkronA( dimI, A, dimA, row, u )

  use precision, only: dp

  implicit none

  integer, intent(in)                                  :: dimI, dimA, row
  real(kind=dp), dimension(1:dimA, 1:dimA), intent(in) :: A
  real(kind=dp), dimension(1:dimA*dimI), intent(out)   :: u
  integer                                              :: block, ndx, iistart, iiend

  ! Which block is this row in?
  block = floor( (row - 1) / real( dimA , kind=dp ) ) + 1

  ! Which index in this block is this row in?
  ndx = mod( row - 1, dimA ) + 1

  ! Return the row of I \kron A.
  u                = 0.0_dp
  iistart          = (block - 1) * dimA + 1
  iiend            = (block - 1) * dimA + dimA
  u(iistart:iiend) = A(ndx, :)

end subroutine getrow_IkronA

! Get a column of a matrix that is (I \kron A ).
subroutine getcol_IkronA( dimI, A, dimA, col, u )

  use precision, only: dp

  implicit none

  integer, intent(in)                                  :: dimI, dimA, col
  real(kind=dp), dimension(1:dimA, 1:dimA), intent(in) :: A
  real(kind=dp), dimension(1:dimA*dimI), intent(out)   :: u
  integer                                              :: block, ndx, iistart, iiend

  ! Which block is this row in?
  block = floor( (col - 1) / real( dimA, kind=dp ) ) + 1

  ! Which index in this block is this row in?
  ndx = mod( col - 1, dimA ) + 1

  ! Return the row of I \kron A.
  u                = 0.0_dp
  iistart          = (block - 1) * dimA + 1
  iiend            = (block - 1) * dimA + dimA
  u(iistart:iiend) = A(:, ndx)

end subroutine getcol_IkronA

! Get a row of a matrix that is (A \kron I).
subroutine getrow_AkronI( dimI, A, dimA, row, u )

  use precision, only: dp

  implicit none

  integer, intent(in)                                  :: dimI, dimA, row
  real(kind=dp), dimension(1:dimA, 1:dimA), intent(in) :: A
  real(kind=dp), dimension(1:dimA*dimI), intent(out)   :: u
  integer                                              :: block, ndx

  ! Which block is this row in?
  block = floor( (row - 1) / real( dimI, kind=dp ) ) + 1

  ! Which index in this block is this row in?
  ndx = mod( row - 1, dimI ) + 1

  ! Return the row of A \kron I.
  u                       = 0.0_dp
  u(ndx:dimI * dimA:dimI) = A(block, :)

end subroutine getrow_AkronI

! Get a column of a matrix that is (A \kron I).
subroutine getcol_AkronI( dimI, A, dimA, col, u )

  use precision, only: dp

  implicit none

  integer, intent(in)                                  :: dimI, dimA, col
  real(kind=dp), dimension(1:dimA, 1:dimA), intent(in) :: A
  real(kind=dp), dimension(1:dimA*dimI), intent(out)   :: u
  integer                                              :: block, ndx

  ! Which block is this row in?
  block = floor( (col - 1) / real( dimI , kind=dp ) ) + 1

  ! Which index in this block is this row in?
  ndx = mod( col - 1, dimI ) + 1

  ! Return the row of I \kron A.
  u                       = 0.0_dp
  u(ndx:dimI * dimA:dimI) = A(:, block)

end subroutine getcol_AkronI

! Get a row of the dense matrix (I \kron A)*(A \kron I).
subroutine getrow_cross( dimI, A, dimA, row, u )

  use precision, only: dp

  implicit none

  integer, intent(in)                                  :: dimI, dimA, row
  real(kind=dp), dimension(1:dimA, 1:dimA), intent(in) :: A
  real(kind=dp), dimension(1:dimA*dimI), intent(out)   :: u
  integer                                              :: r, ii
  real(kind=dp), dimension(1:dimA*dimI)                :: IkronA_row, AkronI_col

  real(kind=dp), external                              :: DDOT

  ! Get the row we need of IkronA.
  call getrow_IkronA( dimI, A, dimA, row, IkronA_row )

  ! Loop over the columns of AkronI to get the row we want.
  r = dimI * dimA
  do ii = 1, r
     call getcol_AkronI( dimI, A, dimA, ii, AkronI_col )
     u(ii) = DDOT( size( IkronA_row, 1 ), IkronA_row, 1, AkronI_col, 1 )
  enddo

end subroutine getrow_cross
