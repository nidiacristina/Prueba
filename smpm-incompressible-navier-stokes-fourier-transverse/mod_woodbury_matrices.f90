module woodbury_matrices
! Contains the matrices for the Woodbury solver for the Poisson equation.

  use precision, only: dp

  implicit none
  save

  ! Sub process-level.

  ! The block diagonal component.
  real(kind=dp), allocatable, dimension(:, :, :)    :: A_poisson            ! The block-diagonal component.
  complex(kind=dp), allocatable, dimension(:, :, :) :: U_poisson            ! The unitary matrix in the Schur factorization of A.
  complex(kind=dp), allocatable, dimension(:, :, :) :: T_poisson            ! The upper triangular part of the Schur factorization A.
  integer, allocatable, dimension(:, :)             :: A_poisson_pivot      ! The pivot permutations for the LU-factorizations.
  complex(kind=dp), allocatable, dimension(:, :)    :: R_factors            ! Coefficients multipliying the Householder reflector vectors.
  real(kind=dp), allocatable, dimension(:)          :: QR_work_array        ! Work array for a variety of QR routines.
  complex(kind=dp), allocatable, dimension(:)       :: ZQR_work_array       ! Work array for a variety of LAPACK routines.
  integer, allocatable, dimension(:)                :: QR_lo, QR_hi         ! Paramters for balancing the matrix (c.f. DGEBAL).
  complex(kind=dp), allocatable, dimension(:, :)    :: QR_scale             ! Scaling factor for the QR factorization.
  complex(kind=dp), allocatable, dimension(:, :)    :: eig_A                ! Complex eigenvalues of A.

  ! Set some constants.
  integer                                           :: sub_nsubx_block      ! Number of subdomain columns per sub-block.
  integer                                           :: sub_s                ! Number of rows in the off-diagonal sub-block.
  integer                                           :: sub_k                ! Dimension of the sub-block capacitance matrix.
  integer                                           :: sub_nsg              ! Number of grid-points in a sub-block.

  ! Process-level.

  ! Set some constants.
  integer                                           :: dimA
  integer                                           :: dimC
  integer                                           :: dimCk                ! The dimension of portion of the capacitance grid the local block owns.
  integer                                           :: dimLk                ! The dimension of the portion of the total grid the local block owns.
  integer                                           :: numA_per_rank
  integer                                           :: rpk                  ! Number of grid-points per rank.
  integer                                           :: arank, zrank         ! Ownership range of each rank's portion of the Poisson solve.
  integer                                           :: arankC, zrankC       ! Ownership range of each rank's portion of the capacitance solve.
  integer, allocatable, dimension(:)                :: aCk, zCk             ! Ownership range of each block in the full capacitance grid.

  ! Arrays for unbalanced MPI communication where each rank holds a different
  ! sized piece of the grid.
  integer, allocatable, dimension(:)                :: displacements, recv_counts

  ! The three-dimensional operator matrices.
  real(kind=dp), allocatable, dimension(:, :, :, :) :: S_poisson            ! (1:4s, 1:2s, 1:numA_per_rank, 1:nsuby).
  real(kind=dp), allocatable, dimension(:, :, :, :) :: BJS                  !
  integer, allocatable, dimension(:, :, :)          :: BJS_pivot            ! The block-jacobi preconditioner for the Schur matrix.
  real(kind=dp), allocatable, dimension(:, :, :)    :: S_coarse             ! The coarse matrix version of the Schur matrix.
  integer, allocatable, dimension(:, :)             :: S_coarse_pivot       ! The pivots for the coarse matrix.

  ! The sparse matrix version of the capacitance matrix.
  real(kind=dp), allocatable, dimension(:, :)       :: Eblock, Cblock
  real(kind=dp), allocatable, dimension(:, :, :)    :: C_poisson_block

  ! The preconditioner for the sparse-matrix version of the Poisson capacitance matrix.
  real(kind=dp), allocatable, dimension(:, :, :)    :: preC
  real(kind=dp), allocatable, dimension(:, :)       :: preC_pivot

  ! The left kernel vectors.
  real(kind=dp), allocatable, dimension(:)          :: uL                   ! The left kernel vector of the poisson operator.
  real(kind=dp), allocatable, dimension(:)          :: uC                   ! The left kernel vector of the capacitance matrix for the poisson operator.
  real(kind=dp), allocatable, dimension(:)          :: uC_right             ! The right kernel vector of the capacitance matrix.

  ! The solution to the capacitance matrix.
  real(kind=dp), allocatable, dimension(:,:)        :: xC_previous_real     ! The last real solution to the Schur matrix problem.
  real(kind=dp), allocatable, dimension(:,:)        :: xC_previous_imag     ! The last imag solution to the Schur matrix problem.

  ! Set some constants.
  integer                                           :: numblocks            ! Number of CPU processes allocated at run-time.
  integer                                           :: nsubx_block          ! Number of subdomain columns per process.
  integer                                           :: dimblock             ! Number of grid points per process.
  integer                                           :: s                    ! Number of rows in the off-diagonal block.
  integer                                           :: k                    ! Dimension of the capacitance matrix.

  ! Deflation preconditioner stuff.
  real(kind=dp), allocatable, dimension(:, :)       :: C_coarse             ! Coarse version of the capacitance matrix.
  real(kind=dp), allocatable, dimension(:)          :: C_coarse_pivot       ! Pivot for the factored version of C.
  real(kind=dp), allocatable, dimension(:, :)       :: C_coarse_unfactored
  real(kind=dp), allocatable, dimension(:)          :: uC_coarse            ! Null vector of the coarse problem.

end module woodbury_matrices
