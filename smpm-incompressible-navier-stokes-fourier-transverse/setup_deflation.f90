subroutine setup_deflation
! This implements a block-Jacobi preconditioner for the capacitance problem.
! This is the setup routine for the preconditioner.

  use constants, only:         nsubx
  use precision, only:         dp
  use woodbury_matrices, only: C_coarse, C_coarse_pivot, dimCk, uC_coarse

  implicit none

  ! Internal variables.
  integer                                     :: ii
  integer                                     :: info
  real(kind=dp), dimension(1:nsubx - 1)       :: ei
  real(kind=dp)                               :: sigmaC, tolC
  integer                                     :: maxitC
  real(kind=dp), dimension(1:nsubx - 1)       :: coarsei
  real(kind=dp), dimension(1:dimck)           :: Zei
  real(kind=dp), dimension(1:dimck)           :: SZei
  real(kind=dp), allocatable, dimension(:, :) :: L_coarse
  real(kind=dp), allocatable, dimension(:)    :: en
  real(kind=dp), allocatable, dimension(:, :) :: coarse_matrix

  call notify( 'Setting up deflation matrices.' )

  ! Build the coarse capacitance matrix column-by-column.
  do ii = 1, nsubx - 1

     ! Get the basis vector.
     ei     = 0.0_dp
     ei(ii) = 1.0_dp

     ! Apply the sequence of operators to get the matrix.
     call apply_Z( ei, Zei, 1 )
     call apply_poisson_capacitance_sparse( SZei, Zei )
     call apply_Z_transpose( SZei, coarsei, 1 )

     ! Store this column.
     C_coarse( :, ii ) = coarsei

  enddo
  allocate( coarse_matrix(1:nsubx - 1, 1:nsubx - 1) )
  coarse_matrix = C_coarse

  ! Set some constants.
  maxitC = 10
  sigmaC = 1.0e-7_dp
  tolC   = 1.0e-10_dp ! Just assume a small number as a tolerance for now.

  ! Compute the left and right kernel vectors directly.

  ! Factor the coarse matrix.
  call DGETRF( nsubx - 1, nsubx - 1, C_coarse, nsubx - 1, C_coarse_pivot, info )

  ! Build the lower triangular part of C_coarse.
  allocate( L_coarse(1:nsubx - 1, 1:nsubx - 1) )
  do ii = 1, nsubx - 1
     L_coarse(ii, ii)             = 1.0_dp
     L_coarse(ii+1:nsubx - 1, ii) = C_coarse(ii+1:nsubx - 1, ii)
  enddo

  ! Build a Euclidean basis vector en.
  allocate( en(1:nsubx-1) )
  en            = 0.0_dp
  en(nsubx - 1) = 1.0_dp

  ! Solve L^T into en.
  call DTRTRS( 'L', 'T', 'U', nsubx - 1, 1, L_coarse, nsubx - 1, en, nsubx - 1, info )

  ! Apply the permutation.
  uC_coarse = en
  call DLASWP( 1, uC_coarse, nsubx - 1, 1, nsubx - 1, C_coarse_pivot, 1 )
  uC_coarse = uC_coarse / norm2( uC_coarse )

end subroutine setup_deflation
