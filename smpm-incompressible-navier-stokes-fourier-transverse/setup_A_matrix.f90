subroutine setup_A_matrix
! This routine assembles the A matrix, which represents the block-diagonal
! part of the Poisson matrix.  This routine then factors each A block into its
! Schur decomposition.

  use constants, only:         rank
  use precision, only:         dp
  use woodbury_matrices, only: A_poisson, dimA, dimCk, eig_A, &
                               numA_per_rank, rpk, R_factors, QR_work_array, QR_lo, QR_hi, &
                               QR_scale, T_poisson, U_poisson, ZQR_work_array

  implicit none

  ! Internal variables.
  integer                                     :: ii
  integer                                     :: rankndx
  integer                                     :: info
  real(kind=dp), allocatable, dimension(:)    :: jjE
  real(kind=dp), allocatable, dimension(:)    :: Ej, Ej_x, Ej_z, BEj
  character(len=32)                           :: caststr
  real(kind=dp)                               :: percent

  ! Allocate the Woodbury matrices.
  allocate( R_factors(1:dimA, 1:numA_per_rank) )
  allocate( eig_A(1:dimA, 1:numA_per_rank) )
  allocate( QR_lo(1:numA_per_rank), QR_hi(1:numA_per_rank) );
  allocate( QR_scale(1:dimA, 1:numA_per_rank) )
  allocate( QR_work_array(1:dimA) )
  allocate( ZQR_work_array(1:dimA) )
  allocate( Ej(1:rpk), Ej_x(1:rpk), EJ_z(1:rpk), BEj(1:dimCk) )
  allocate( jjE(1:dimA) )

  ! Set this processes' rankndx.
  rankndx = rank + 1
  call notify( '   Computing Schur factorization of the A-matrix. ' )

  ! Factor each of the Ak blocks this rank is responsible for.
  do ii = 1, numA_per_rank

     ! Notify the user of progress.
     percent = 100.0_dp * real( ii, kind=dp ) / real( numA_per_rank, kind=dp )
     write( caststr, '(f6.1)' ) percent
     call notify('      Completed ' // trim(caststr) // ' percent of factorization.' )

     ! Factor the block diagonal component.

     ! Compute the Schur factorization of this block diagonal component.

     ! Balance A for accurate computatation of eigenvalues/vectors.
     QR_lo    = 1
     QR_hi    = dimA
     QR_scale = 0.0_dp

     T_poisson(:, :, ii) = cmplx( A_poisson(:, :, ii), kind=dp )
     call ZGEBAL( 'N', dimA, T_poisson(:, :, ii), dimA, &
                  QR_lo(ii), QR_hi(ii), QR_scale(:, ii), info )

     ! Reduce A to upper Hessenberg form.
     R_factors = 0.0_dp
     call ZGEHRD( dimA, QR_lo(ii), QR_hi(ii), T_poisson(:, :, ii), dimA, &
                  R_factors(:, ii), zQR_work_array, dimA, info )

     ! Compute the orthogonal matrix that is the product of the factorizations
     ! of the transformations that lead to the upper Hessenberg form.
     U_poisson(:, :, ii) = T_poisson(:, :, ii)
     call ZUNGHR( dimA, QR_lo(ii), QR_hi(ii), U_poisson(:, :, ii), dimA, R_factors(:, ii), ZQR_work_array, dimA, info )

     ! Compute the Schur factorization of the upper Hessenberg form of A.
     call ZHSEQR( 'S', 'V', dimA, QR_lo(ii), QR_hi(ii), &
                  T_poisson(:, :, ii), dimA, &
                  eig_A(:, ii), &
                  U_poisson(:, :, ii), dimA, ZQR_work_array, dimA, info )

     ! NOTE: At the end of this tortured series of computations, we should
     !       have the Schur decomposition of our operator A.  U_poisson
     !       contains the unitary matrix, and A_poisson contains the
     !       upper-triangular matrix in the Schur decomposition.
  enddo

end subroutine setup_A_matrix
