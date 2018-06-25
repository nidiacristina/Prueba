subroutine setup_3D_deflation!
! This implements a block-Jacobi preconditioner
! for the capacitance problem.  This is the setup
! routine for the preconditioner.

  use constants, only:             nky, nsubx 
  use precision, only:             dp
  use woodbury_matrices, only:     dimCk, S_coarse

  implicit none

  ! Internal variables.
  integer                                            :: ii
  real(kind=dp), dimension( 1:nsubx - 1, 1:nky )     :: ei, coarsei
  real(kind=dp), dimension( 1:dimCk,     1:nky )     :: Zei, SZei

  call notify( 'Setting up the 3D deflation matrices.' )

  ! Build the coarse schur matrix column-by-column for each wavenumber.
  do ii = 1, nsubx - 1

     ! Get the basis vector.
     ei       = 0.0_dp
     ei(ii,:) = 1.0_dp

     ! Apply the sequence of operators to get the matrix.
     call apply_Z( ei, Zei, nky )
     call apply_3D_poisson_schur_sparse( SZei, Zei )
     call apply_Z_transpose( SZei, coarsei, nky )

     ! Store this column for all the wavenumbers.
     S_coarse( :, ii, : ) = coarsei

  enddo

end subroutine setup_3D_deflation




