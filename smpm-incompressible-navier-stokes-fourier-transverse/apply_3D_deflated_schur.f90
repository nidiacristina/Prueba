subroutine apply_3D_deflated_schur( PSMx, x_in )
! Apply the deflated Schur problem, PS(ky)M(ky)^-1, for each wavenumber ky.

  use constants, only:         nky
  use precision, only:         dp
  use woodbury_matrices, only: dimCk

  implicit none

  real(kind=dp), dimension(dimCk, nky), intent(out) :: PSMx
  real(kind=dp), dimension(dimCk, nky), intent(in)  :: x_in
  real(kind=dp), dimension(dimCk, nky)              :: x, y

  ! We want to apply PSM^-1.
  PSMx = 0.0_dp
  call apply_3D_schur_preconditioner( x, x_in )
  call apply_3D_poisson_schur_sparse( y, x )
  call apply_3D_P( PSMx, y )

end subroutine apply_3D_deflated_schur



