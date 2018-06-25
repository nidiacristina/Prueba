subroutine setup_deformation_maps
! Configures the deformation maps used by the solver.  This computes the
! collocation points, differentiation maps, quadrature weights, and filtering
! matrix before computing the deformed derivatives and spacing parameters needed
! by CFL calculations.

  use constants, only: n, pd
  use legendre, only:  D, D2, D3, filter_xz, points, wg
  use options, only:   filter_order_xz
  use precision, only: dp

  implicit none

  ! degree of the polynomial requested.
  pd = n - 1

  ! Generate array with collocation points.
  call notify( 'Generating collocation points.' )
  call jacobl( pd, 0.0_dp, 0.0_dp, points, n )

  ! Generate differentiation matrices.
  call notify( 'Generating differentiation matrices.' )
  call derv( pd, points, D, D2, D3, pd )

  ! Generate weights for Legendre polynomials and filter matrix.
  call notify( 'Generating quadrature weights.' )
  call quad( pd, points, wg, pd )

  ! Setup the filtering matrix.
  call setup_filter_xz( n, points, filter_xz, wg, filter_order_xz )

  ! Compute the deformation maps.
  call notify( 'Computing the deformation maps.' )
  call setup_deformed_derivatives()

  ! Do some precomputations for future CFL calculations.
  call notify( 'Computing local grid spacing.' )
  call setup_cfl_calculation()

end subroutine setup_deformation_maps
