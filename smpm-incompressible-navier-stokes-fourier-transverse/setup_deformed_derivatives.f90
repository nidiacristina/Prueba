subroutine setup_deformed_derivatives
! Computes metric terms for mapping derivatives of grid functions from
! physical space into the master element.

  use constants, only:             metric_tol, n, nsubz, rank
  use geom, only:                  cx, cz
  use legendre, only:              D, D2
  use mesh_deformation_maps, only: detJ, detJ_eta, detJ_xi, d_eta_to_laplacian, d_etaeta_to_laplacian, &
                                   d_xi_to_laplacian, d_xieta_to_laplacian, d_xixi_to_laplacian, &
                                   eta_x, eta_z, x_eta, x_etaeta, x_xi, x_xieta, x_xixi, xi_x, xi_z, &
                                   z_eta, z_etaeta, z_xi, z_xieta, z_xixi, &
                                   x_xi_n, x_eta_n, z_xi_n, z_eta_n, &
                                   x_xi_n_x, x_xi_n_z, x_eta_n_x, x_eta_n_z, &
                                   z_xi_n_x, z_xi_n_z, z_eta_n_x, z_eta_n_z, &
                                   nx, tx, nz, tz, tx_x, tx_z, tz_x, tz_z
  use precision, only:             dp
  use woodbury_matrices, only:     dimA, dimblock, numA_per_rank, rpk

  implicit none

  ! Block variables.
  real(kind=dp), dimension(1:dimblock) :: iix_xi, iix_eta, iiz_xi, iiz_eta
  real(kind=dp), dimension(1:dimblock) :: iix_xieta, iiz_xieta
  real(kind=dp), dimension(1:dimblock) :: iix_xixi, iiz_xixi, iix_etaeta, iiz_etaeta
  real(kind=dp), dimension(1:dimblock) :: iidetJ, iidetJ_xi, iidetJ_eta
  real(kind=dp), dimension(1:dimblock) :: iix, iiz

  ! Internal variables.
  integer                              :: iistart, iiend
  integer                              :: a, z, g

  ! Get the start and stop indices.
  iistart = 1
  iiend   = dimA * numA_per_rank

  ! Grab the data for this block.
  iix     = cx(iistart:iiend, 1)
  iiz     = cz(iistart:iiend, 1)

  ! Calculate the derivatives in xi.
  call DGEMM( 'N', 'N', n, dimblock / n, n, 1.0_dp, D,  n, iix, n, 0.0_dp, iix_xi,   n )
  call DGEMM( 'N', 'N', n, dimblock / n, n, 1.0_dp, D,  n, iiz, n, 0.0_dp, iiz_xi,   n )
  call DGEMM( 'N', 'N', n, dimblock / n, n, 1.0_dp, D2, n, iix, n, 0.0_dp, iix_xixi, n )
  call DGEMM( 'N', 'N', n, dimblock / n, n, 1.0_dp, D2, n, iiz, n, 0.0_dp, iiz_xixi, n )

  ! Permute into x/eta-first indexing.
  call perfect_shuffle( n * numA_per_rank, n * nsubz, iix )
  call perfect_shuffle( n * numA_per_rank, n * nsubz, iiz )
  call perfect_shuffle( n * numA_per_rank, n * nsubz, iix_xi )
  call perfect_shuffle( n * numA_per_rank, n * nsubz, iiz_xi )

  ! Calculate derivatives in eta.
  call DGEMM( 'N', 'N', n, dimblock / n, n, 1.0_dp, D,  n, iix,    n, 0.0_dp, iix_eta,     n )
  call DGEMM( 'N', 'N', n, dimblock / n, n, 1.0_dp, D,  n, iiz,    n, 0.0_dp, iiz_eta,     n )
  call DGEMM( 'N', 'N', n, dimblock / n, n, 1.0_dp, D2, n, iix,    n, 0.0_dp, iix_etaeta,  n )
  call DGEMM( 'N', 'N', n, dimblock / n, n, 1.0_dp, D2, n, iiz,    n, 0.0_dp, iiz_etaeta,  n )
  call DGEMM( 'N', 'N', n, dimblock / n, n, 1.0_dp, D,  n, iix_xi, n, 0.0_dp, iix_xieta,   n )
  call DGEMM( 'N', 'N', n, dimblock / n, n, 1.0_dp, D,  n, iiz_xi, n, 0.0_dp, iiz_xieta,   n )

  ! Permute back into xi/z-indexing.
  call perfect_shuffle( n * nsubz, n * numA_per_rank, iix )
  call perfect_shuffle( n * nsubz, n * numA_per_rank, iiz )
  call perfect_shuffle( n * nsubz, n * numA_per_rank, iix_xi )
  call perfect_shuffle( n * nsubz, n * numA_per_rank, iiz_xi )
  call perfect_shuffle( n * nsubz, n * numA_per_rank, iix_xieta )
  call perfect_shuffle( n * nsubz, n * numA_per_rank, iiz_xieta )
  call perfect_shuffle( n * nsubz, n * numA_per_rank, iix_eta )
  call perfect_shuffle( n * nsubz, n * numA_per_rank, iix_etaeta )
  call perfect_shuffle( n * nsubz, n * numA_per_rank, iiz_eta )
  call perfect_shuffle( n * nsubz, n * numA_per_rank, iiz_etaeta )

  ! Build the determinant of the Jacobian matrix and its derivatives.
  iidetJ     = iix_eta * iiz_xi - iiz_eta * iix_xi
  iidetJ_xi  = iix_xieta  * iiz_xi + iix_eta * iiz_xixi  - iiz_xieta  * iix_xi - iiz_eta * iix_xixi
  iidetJ_eta = iix_etaeta * iiz_xi + iix_eta * iiz_xieta - iiz_etaeta * iix_xi - iiz_eta * iix_xieta

  if ( any( abs( iidetJ ) < metric_tol ) ) then
     write(*,*) 'Warning triggered by setup_deformed_derivatives:'
     write(*,*) '   Rank ', rank, ' has a metric term with determinant less than the tolerance of ', metric_tol
  endif

  ! Build some metric terms that are needed to compute the gradient and
  ! laplacian.
  d_eta_to_laplacian(iistart:iiend) = &
       (-iiz_xi**2 * iidetJ_eta / iidetJ**3 + iiz_eta * iiz_xi * iidetJ_xi / iidetJ**3 + &
        iiz_xi * iiz_xieta / iidetJ**2 - iiz_eta * iiz_xixi / iidetJ**2 - &
        iix_xi**2 * iidetJ_eta / iidetJ**3 + iix_eta * iix_xi * iidetJ_xi / iidetJ**3 + &
        iix_xi * iix_xieta / iidetJ**2 - iix_eta * iix_xixi / iidetJ**2)

  d_xi_to_laplacian(iistart:iiend) = &
       (-iiz_eta**2 * iidetJ_xi / iidetJ**3 + iiz_eta * iiz_xi * iidetJ_eta / iidetJ**3 + &
        iiz_eta * iiz_xieta / iidetJ**2 - iiz_xi * iiz_etaeta / iidetJ**2 &
        -iix_eta**2 * iidetJ_xi / iidetJ**3 + iix_eta * iix_xi * iidetJ_eta / iidetJ**3 + &
        iix_eta * iix_xieta / iidetJ**2 - iix_xi * iix_etaeta / iidetJ**2)

  d_etaeta_to_laplacian(iistart:iiend) = (iiz_xi**2 / iidetJ**2 + iix_xi**2 / iidetJ**2)
  d_xixi_to_laplacian(iistart:iiend)   = (iiz_eta**2 / iidetJ**2 + iix_eta**2 / iidetJ**2)
  d_xieta_to_laplacian(iistart:iiend)  = (-2 * iiz_eta * iiz_xi - 2 * iix_eta * iix_xi) / iidetJ**2

  ! Store the data in the global vectors.
  x_xi(iistart:iiend)     = iix_xi
  z_xi(iistart:iiend)     = iiz_xi
  x_eta(iistart:iiend)    = iix_eta
  z_eta(iistart:iiend)    = iiz_eta
  x_xixi(iistart:iiend)   = iix_xixi
  z_xixi(iistart:iiend)   = iiz_xixi
  x_etaeta(iistart:iiend) = iix_etaeta
  z_etaeta(iistart:iiend) = iiz_etaeta
  x_xieta(iistart:iiend)  = iix_xieta
  z_xieta(iistart:iiend)  = iiz_xieta
  detJ(iistart:iiend)     = iidetJ
  detJ_xi(iistart:iiend)  = iidetJ_xi
  detJ_eta(iistart:iiend) = iidetJ_eta

  ! Compute the forward jacobians from the backwards jacobians.
  eta_x =  z_xi  / detJ
  xi_x  = -z_eta / detJ
  eta_z = -x_xi  / detJ
  xi_z  =  x_eta / detJ

  ! Compute x and z derivatives of the left/right/top/bottom boundary
  ! normal/tangential vectors.
  x_xi_n  = x_xi  / hypot( z_xi, x_xi )
  z_xi_n  = z_xi  / hypot( z_xi, x_xi )
  x_eta_n = x_eta / hypot( z_eta, x_eta ) ! tangential/normal vectors of each of the four sides
  z_eta_n = z_eta / hypot( z_eta, x_eta ) ! as needed.

  ! And their derivatives _x and _z.
  call compute_gradient( x_xi_n,   x_xi_n_x, x_xi_n_z  )
  call compute_gradient( z_xi_n,   z_xi_n_x, z_xi_n_z  )
  call compute_gradient( x_eta_n, x_eta_n_x, x_eta_n_z )
  call compute_gradient( z_eta_n, z_eta_n_x, z_eta_n_z )

  ! Compute normal and tangential vectors for use on the boundary (and _only_
  ! on the boundary).

  ! XXX: There is a problem here.  Corner points lie on two boundaries, and
  !      consequently have different normal/tangential vectors depending on
  !      which boundary you view them as being a part of.  To account for this
  !      I may need to store four copies of the normal/tangential vectors and
  !      their derivatives, or somehow otherwise come up with a way to adjust
  !      for the corners.
  !
  !      The solution here is to make (e.g.) nx(1:rpk, 1:4), nz(1:rpk, 1:4),
  !      etc. for each of the four boundaries and store them seperately.
  !      Annoying but necessary.  One way to partially alleviate the problem
  !      is to use nx1(1:n*nsubz,1:2) for the left/right boundaries and
  !      nx2(1:n*nsubx,1:2) for the top/bottom boundaries.  This will save a
  !      lot in terms of storage at the cost of adding another variable to
  !      understand and keep track of.

  ! Bottom boundary.
  a              = 1
  z              = rpk
  g              = n * nsubz
  nx(a:z:g, 1)   =  z_eta_n(a:z:g)
  nz(a:z:g, 1)   = -x_eta_n(a:z:g)
  tx(a:z:g, 1)   = -nz(a:z:g, 1)
  tz(a:z:g, 1)   =  nx(a:z:g, 1)
  tx_x(a:z:g, 1) = x_eta_n_x(a:z:g)
  tx_z(a:z:g, 1) = x_eta_n_z(a:z:g)
  tz_x(a:z:g, 1) = z_eta_n_x(a:z:g)
  tz_z(a:z:g, 1) = z_eta_n_z(a:z:g)

  ! Right boundary.
  a              = rpk - (n * nsubz - 1)
  z              = rpk - 0
  g              = 1
  nx(a:z:g, 2)   =  z_xi_n(a:z:g)
  nz(a:z:g, 2)   = -x_xi_n(a:z:g)
  tx(a:z:g, 2)   = -nz(a:z:g, 2)
  tz(a:z:g, 2)   =  nx(a:z:g, 2)
  tx_x(a:z:g, 2) = x_xi_n_x(a:z:g)
  tx_z(a:z:g, 2) = x_xi_n_z(a:z:g)
  tz_x(a:z:g, 2) = z_xi_n_x(a:z:g)
  tz_z(a:z:g, 2) = z_xi_n_z(a:z:g)

  ! Top boundary.
  a              = n * nsubz
  z              = rpk
  g              = n * nsubz
  nx(a:z:g, 3)   = -z_eta_n(a:z:g)
  nz(a:z:g, 3)   =  x_eta_n(a:z:g)
  tx(a:z:g, 3)   = -nz(a:z:g, 3)
  tz(a:z:g, 3)   =  nx(a:z:g, 3)
  tx_x(a:z:g, 3) = -x_eta_n_x(a:z:g)
  tx_z(a:z:g, 3) = -x_eta_n_z(a:z:g)
  tz_x(a:z:g, 3) = -z_eta_n_x(a:z:g)
  tz_z(a:z:g, 3) = -z_eta_n_z(a:z:g)

  ! Left boundary.
  a              = 1
  z              = n * nsubz
  g              = 1
  nx(a:z:g, 4)   = -z_xi_n(a:z:g)
  nz(a:z:g, 4)   =  x_xi_n(a:z:g)
  tx(a:z:g, 4)   = -nz(a:z:g, 4)
  tz(a:z:g, 4)   =  nx(a:z:g, 4)
  tx_x(a:z:g, 4) = -x_xi_n_x(a:z:g)
  tx_z(a:z:g, 4) = -x_xi_n_z(a:z:g)
  tz_x(a:z:g, 4) = -z_xi_n_x(a:z:g)
  tz_z(a:z:g, 4) = -z_xi_n_z(a:z:g)

end subroutine setup_deformed_derivatives
