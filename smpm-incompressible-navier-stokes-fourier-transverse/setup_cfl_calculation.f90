subroutine setup_cfl_calculation
! Pre-computes some local deformation terms for helping the computation of the
! CFL number for establishing stability and setting an adaptive time-step.

  use constants, only:             n, nsubz
  use legendre, only:              points
  use mesh_deformation_maps, only: delta_x, delta_z, x_eta, x_xi, z_eta, z_xi
  use precision, only:             dp
  use woodbury_matrices, only:     dimA, dimblock, numA_per_rank

  implicit none

  ! Master element grid spacing terms.
  real(kind=dp), allocatable, dimension(:) :: delta_eta, delta_xi
  real(kind=dp), allocatable, dimension(:) :: dpoints

  ! Internal variables.
  integer                                  :: ii, iistart, iiend

  ! Allocate the CFL local deformation terms.
  allocate( delta_xi(1:n * n * nsubz), delta_eta(1:n * n * nsubz) )
  allocate( delta_x(1:dimblock), delta_z(1:dimblock) )
  allocate( dpoints(1:n) )

  ! Compute the first-order central difference of the quadrature knots.
  dpoints(1) = points(2) - points(1)
  do ii = 2, n-1
     dpoints(ii) = points(ii+1) - points(ii-1)
  enddo
  dpoints(n) = points(n) - points(n-1)

  ! Set up the metric delta terms.
  do ii = 1, n * nsubz

     iistart                 = (ii - 1) * n + 1
     iiend                   = (ii - 1) * n + n
     delta_xi(iistart:iiend) = dpoints

  enddo
  delta_eta = delta_xi
  call perfect_shuffle( n * nsubz, n, delta_eta )

  ! Loop over vertical strips, computing the local deformations.
  do ii = 1, numA_per_rank

     ! Get the range.
     iistart = (ii - 1) * dimA + 1
     iiend   = (ii - 1) * dimA + dimA

     ! Compute this block's worth of metric terms.
     delta_x(iistart:iiend) = x_eta(iistart:iiend) * delta_eta + x_xi(iistart:iiend) * delta_xi
     delta_z(iistart:iiend) = z_eta(iistart:iiend) * delta_eta + z_xi(iistart:iiend) * delta_xi

  enddo

end subroutine setup_cfl_calculation
