module constants
! Contains constants related to the grid, physics, and numerics.

  use precision, only: dp

  implicit none
  save

  integer                  :: n                    ! Number of collocation points in each direction per subdomain.
  integer                  :: nsubx                ! Number of subdomain in X direction.
  integer                  :: nsuby                ! Number of transverse grid points.
  integer                  :: nky                  ! Number of transverse wavenumbers we store.
  integer                  :: nsubz                ! Number of subdomain in Z direction.
  integer                  :: nsg                  ! Global number of collocation points (n^2 * numsub).
  integer                  :: r                    ! Total number of grid points.
  integer                  :: pd                   ! Polynomial degree (n - 1).
  real(kind=dp)            :: nu                   ! Diffusivity of momentum.
  real(kind=dp)            :: nu_d                 ! Diffusivity of mass in m^2s^-1.
  real(kind=dp)            :: rho_0                ! Density of water in kg / m^3.
  real(kind=dp)            :: sigma                ! Spectral shift for use in inverse iteration for the partial SVD.

  real(kind=dp), parameter :: g = 9.81_dp          ! Gravitational acceleration in m/s^2.
  real(kind=dp), parameter :: metric_tol = 1e-6_dp ! Lowest allowable value for the norm of the determinant of the metric terms.

  ! Some parallelism constants.
  integer                  :: nprocs, rank, root

end module constants
