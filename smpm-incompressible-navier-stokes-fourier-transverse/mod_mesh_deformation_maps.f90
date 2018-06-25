module mesh_deformation_maps

  use precision, only: dp

  implicit none
  save

  ! Mesh deformation maps.
  real(kind=dp), allocatable, dimension(:)    :: x_xi, x_eta, x_xixi, x_etaeta, x_xieta
  real(kind=dp), allocatable, dimension(:)    :: z_xi, z_eta, z_xixi, z_etaeta, z_xieta
  real(kind=dp), allocatable, dimension(:)    :: detJ, detJ_xi, detJ_eta

  ! Element level inverse maps (for convenience).
  real(kind=dp), allocatable, dimension(:)    :: eta_x, xi_x, eta_z, xi_z

  ! Deformation angle vectors and their derivatives.
  real(kind=dp), allocatable, dimension(:)    :: x_xi_n,    x_eta_n,    z_xi_n,    z_eta_n
  real(kind=dp), allocatable, dimension(:)    :: x_xi_n_x,  x_xi_n_z,   x_eta_n_x, x_eta_n_z
  real(kind=dp), allocatable, dimension(:)    :: z_xi_n_x,  z_xi_n_z,   z_eta_n_x, z_eta_n_z

  ! Local distortion terms for CFL computation.
  real(kind=dp), allocatable, dimension(:)    :: delta_x, delta_z

  ! Normal and tangential vectors.
  real(kind=dp), allocatable, dimension(:, :) :: nx, tx, nz, tz
  real(kind=dp), allocatable, dimension(:, :) :: tx_x, tx_z, tz_x, tz_z

  ! Some useful metric terms to simplify the computation of the laplacian.
  real(kind=dp), allocatable, dimension(:)    :: d_eta_to_laplacian
  real(kind=dp), allocatable, dimension(:)    :: d_xi_to_laplacian
  real(kind=dp), allocatable, dimension(:)    :: d_etaeta_to_laplacian
  real(kind=dp), allocatable, dimension(:)    :: d_xixi_to_laplacian
  real(kind=dp), allocatable, dimension(:)    :: d_xieta_to_laplacian

end module mesh_deformation_maps
