module field_variables
! Contains the field variables (the unknowns) for the INS solver including the
! background current parameters

  use precision, only: dp

  implicit none
  save

  ! The field variables to be solved for.
  real(kind=dp), allocatable, dimension(:, :)    :: ux
  real(kind=dp), allocatable, dimension(:, :)    :: uy
  real(kind=dp), allocatable, dimension(:, :)    :: uz
  real(kind=dp), allocatable, dimension(:, :)    :: p
  real(kind=dp), allocatable, dimension(:, :)    :: rho
  real(kind=dp), allocatable, dimension(:, :)    :: rho_bar
  real(kind=dp), allocatable, dimension(:, :)    :: rho_bar_z
  real(kind=dp), allocatable, dimension(:, :)    :: rho_bar_zz

  ! The field variables containing the background current.
  real(kind=dp), allocatable, dimension(:,:)     :: ubc         ! Read from initial conditions/restart file.
  real(kind=dp), allocatable, dimension(:,:)     :: dubcdz      ! Read from initial conditions/restart file.
  real(kind=dp), allocatable, dimension(:,:)     :: dudx        ! Computed via gradient function.
  real(kind=dp), allocatable, dimension(:,:)     :: dudy        ! Computed via gradient function.
  real(kind=dp), allocatable, dimension(:,:)     :: dudz        ! Computed via gradient function.
  real(kind=dp), allocatable, dimension(:,:)     :: dvdx        ! Computed via gradient function.
  real(kind=dp), allocatable, dimension(:,:)     :: dvdy        ! Computed via gradient function.
  real(kind=dp), allocatable, dimension(:,:)     :: dvdz        ! Computed via gradient function.
  real(kind=dp), allocatable, dimension(:,:)     :: dwdx        ! Computed via gradient function.
  real(kind=dp), allocatable, dimension(:,:)     :: dwdy        ! Computed via gradient function.
  real(kind=dp), allocatable, dimension(:,:)     :: dwdz        ! Computed via gradient function.
  real(kind=dp), allocatable, dimension(:,:)     :: drhodx      ! Computed via gradient function.
  real(kind=dp), allocatable, dimension(:,:)     :: drhody
  real(kind=dp), allocatable, dimension(:,:)     :: drhodz

  ! The vectors to store the boundary values.
  real(kind=dp), allocatable, dimension(:, :)    :: ux_b
  real(kind=dp), allocatable, dimension(:, :)    :: uy_b
  real(kind=dp), allocatable, dimension(:, :)    :: uz_b
  real(kind=dp), allocatable, dimension(:, :)    :: rho_b

  ! The Fourier transforms of the field variables.
  complex(kind=dp), allocatable, dimension(:, :) :: Fp

  ! The intermediate variables for the KIO time-splitting.
  real(kind=dp), allocatable, dimension(:, :)    :: ux0, ux1, ux2, uy0, uy1, uy2, uz0, uz1, uz2
  real(kind=dp), allocatable, dimension(:, :)    :: Nux0, Nux1, Nux2, Nuy0, Nuy1, Nuy2, Nuz0, Nuz1, Nuz2
  real(kind=dp), allocatable, dimension(:, :)    :: Lux0, Lux1, Lux2, Luy0, Luy1, Luy2, Luz0, Luz1, Luz2
  real(kind=dp), allocatable, dimension(:, :)    :: Cux0, Cux1, Cux2, Cuy0, Cuy1, Cuy2, Cuz0, Cuz1, Cuz2
  real(kind=dp), allocatable, dimension(:, :)    :: rho0, rho1, rho2, Nrho0, Nrho1, Nrho2
  real(kind=dp), allocatable, dimension(:, :)    :: px, py, pz, rho_int, ux_int, uy_int, uz_int, div_u_int

  ! Variables for post-processing.
  real(kind=dp), allocatable, dimension(:,:)     :: divergence, enstrophy, kinetic_energy, streamfunction
  real(kind=dp), allocatable, dimension(:,:)     :: omega_x, omega_y, omega_z ! Vorticity components

end module field_variables
