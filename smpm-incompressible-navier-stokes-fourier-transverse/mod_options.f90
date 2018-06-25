module options
! Contains execution options for SMPM.  The values set herein are the
! defaults.

  use precision, only: dp

  implicit none
  save

  real(kind=dp)           :: filter_order_xz = 0.0_dp                  ! Filter order in x-z plane.
  real(kind=dp)           :: filter_order_y  = 0.0_dp                  ! Filter order in y direction.
  real(kind=dp)           :: facrobin                                  ! Penalty factor for patching (Helmholtz).
  real(kind=dp)           :: facrobin_PPE                              ! Penalty factor for patching (Poisson).

  ! Boundary condition flags.
  integer, dimension(1:4) :: bc_flag_viscous_x
  integer, dimension(1:4) :: bc_flag_viscous_y
  integer, dimension(1:4) :: bc_flag_viscous_z

  ! Boundary condition flags: [bottom, right, top, left]
  !                           1 means no slip.  2 means free-slip.
  integer, dimension(1:4) :: bc_flag_diffusion
  integer, dimension(1:4) :: bc_flag_lhsgmres               = (/ 2, 2, 2, 2/)
  integer, dimension(1:4) :: bc_flag_viscous

  ! Default CFL Limits
  real(kind=dp)           :: cflzmin_lim = 0.0_dp, cflzmax_lim = 0.50_dp
  real(kind=dp)           :: cflymin_lim = 0.0_dp, cflymax_lim = 0.20_dp
  real(kind=dp)           :: cflxmin_lim = 0.0_dp, cflxmax_lim = 0.50_dp

  ! Execution options (i.e. not related to physics).
  logical                 :: adaptive_timestep              = .false.  ! Apply adaptive timestep.
  logical                 :: check_null_error               = .false.  ! Check error in null space computation.
  logical                 :: check_numerical_error          = .false.  ! Check numerical error for each piece of the solver.
  logical                 :: do_interfacial_averaging       = .false.  ! Apply interfacial averaging as a post-processing step.
  logical                 :: use_capacitance_preconditioner = .true.   ! Use block-Jacobi preconditioner for the capacitance problem.
  logical                 :: use_deflation                  = .true.   ! Use constant along interface deflation vectors in the Poisson solve.
  logical                 :: enforce_strong_velocity_bc     = .false.  ! Enforce strong velocity boundary conditions.
  integer                 :: timesteps_between_writes       = 50       ! Number of timesteps between writes of the data to disk.
  integer                 :: timesteps_between_logs         = 1        ! Number of timesteps between writing to the log file.
  logical                 :: read_bcs_from_initfile         = .false.  ! For non-homogenous boundary conditions, read the values from initfile.
  logical                 :: exact_nullspace_projection     = .false.  ! Apply the exact nullspace projection?
  logical                 :: apply_sponge_layer             = .false.  ! Apply sponge layer.

  ! Restart Options
  logical                 :: apply_restart                  = .false.  ! Specifies whether current simulation is a restart.
  integer                 :: timesteps_between_restarts     = 10       ! Number of timesteps between restart file output.

  ! Options related to reading and writing from a setup file.
  logical                 :: read_from_setupfile            = .false.
  logical                 :: write_to_setupfile             = .false.
  logical                 :: setup_and_stop                 = .false.

  ! Iterative solver options.
  real(kind=dp)           :: gmres_tol_viscous              = 1e-6_dp
  integer                 :: gmres_maxit_viscous            = 100
  real(kind=dp)           :: gmres_tol_poisson              = 1e-6_dp
  integer                 :: gmres_maxit_poisson            = 100
  integer                 :: gmres_restart_poisson          = -1
  integer                 :: gmres_restart_viscous          = -1
  integer                 :: gmres_maxit_kernel             = -1

  ! Iterative Solver options to be saved
  integer                 :: gmres_viscous_iterations_real
  integer                 :: gmres_viscous_iterations_imag
  integer                 :: gmres_diffusion_iterations_real
  integer                 :: gmres_diffusion_iterations_imag
  real(kind=dp)           :: linf_viscous_error_real
  real(kind=dp)           :: linf_viscous_error_imag
  real(kind=dp)           :: linf_diffusion_error_real
  real(kind=dp)           :: linf_diffusion_error_imag

  ! INSE Treatment options
  logical                 :: solve_momentum_equation        = .true.
  logical                 :: use_scalar_transport           = .true.
  logical                 :: use_scalar_diffusion           = .true.
  logical                 :: use_gravity_force              = .true.

  ! Files.
  character(len=100)      :: fname_init    = ''
  character(len=100)      :: fname_restart = ''
  character(len=100)      :: fname_runname = ''
  character(len=100)      :: fname_setup   = ''

end module options
