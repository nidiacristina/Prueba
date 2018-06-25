module profiling
! Contains elapsed times for sections of the solver.

  use precision, only: dp

  implicit none
  save

  ! Wall time spent in the setup, in aggregate and by dominant sections.
  real(kind=dp)                            :: time_setup             = 0.0_dp
  real(kind=dp)                            :: time_setup_null_basis  = 0.0_dp
  real(kind=dp)                            :: time_setup_null_error  = 0.0_dp

  ! Wall time spent in each timestep, in aggregate and by dominant sections.
  real(kind=dp)                            :: time_timesteps         = 0.0_dp
  real(kind=dp)                            :: time_field_io          = 0.0_dp
  real(kind=dp)                            :: time_restart_io        = 0.0_dp
  real(kind=dp)                            :: time_compute_diffusion = 0.0_dp
  real(kind=dp)                            :: time_compute_poisson   = 0.0_dp
  real(kind=dp)                            :: time_compute_viscous   = 0.0_dp

  ! CPU time spent during the solvers execution (both setup and timesteps) for
  ! each MPI rank.
  real(kind=dp), allocatable, dimension(:) :: time_cpu

end module profiling
