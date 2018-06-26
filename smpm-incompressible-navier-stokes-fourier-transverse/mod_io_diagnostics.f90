module io_diagnostics
! Contains global state for writing the diagnostics file.
!
! Adapted from mod_io.f90
!
! May 2018
! Gustavo Rivera

  use HDF5, only:                 hid_t

  implicit none
  save

  ! Identifiers that are used throughout the lifetime of the solver when
  ! writing to the diagnostics file.
  integer(hid_t)               :: diagnostics_file_id
  integer(hid_t)               :: diagnostics_length1_vector_dataspace_id
  integer(hid_t)               :: diagnostics_grid_2D_dataspace_id
  integer(hid_t)               :: diagnostics_scalar_dataspace_id

  integer(hid_t)               :: diagnostics_number_steps_dataset_id
  integer(hid_t)               :: diagnostics_gmres_diffusion_iterations_real_dataset_id
  integer(hid_t)               :: diagnostics_gmres_diffusion_iterations_imag_dataset_id
  integer(hid_t)               :: diagnostics_gmres_poisson_iterations_real_dataset_id
  integer(hid_t)               :: diagnostics_gmres_poisson_iterations_imag_dataset_id
  integer(hid_t)               :: diagnostics_gmres_viscous_iterations_real_dataset_id
  integer(hid_t)               :: diagnostics_gmres_viscous_iterations_imag_dataset_id
  integer(hid_t)               :: diagnostics_l2_poisson_error_dataset_id
  integer(hid_t)               :: diagnostics_l2_poisson_kernel_error_dataset_id
  integer(hid_t)               :: diagnostics_l2_poisson_schur_real_error_dataset_id
  integer(hid_t)               :: diagnostics_l2_poisson_schur_imag_error_dataset_id
  integer(hid_t)               :: diagnostics_l2_schur_kernel_error_dataset_id
  integer(hid_t)               :: diagnostics_linf_diffusion_real_error_dataset_id
  integer(hid_t)               :: diagnostics_linf_diffusion_imag_error_dataset_id
  integer(hid_t)               :: diagnostics_linf_divergence_error_dataset_id
  integer(hid_t)               :: diagnostics_linf_viscous_real_error_dataset_id
  integer(hid_t)               :: diagnostics_linf_viscous_imag_error_dataset_id
  integer(hid_t)               :: diagnostics_cfl_x_dataset_id
  integer(hid_t)               :: diagnostics_cfl_y_dataset_id
  integer(hid_t)               :: diagnostics_cfl_z_dataset_id
  integer(hid_t)               :: diagnostics_time_step_size_dataset_id
  integer(hid_t)               :: diagnostics_simulation_time_dataset_id
  integer(hid_t)               :: diagnostics_wall_time_step_dataset_id

  ! Configuration parameters.
  character(len=14), parameter :: diagnostics_name_config                        = "/configuration"
  character(len=22), parameter :: diagnostics_name_dt                            = "/configuration/delta_t"
  character(len=23), parameter :: diagnostics_name_facrobin                      = "/configuration/facrobin"
  character(len=27), parameter :: diagnostics_name_facrobin_ppe                  = "/configuration/facrobin_ppe"
  character(len=30), parameter :: diagnostics_name_filter_order_xz               = "/configuration/filter_order_xz"
  character(len=29), parameter :: diagnostics_name_filter_order_y                = "/configuration/filter_order_y"
  character(len=24), parameter :: diagnostics_name_nu                            = "/configuration/viscosity"
  character(len=29), parameter :: diagnostics_name_tend                          = "/configuration/time_threshold"
  character(len=30), parameter :: diagnostics_name_timesteps_per_write           = "/configuration/steps_per_write"
  character(len=17), parameter :: diagnostics_name_uL                            = "/configuration/uL"
  character(len=17), parameter :: diagnostics_name_uC                            = "/configuration/uC"

  ! GMRES-specific configuration parameters.
  character(len=20), parameter :: diagnostics_name_config_gmres                  = "/configuration/gmres"
  character(len=38), parameter :: diagnostics_name_poisson_tolerance             = "/configuration/gmres/poisson_tolerance"
  character(len=43), parameter :: diagnostics_name_poisson_max_iters             = "/configuration/gmres/poisson_max_iterations"
  character(len=36), parameter :: diagnostics_name_poisson_restart               = "/configuration/gmres/poisson_restart"
  character(len=38), parameter :: diagnostics_name_viscous_tolerance             = "/configuration/gmres/viscous_tolerance"
  character(len=43), parameter :: diagnostics_name_viscous_max_iters             = "/configuration/gmres/viscous_max_iterations"

  ! Execution statistics.
  character(len=10), parameter :: diagnostics_name_execution                     = "/execution"

  ! CFL Subgroup
  character(len=14), parameter :: diagnostics_name_cfl                           = "/execution/cfl"
  character(len=20), parameter :: diagnostics_name_cfl_x                         = "/execution/cfl/cfl_x"
  character(len=20), parameter :: diagnostics_name_cfl_y                         = "/execution/cfl/cfl_y"
  character(len=20), parameter :: diagnostics_name_cfl_z                         = "/execution/cfl/cfl_z"
  character(len=30), parameter :: diagnostics_name_simulation_time               = "/execution/cfl/simulation_time"
  character(len=29), parameter :: diagnostics_name_time_step_size                = "/execution/cfl/time_step_size"

  ! Error Subgroup
  character(len=16), parameter :: diagnostics_name_error                       = "/execution/error"
  character(len=48), parameter :: diagnostics_name_gmres_diffusion_real_iters  = "/execution/error/gmres_diffusion_iterations_real"
  character(len=48), parameter :: diagnostics_name_gmres_diffusion_imag_iters  = "/execution/error/gmres_diffusion_iterations_imag"
  character(len=46), parameter :: diagnostics_name_gmres_poisson_real_iters    = "/execution/error/gmres_poisson_iterations_real"
  character(len=46), parameter :: diagnostics_name_gmres_poisson_imag_iters    = "/execution/error/gmres_poisson_iterations_imag"
  character(len=46), parameter :: diagnostics_name_gmres_viscous_real_iters    = "/execution/error/gmres_viscous_iterations_real"
  character(len=46), parameter :: diagnostics_name_gmres_viscous_imag_iters    = "/execution/error/gmres_viscous_iterations_imag"
  character(len=27), parameter :: diagnostics_name_l2_poisson_error            = "/execution/error/l2_poisson"
  character(len=34), parameter :: diagnostics_name_l2_poisson_kernel_error     = "/execution/error/l2_poisson_kernel"
  character(len=38), parameter :: diagnostics_name_l2_poisson_schur_real_error = "/execution/error/l2_poisson_schur_real"
  character(len=38), parameter :: diagnostics_name_l2_poisson_schur_imag_error = "/execution/error/l2_poisson_schur_imag"
  character(len=40), parameter :: diagnostics_name_l2_poisson_schur_kernel_error = "/execution/error/l2_poisson_schur_kernel"
  character(len=36), parameter :: diagnostics_name_linf_diffusion_real_error   = "/execution/error/linf_diffusion_real"
  character(len=36), parameter :: diagnostics_name_linf_diffusion_imag_error   = "/execution/error/linf_diffusion_imag"
  character(len=32), parameter :: diagnostics_name_linf_divergence_error       = "/execution/error/linf_divergence"
  character(len=34), parameter :: diagnostics_name_linf_viscous_real_error     = "/execution/error/linf_viscous_real"
  character(len=34), parameter :: diagnostics_name_linf_viscous_imag_error     = "/execution/error/linf_viscous_imag"
  character(len=41), parameter :: diagnostics_name_total_numeric_error         = "/execution/error/total_norm_numeric_error"

  ! Wall Time subgroup
  character(len=20), parameter :: diagnostics_name_wall_time                     = "/execution/wall_time"
  character(len=29), parameter :: diagnostics_name_field_io_time                 = "/execution/wall_time/field_io"
  character(len=33), parameter :: diagnostics_name_restart_io_time               = "/execution/wall_time/restart_io"
  character(len=26), parameter :: diagnostics_name_setup_wall_time               = "/execution/wall_time/setup"
  character(len=26), parameter :: diagnostics_name_start_time                    = "/execution/wall_time/start"
  character(len=27), parameter :: diagnostics_name_step_wall_time                = "/execution/wall_time/steps"
  character(len=26), parameter :: diagnostics_name_total_wall_time               = "/execution/wall_time/total"
  character(len=36), parameter :: diagnostics_name_diffusion_solve_time          = "/execution/wall_time/solve_diffusion"
  character(len=34), parameter :: diagnostics_name_poisson_solve_time            = "/execution/wall_time/solve_poisson"
  character(len=34), parameter :: diagnostics_name_viscous_solve_time            = "/execution/wall_time/solve_viscous"
  character(len=31), parameter :: diagnostics_name_null_basis_wall_time          = "/execution/wall_time/null_basis"
  character(len=31), parameter :: diagnostics_name_null_error_wall_time          = "/execution/wall_time/null_error"

  ! Others
  character(len=23), parameter :: diagnostics_name_number_ranks                  = "/execution/number_ranks"
  character(len=25), parameter :: diagnostics_name_number_threads                = "/execution/number_threads"
  character(len=19), parameter :: diagnostics_name_cpu_time                      = "/execution/cpu_time"

end module io_diagnostics
