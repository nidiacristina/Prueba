module io_field
! Contains global state for writing the field file as well as constant field
! file field names.
!
! Adapted from mod_io.f90
!
! May 2018
! Gustavo Rivera

  use HDF5, only:                 hid_t

  implicit none
  save

  ! Identifiers that are used throughout the lifetime of the solver when
  ! writing to the field file.
  integer(hid_t)               :: field_file_id
  integer(hid_t)               :: field_grid_dataspace_id
  integer(hid_t)               :: field_grid_2D_dataspace_id
  integer(hid_t)               :: field_length1_vector_dataspace_id
  integer(hid_t)               :: field_scalar_dataspace_id

  integer(hid_t)               :: field_number_steps_dataset_id
  integer(hid_t)               :: field_time_step_size_dataset_id
  integer(hid_t)               :: field_simulation_time_dataset_id
  integer(hid_t)               :: field_wall_time_step_dataset_id

  ! Initial condition file variables.
  character(len=9),  parameter :: ic_name_grid_x                           = "/grid/x"
  character(len=9),  parameter :: ic_name_grid_y                           = "/grid/y"
  character(len=9),  parameter :: ic_name_grid_z                           = "/grid/z"
  character(len=10), parameter :: ic_name_dubcdz                           = "/ic/dubcdz"
  character(len=9),  parameter :: ic_name_rho                              = "/ic/rho"
  character(len=13), parameter :: ic_name_rho_bar                          = "/ic/rho_bar"
  character(len=15), parameter :: ic_name_rho_bar_z                        = "/ic/rho_bar_z"
  character(len=7),  parameter :: ic_name_ubc                              = "/ic/ubc"
  character(len=6),  parameter :: ic_name_ux                               = "/ic/ux"
  character(len=6),  parameter :: ic_name_uy                               = "/ic/uy"
  character(len=6),  parameter :: ic_name_uz                               = "/ic/uz"

  ! Grid parameters.
  character(len=5),  parameter :: field_name_grid                          = "/grid"
  character(len=7),  parameter :: field_name_n                             = "/grid/n"
  character(len=8),  parameter :: field_name_mx                            = "/grid/mx"
  character(len=8),  parameter :: field_name_my                            = "/grid/my"
  character(len=8),  parameter :: field_name_mz                            = "/grid/mz"
  character(len=7),  parameter :: field_name_x                             = "/grid/x"
  character(len=7),  parameter :: field_name_y                             = "/grid/y"
  character(len=7),  parameter :: field_name_z                             = "/grid/z"

  ! Field variables.
  !
  ! NOTE: rho, ux, uy, and uz are relative names since the group that contains
  !       them is dynamically constructed at run-time from
  !       field_name_step_base.
  !
  ! NOTE: step #0 (/field/step0) represents the initial conditions.  We reuse
  !       the field variable writing code to create the group and write out
  !       rho, ux, uy, and uz, though we need to explicitly write out rho_bar,
  !       ubc, and dubcdz after the fact.  This means we don't have to
  !       explicitly enumerate those variables.
  character(len=6),  parameter :: field_name_field                         = "/field"
  character(len=12), parameter :: field_name_field_initial_conditions      = "/field/step0"
  character(len=19), parameter :: field_name_initial_conditions_dubcdz     = "/field/step0/dubcdz"
  character(len=20), parameter :: field_name_initial_conditions_rho_bar    = "/field/step0/rho_bar"
  character(len=22), parameter :: field_name_initial_conditions_rho_bar_z  = "/field/step0/rho_bar_z"
  character(len=16), parameter :: field_name_initial_conditions_ubc        = "/field/step0/ubc"
  character(len=19), parameter :: field_name_number_steps                  = "/field/number_steps"
  character(len=3),  parameter :: field_name_rho                           = "rho"
  character(len=2),  parameter :: field_name_ux                            = "ux"
  character(len=2),  parameter :: field_name_uy                            = "uy"
  character(len=2),  parameter :: field_name_uz                            = "uz"
  character(len=20), parameter :: field_name_step_base                     = "/field/step" ! Additional length since we append the timestep count to it.
  character(len=4),  parameter :: field_name_time                          = "time"
  character(len=9),  parameter :: field_name_time_step                     = "time_step"

end module io_field
