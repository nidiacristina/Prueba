module io_restart
! Contains global state for writing the restart file.

  use HDF5, only:                 hid_t

  implicit none
  save

  ! Identifiers that are used throughout the lifetime of the solver when
  ! writing to the post file.
  integer(hid_t)               :: restart_file_id

  integer(hid_t)               :: restart_grid_dataspace_id
  integer(hid_t)               :: restart_grid_2D_dataspace_id
  integer(hid_t)               :: restart_scalar_dataspace_id

  integer(hid_t)               :: restart_number_steps_dataset_id


  ! Grid parameters.
  character(len=5),  parameter :: restart_name_grid                        = "/grid"
  character(len=7),  parameter :: restart_name_n                           = "/grid/n"
  character(len=8),  parameter :: restart_name_mx                          = "/grid/mx"
  character(len=8),  parameter :: restart_name_my                          = "/grid/my"
  character(len=8),  parameter :: restart_name_mz                          = "/grid/mz"
  character(len=7),  parameter :: restart_name_x                           = "/grid/x"
  character(len=7),  parameter :: restart_name_y                           = "/grid/y"
  character(len=7),  parameter :: restart_name_z                           = "/grid/z"

  ! Group containing all of the restart variables.
  character(len=8),  parameter :: restart_name_field                       = "/restart"

  ! Restart variables for the streamwise velocity and advective terms.
  character(len=12),  parameter :: restart_name_ux0                        = "/restart/ux0"
  character(len=12),  parameter :: restart_name_ux1                        = "/restart/ux1"
  character(len=12),  parameter :: restart_name_ux2                        = "/restart/ux2"

  ! Restart variables for the transverse velocity and advective terms.
  character(len=12),  parameter :: restart_name_uy0                        = "/restart/uy0"
  character(len=12),  parameter :: restart_name_uy1                        = "/restart/uy1"
  character(len=12),  parameter :: restart_name_uy2                        = "/restart/uy2"

  ! Restart variables for the vertical velocity and advective terms.
  character(len=12),  parameter :: restart_name_uz0                        = "/restart/uz0"
  character(len=12),  parameter :: restart_name_uz1                        = "/restart/uz1"
  character(len=12),  parameter :: restart_name_uz2                        = "/restart/uz2"

  ! Restart variables for the density perturbation and advective terms.
  character(len=13),  parameter :: restart_name_rho0                       = "/restart/rho0"
  character(len=13),  parameter :: restart_name_rho1                       = "/restart/rho1"
  character(len=13),  parameter :: restart_name_rho2                       = "/restart/rho2"

  ! Restart variables for the time.
  character(len=11),  parameter :: restart_name_dt                         = "/restart/dt"
  character(len=12),  parameter :: restart_name_dt1                        = "/restart/dt1"
  character(len=12),  parameter :: restart_name_dt2                        = "/restart/dt2"
  character(len=13),  parameter :: restart_name_time                       = "/restart/time"
  character(len=18),  parameter :: restart_name_time_step                  = "/restart/time_step"

  ! Restart fields that are constant in space
  character(len=16), parameter :: restart_name_constant_fields             = "/constant_fields"

  ! Restart variable for the background density and its vertical derivative
  character(len=24), parameter :: restart_name_rho_bar                     = "/constant_fields/rho_bar"
  character(len=26), parameter :: restart_name_rho_bar_z                   = "/constant_fields/rho_bar_z"

  ! Restart variables for the background current.
  character(len=20), parameter :: restart_name_ubc                         = "/constant_fields/ubc"
  character(len=23), parameter :: restart_name_dubcdz                      = "/constant_fields/dubcdz"

  ! Restart variables for the boundary conditions
  character(len=21), parameter :: restart_name_ux_b                        = "/constant_fields/ux_b"
  character(len=21), parameter :: restart_name_uy_b                        = "/constant_fields/uy_b"
  character(len=21), parameter :: restart_name_uz_b                        = "/constant_fields/uz_b"

  ! The address for the fields that are read when restarting. Because
  ! the user determines which checkpoint to restart from, and the present
  ! restart configuration reads only a single checkpoint, we must generate
  ! a new restart file with the data from checkpoint of interest.
  ! This code is not programmed to read in the full restart file with all checkpoints.
  !
  character(len=8), parameter :: restart_ic_name_rho0                      = "/ic/rho0"
  character(len=7), parameter :: restart_ic_name_ux0                       = "/ic/ux0"
  character(len=7), parameter :: restart_ic_name_uy0                       = "/ic/uy0"
  character(len=7), parameter :: restart_ic_name_uz0                       = "/ic/uz0"

  character(len=8), parameter :: restart_ic_name_rho1                      = "/ic/rho1"
  character(len=7), parameter :: restart_ic_name_ux1                       = "/ic/ux1"
  character(len=7), parameter :: restart_ic_name_uy1                       = "/ic/uy1"
  character(len=7), parameter :: restart_ic_name_uz1                       = "/ic/uz1"

  character(len=8), parameter :: restart_ic_name_rho2                      = "/ic/rho2"
  character(len=7), parameter :: restart_ic_name_ux2                       = "/ic/ux2"
  character(len=7), parameter :: restart_ic_name_uy2                       = "/ic/uy2"
  character(len=7), parameter :: restart_ic_name_uz2                       = "/ic/uz2"

  character(len=6), parameter :: restart_ic_name_dt                        = "/ic/dt"
  character(len=7), parameter :: restart_ic_name_dt1                       = "/ic/dt1"
  character(len=7), parameter :: restart_ic_name_dt2                       = "/ic/dt2"
  character(len=8), parameter :: restart_ic_name_time                      = "/ic/time"

end module io_restart
