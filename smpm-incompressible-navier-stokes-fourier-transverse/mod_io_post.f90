module io_post
! Contains global state for writing the post file as well as constant post
! file post names.

  use HDF5, only:                 hid_t

  implicit none
  save

  ! Identifiers that are used throughout the lifetime of the solver when
  ! writing to the post file.
  integer(hid_t)               :: post_file_id

  integer(hid_t)               :: post_grid_dataspace_id
  integer(hid_t)               :: post_scalar_dataspace_id

  ! Initial condition file variables.
  character(len=9),  parameter :: ic_name_grid_x                           = "/grid/x"
  character(len=9),  parameter :: ic_name_grid_y                           = "/grid/y"
  character(len=9),  parameter :: ic_name_grid_z                           = "/grid/z"

  ! Paths to the grid's datasets.

  ! Grid parameters.
  character(len=5),  parameter :: post_name_grid                          = "/grid"
  character(len=7),  parameter :: post_name_n                             = "/grid/n"
  character(len=8),  parameter :: post_name_mx                            = "/grid/mx"
  character(len=8),  parameter :: post_name_my                            = "/grid/my"
  character(len=8),  parameter :: post_name_mz                            = "/grid/mz"
  character(len=7),  parameter :: post_name_x                             = "/grid/x"
  character(len=7),  parameter :: post_name_y                             = "/grid/y"
  character(len=7),  parameter :: post_name_z                             = "/grid/z"

  ! Field variables.
  !
  ! NOTE: rho, ux, uy, and uz are relative names since the group that contains
  !       them is dynamically constructed at run-time from
  !       post_name_step_base.
  !
  ! NOTE: step #0 (/post/step0) represents the initial conditions.  We reuse the
  !       post variable writing code to create the group and write out rho, ux, and
  !       uz, though we need to explicitly write out rho_bar after the fact.  This
  !       means we don't have to explicitly enumerate those variables.
  character(len=5),  parameter :: post_name_post                         = "/post"
  character(len=9),  parameter :: post_name_div                          = "/post/div"
  character(len=15), parameter :: post_name_enstrophy                    = "/post/enstrophy"
  character(len=20), parameter :: post_name_kinetic_energy               = "/post/kinetic_energy"
  character(len=13), parameter :: post_name_omega_x                      = "/post/omega_x"
  character(len=13), parameter :: post_name_omega_y                      = "/post/omega_y"
  character(len=13), parameter :: post_name_omega_z                      = "/post/omega_z"
  character(len=12), parameter :: post_name_stream                       = "/post/stream"

  ! Velocity Group
  character(len=8),  parameter :: post_name_ux                           = "/post/ux"
  character(len=8),  parameter :: post_name_uy                           = "/post/uy"
  character(len=8),  parameter :: post_name_uz                           = "/post/uz"

  character(len=10), parameter :: post_name_dudx                         = "/post/dudx"
  character(len=10), parameter :: post_name_dudy                         = "/post/dudy"
  character(len=10), parameter :: post_name_dudz                         = "/post/dudz"

  character(len=10), parameter :: post_name_dvdx                         = "/post/dvdx"
  character(len=10), parameter :: post_name_dvdy                         = "/post/dvdy"
  character(len=10), parameter :: post_name_dvdz                         = "/post/dvdz"

  character(len=10), parameter :: post_name_dwdx                         = "/post/dwdx"
  character(len=10), parameter :: post_name_dwdy                         = "/post/dwdy"
  character(len=10), parameter :: post_name_dwdz                         = "/post/dwdz"

  ! Pressure group
  character(len=9),  parameter :: post_name_p                            = "/post/p"

  ! Density Group
  character(len=9),  parameter :: post_name_rho                          = "/post/rho"

  character(len=13), parameter :: post_name_drhodx                       = "/post/drhodx"
  character(len=13), parameter :: post_name_drhody                       = "/post/drhody"
  character(len=13), parameter :: post_name_drhodz                       = "/post/drhodz"

  character(len=13), parameter :: post_name_rho_bar                      = "/post/rho_bar"

  ! Field time: physical time of the simulation
  character(len=16), parameter :: post_name_field_time                   = "/post/field_time"

end module io_post
