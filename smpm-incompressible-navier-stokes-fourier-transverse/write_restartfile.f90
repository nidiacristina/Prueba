! NOTE: The calls to h5fflush_f() do not guarantee correct operation in the
!       case of a single writer, multiple reader scenario - it merely reduces
!       the window where the readers will see something different than they
!       expect.
!
!       THIS IS JUST A BANDAID ATTEMPTING TO HIDE THE FACT THAT FUNDAMENTALLY
!       HDF5 CANNOT HANDLE A SINGLE WRITER WITH MULTIPLE READERS IN ALL CASES.
!
!       It is quite possible that the reader can see an invalid view of the
!       file as it is being written, which can be seen when the time per
!       timestep is small and the number of timesteps between writes is
!       tiny.  In very limited testing, the invalid view was simply that
!       the /restart/ group appeared to not contain any sub-groups, nor the
!       timestep count dataset.
!
!       The Single Writer, Multiple Reader (SWMR) interface being developed
!       by the HDF5 group (as of spring 2015) will not solve this problem.
!       According to section 3 of the most recent version of their user's
!       guide:
!
!     https://www.hdfgroup.org/HDF5/docNewFeatures/UG-HDF5-SWMR-20130629-v3.pdf
!
!       Groups and datasets cannot be added to a file accessed for SWMR.
!       Currently, the restart file is written such that each timestep written
!       lives within its own group.

subroutine open_restart_file( restart_file_name, success_flag )
! Opens the restart file for writing throughout the solver's life time.  If
! the specified restart file already exists, it will be overwritten.  Nothing
! is written to the restart file.
!
! This routine opens both the scalar and grid data spaces.

  use constants, only:  n, nsubx, nsuby, nsubz, rank, root
  use HDF5, only:       H5F_ACC_TRUNC_F, H5S_SCALAR_F, &
                        hid_t, hsize_t, &
                        h5fcreate_f, &
                        h5screate_f, h5screate_simple_f
  use io_restart, only: restart_file_id, &
                        restart_grid_dataspace_id, restart_grid_2D_dataspace_id, restart_scalar_dataspace_id

  implicit none

  character(len=*), intent(in)        :: restart_file_name
  logical, intent(out)                :: success_flag

  integer                             :: data_rank
  integer(hsize_t), dimension(3)      :: data_dimensions
  integer                             :: hdf5_status

  ! We assume that since we were able to open the file, everything is fine.
  !
  ! NOTE: Checking each HDF5 operation would *greatly* bloat the code and
  !       isn't likely worth the effort since a non-programmatic failure in
  !       the middle of a series of HDF5 operations is almost certainly
  !       indicative of an environmental problem that is beyond our
  !       control.  Spiraling out of control in a glorious ball of flames
  !       is perfectly acceptable in that situation.


  if (rank .eq. root) then

     ! Create our output file, overwriting any existing one.
     call H5FCREATE_F( trim (restart_file_name ), &
                       H5F_ACC_TRUNC_F, restart_file_id, hdf5_status )

     ! Check if file was created succesfully. If not, stop execution.
     success_flag = (hdf5_status >= 0)
     if (success_flag .eqv. .false. ) then
        return
     end if

     ! Create a scalar dataspace that we'll use to create scalar datatypes
     ! with.  each of the non-grid variables are tied to this dataspace.
     call H5SCREATE_F( H5S_SCALAR_F, restart_scalar_dataspace_id, hdf5_status )

     ! Create a 3D dataspace that we'll use for grid-interface variables
     ! (coordinates, velocities, density, etc).
     data_rank          = 3
     data_dimensions(1) = n * nsubz
     data_dimensions(2) = n * nsubx
     data_dimensions(3) = nsuby
     call H5SCREATE_SIMPLE_F( data_rank, data_dimensions, restart_grid_dataspace_id, &
                              hdf5_status )

     ! Create a 2D dataspace that we'll use for transverse-independent grid
     ! variable (background current and derivatives).
     data_rank          = 2
     call H5SCREATE_SIMPLE_F( data_rank, data_dimensions, restart_grid_2D_dataspace_id, &
                              hdf5_status )

   endif

end subroutine open_restart_file

subroutine write_restart_grid
! Writes out the solver's grid to the restart file's header.

  use constants, only:               n, nsg, nsubx, nsuby, nsubz, rank, root
  use geom, only:                    cx, cz
  use HDF5, only:                    H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER, H5T_STD_I32LE, &
                                     hid_t, hsize_t, &
                                     h5dclose_f, h5dcreate_f, h5dwrite_f, &
                                     h5gclose_f, h5gcreate_f, &
                                     h5screate_f, h5screate_simple_f
  use io_restart, only:              restart_file_id, &
                                     restart_grid_dataspace_id, restart_scalar_dataspace_id, &
                                     restart_name_grid, restart_name_n, restart_name_mx, restart_name_my, &
                                     restart_name_mz, restart_name_x, restart_name_y, restart_name_z
  use precision, only:               dp
  use transverse, only:              cy

  implicit none

  real(kind=dp), allocatable, dimension(:, :) :: full_cx, full_cy, full_cz

  integer(hsize_t), dimension(3)              :: data_dimensions

  ! HDF5 identifiers and status code.
  integer(hid_t)                              :: dataset_id
  integer(hid_t)                              :: group_id
  integer                                     :: hdf5_status

  ! Allocate space for a consolidated grid on the root.
  if (rank == root )then
     allocate( full_cx(1:nsg, 1:nsuby), full_cy(1:nsg, 1:nsuby), full_cz(1:nsg, 1:nsuby) )
  else
     ! Avoid warnings about using uninitialized memory.  Non-root ranks wont
     ! touch the full grid during the collective operation.
     allocate( full_cx(1, 1), full_cy(1, 1), full_cz(1, 1) )
  end if

  ! Pull the entirety of the grid onto the root rank for writing.
  call gather_3D_array( cx, full_cx )
  call gather_3D_array( cy, full_cy )
  call gather_3D_array( cz, full_cz )

  if (rank == root) then
     ! Specify the dimensions of our grid variables, full_cx, full_cy,
     ! full_cz.  Note that these are ignored when specified with our
     ! scalar variables.
     data_dimensions(1) = n * nsubz
     data_dimensions(2) = n * nsubx
     data_dimensions(3) = nsuby

     ! Create the grid group that holds all of the variables we're
     ! writing out.
     call H5GCREATE_F( restart_file_id, restart_name_grid, group_id, hdf5_status )
     call H5GCLOSE_F( group_id, hdf5_status )

     ! Create a scalar for the collocation dimension, n.
     call H5DCREATE_F( restart_file_id, restart_name_n, H5T_STD_I32LE, &
                       restart_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_INTEGER, n, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar for the number of x sub-domains, nsubx.
     call H5DCREATE_F( restart_file_id, restart_name_mx, H5T_STD_I32LE, &
                       restart_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_INTEGER, nsubx, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's x positions,
     ! full_cx.
     call H5DCREATE_F( restart_file_id, restart_name_x, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, full_cx, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar for the number of z sub-domains, nsubz.
     call H5DCREATE_F( restart_file_id, restart_name_mz, H5T_STD_I32LE, &
                       restart_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_INTEGER, nsubz, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's z positions,
     ! cz4disk.
     call H5DCREATE_F( restart_file_id, restart_name_z, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, full_cz, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar for the number of y grid points, nsuby.
     call H5DCREATE_F( restart_file_id, restart_name_my, H5T_STD_I32LE, &
                       restart_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_INTEGER, nsuby, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's y positions,
     ! full_cy.
     call H5DCREATE_F( restart_file_id, restart_name_y, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, full_cy, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

  end if

  deallocate( full_cx, full_cy, full_cz )

end subroutine write_restart_grid


subroutine write_restart_constant_fields
! Writes out the solver's constants fields to the restart file's header.

  use constants, only:               n, nsg, nsubx, nsuby, nsubz, rank, root
  use field_variables, only:         dubcdz, rho_bar, rho_bar_z, ubc, ux_b, uy_b, uz_b
  use HDF5, only:                    H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER, H5T_STD_I32LE, &
                                     hid_t, hsize_t, &
                                     h5dclose_f, h5dcreate_f, h5dwrite_f, &
                                     h5gclose_f, h5gcreate_f, &
                                     h5screate_f, h5screate_simple_f
  use mpi, only:                     MPI_COMM_WORLD, MPI_DOUBLE_PRECISION
  use io_restart, only:              restart_file_id, &
                                     restart_grid_dataspace_id, &
                                     restart_name_constant_fields, &
                                     restart_name_dubcdz, &
                                     restart_name_rho_bar, restart_name_rho_bar_z, &
                                     restart_name_ubc, &
                                     restart_name_ux_b, restart_name_uy_b, restart_name_uz_b
  use precision, only:               dp

  implicit none

  real(kind=dp), allocatable, dimension(:, :) :: dubcdz4disk
  real(kind=dp), allocatable, dimension(:, :) :: rhobar4disk, rhobarz4disk
  real(kind=dp), allocatable, dimension(:, :) :: ubc4disk, uxb4disk, uyb4disk, uzb4disk
  integer(hsize_t), dimension(3)              :: data_dimensions

  ! HDF5 identifiers and status code.
  integer(hid_t)                              :: dataset_id
  integer(hid_t)                              :: group_id
  integer                                     :: hdf5_status

  ! Allocate space for a consolidated grid on the root.
  if (rank == root )then
     allocate( dubcdz4disk(1:nsg, 1:nsuby) )
     allocate( rhobar4disk(1:nsg, 1:nsuby), rhobarz4disk(1:nsg, 1:nsuby) )
     allocate( ubc4disk(1:nsg, 1:nsuby) )
     allocate( uxb4disk(1:nsg, 1:nsuby), uyb4disk(1:nsg, 1:nsuby), uzb4disk(1:nsg,1:nsuby) )
  else
     ! Avoid warnings about using uninitialized memory.  Non-root ranks wont
     ! touch the full grid during the collective operation.
     allocate( dubcdz4disk(1:1, 1:1) )
     allocate( rhobar4disk(1:1, 1:1), rhobarz4disk(1:1, 1:1) )
     allocate( ubc4disk(1:1, 1:1) )
     allocate( uxb4disk(1:1, 1:1), uyb4disk(1:1, 1:1), uzb4disk(1:1,1:1) )
  end if

  ! The background density, background current, and it's derivative are constant in the transverse
  ! dimension.
  call gather_3D_array( dubcdz, dubcdz4disk )
  call gather_3D_array( rho_bar, rhobar4disk )
  call gather_3D_array( rho_bar_z, rhobarz4disk )
  call gather_3D_array( ubc, ubc4disk )
  call gather_3D_array( ux_b, uxb4disk )
  call gather_3D_array( uy_b, uyb4disk )
  call gather_3D_array( uz_b, uzb4disk )


  if (rank == root) then
     ! Specify the dimensions of our grid variables.
     data_dimensions(1) = n * nsubz
     data_dimensions(2) = n * nsubx
     data_dimensions(3) = nsuby

     ! Create the grid group that holds all of the variables we're
     ! writing out.
     call H5GCREATE_F( restart_file_id, restart_name_constant_fields, group_id, hdf5_status )

     ! Create a scalar grid variable for the derivative of the background current,
     ! dubcdz4disk.
     call H5DCREATE_F( restart_file_id, restart_name_dubcdz, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, dubcdz4disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's background density profile
     ! rhobar4disk.
     call H5DCREATE_F( restart_file_id, restart_name_rho_bar, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, rhobar4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's vertical derivative of the background
     ! density profile, rhobarz4disk.
     call H5DCREATE_F( restart_file_id, restart_name_rho_bar_z, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, rhobarz4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the background current,
     ! ubc4disk.
     call H5DCREATE_F( restart_file_id, restart_name_ubc, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, ubc4disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the horizontal velocity boundary condition,
     ! uxb4disk.
     call H5DCREATE_F( restart_file_id, restart_name_ux_b, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, uxb4disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the transverse velocity boundary condition,
     ! uyb4disk.
     call H5DCREATE_F( restart_file_id, restart_name_uy_b, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, uyb4disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the vertical velocity boundary condition,
     ! uzb4disk.
     call H5DCREATE_F( restart_file_id, restart_name_uz_b, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, uzb4disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )


     call H5GCLOSE_F( group_id, hdf5_status )

  end if

  deallocate( dubcdz4disk )
  deallocate( rhobar4disk, rhobarz4disk )
  deallocate( ubc4disk )
  deallocate( uxb4disk, uyb4disk, uzb4disk )

end subroutine write_restart_constant_fields


subroutine write_restart_data
! Writes out the solver's data necessary for restarting run

  use constants, only:               n, nsg, nsubx, nsuby, nsubz, rank, root
  use field_variables, only:         rho0, rho1, rho2, &
                                     ux0, ux1, ux2, &
                                     uy0, uy1, uy2, &
                                     uz0, uz1, uz2
  use HDF5, only:                    H5T_IEEE_F64LE, H5F_SCOPE_LOCAL_F, H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER, H5T_STD_I32LE, &
                                     hid_t, hsize_t, &
                                     h5dclose_f, h5dcreate_f, h5dwrite_f, &
                                     h5fflush_f, &
                                     h5gclose_f, h5gcreate_f, &
                                     h5screate_f, h5screate_simple_f
  use io_restart, only:              restart_file_id, restart_scalar_dataspace_id, restart_name_field, &
                                     restart_grid_dataspace_id, &
                                     restart_name_dt, restart_name_dt1, restart_name_dt2, &
                                     restart_name_rho0, restart_name_rho1, restart_name_rho2, &
                                     restart_name_ux0, restart_name_ux1, restart_name_ux2, &
                                     restart_name_uy0, restart_name_uy1, restart_name_uy2, &
                                     restart_name_uz0, restart_name_uz1, restart_name_uz2, &
                                     restart_name_time, restart_name_time_step
  use mpi, only:                     MPI_COMM_WORLD, MPI_DOUBLE_PRECISION
  use precision, only:               dp
  use timestepping, only:            dt, dt1, dt2, simulation_time, time_ndx

  implicit none

  ! Vector for the full grid's density stratification variable.
  real(kind=dp), allocatable, dimension(:, :) :: rho04disk, rho14disk, rho24disk
  real(kind=dp), allocatable, dimension(:, :) :: ux04disk,  ux14disk,  ux24disk
  real(kind=dp), allocatable, dimension(:, :) :: uy04disk,  uy14disk,  uy24disk
  real(kind=dp), allocatable, dimension(:, :) :: uz04disk,  uz14disk,  uz24disk

  ! Dummy dimension variable for writing datasets.  All of the configuration
  ! values written are scalars, though a vector of dimensions must be supplied
  ! to each h5dwrite_f() call.  Note that this value is ignored by the HDF5
  ! library.
  integer(hsize_t), dimension(3)              :: data_dimensions

  ! HDF5 identifiers and status code.
  integer(hid_t)                              :: dataset_id
  integer(hid_t)                              :: group_id
  integer                                     :: hdf5_status

  ! Allocate space for the restart processing variables.
  if ( rank == root ) then
     allocate( ux04disk(1:nsg, 1:nsuby), ux14disk(1:nsg, 1:nsuby), ux24disk(1:nsg, 1:nsuby) )
     allocate( uy04disk(1:nsg, 1:nsuby), uy14disk(1:nsg, 1:nsuby), uy24disk(1:nsg, 1:nsuby) )
     allocate( uz04disk(1:nsg, 1:nsuby), uz14disk(1:nsg, 1:nsuby), uz24disk(1:nsg, 1:nsuby) )
     allocate( rho04disk(1:nsg, 1:nsuby), rho14disk(1:nsg, 1:nsuby), rho24disk(1:nsg, 1:nsuby) )
  else
     allocate( ux04disk(1, 1), ux14disk(1, 1), ux24disk(1, 1) )
     allocate( uy04disk(1, 1), uy14disk(1, 1), uy24disk(1, 1) )
     allocate( uz04disk(1, 1), uz14disk(1, 1), uz24disk(1, 1) )
     allocate( rho04disk(1, 1), rho14disk(1, 1), rho24disk(1, 1) )
  endif

  !
  ! XXX: this is horribly inefficient, both wrt memory and comms, and should
  !      be replaced by parallel I/O ASAP...
  !
  call gather_3D_array( ux0, ux04disk )
  call gather_3D_array( ux1, ux14disk )
  call gather_3D_array( ux2, ux24disk )

  call gather_3D_array( uy0, uy04disk )
  call gather_3D_array( uy1, uy14disk )
  call gather_3D_array( uy2, uy24disk )

  call gather_3D_array( uz0, uz04disk )
  call gather_3D_array( uz1, uz14disk )
  call gather_3D_array( uz2, uz24disk )

  call gather_3D_array( rho0, rho04disk )
  call gather_3D_array( rho1, rho14disk )
  call gather_3D_array( rho2, rho24disk )

  ! Write out the restart-processing solution to disk.
  if ( rank == root ) then

     ! Specify the dimensions of our grid variables.
     data_dimensions(1) = n * nsubz
     data_dimensions(2) = n * nsubx
     data_dimensions(3) = nsuby

     ! Create a group that holds the restart field variables.
     call H5GCREATE_F( restart_file_id, restart_name_field, group_id, hdf5_status )

     ! Create a scalar variable for the timestep size, dt.
     call H5DCREATE_F( group_id, restart_name_dt, H5T_IEEE_F64LE, &
                       restart_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, dt, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar variable for the most recent timestep size,
     ! dt1.
     call H5DCREATE_F( group_id, restart_name_dt1, H5T_IEEE_F64LE, &
                       restart_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, dt1, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar variable for the second most recent timestep size,
     ! dt2.
     call H5DCREATE_F( group_id, restart_name_dt2, H5T_IEEE_F64LE, &
                       restart_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, dt2, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar variable for the simulation time, simulation_time.
     call H5DCREATE_F( group_id, restart_name_time, H5T_IEEE_F64LE, &
                       restart_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, simulation_time, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's most recent perturbation
     ! density, rho04disk.
     call H5DCREATE_F( group_id, restart_name_rho0, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, rho04disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's second most recent
     ! perturbation density, rho14disk.
     call H5DCREATE_F( group_id, restart_name_rho1, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, rho14disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's third most recent
     ! perturbation density, rho24disk.
     call H5DCREATE_F( group_id, restart_name_rho2, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, rho24disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's most recent
     ! horizontal velocity, ux04disk.
     call H5DCREATE_F( group_id, restart_name_ux0, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, ux04disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's second most recent
     ! horizontal velocity, ux14disk.
     call H5DCREATE_F( group_id, restart_name_ux1, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, ux14disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's third most recent
     ! horizontal velocity, ux24disk.
     call H5DCREATE_F( group_id, restart_name_ux2, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, ux24disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's most recent transverse
     ! velocity, uy04disk.
     call H5DCREATE_F( group_id, restart_name_uy0, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, uy04disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's second most recent
     ! transverse velocity, uy14disk.
     call H5DCREATE_F( group_id, restart_name_uy1, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, uy14disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's third most recent
     ! transverse velocity, uy24disk.
     call H5DCREATE_F( group_id, restart_name_uy2, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, uy24disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's most recent vertical
     ! velocity, uz04disk.
     call H5DCREATE_F( group_id, restart_name_uz0, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, uz04disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's second most recent
     ! vertical velocity, uz14disk.
     call H5DCREATE_F( group_id, restart_name_uz1, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, uz14disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's third most recent
     ! vertical velocity, uz24disk.
     call H5DCREATE_F( group_id, restart_name_uz2, H5T_IEEE_F64LE, &
                       restart_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, uz24disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar variable for the timestep index, timestep_index.
     call H5DCREATE_F( group_id, restart_name_time_step, H5T_STD_I32LE, &
                       restart_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_INTEGER, time_ndx , data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Flush the step we just constructed to disk.
     call H5FFLUSH_F( group_id, H5F_SCOPE_LOCAL_F, hdf5_status )

    ! Close the restart group.
    call H5GCLOSE_F( group_id, hdf5_status )

    ! Flush the step we just constructed to disk.
    call H5FFLUSH_F( restart_file_id, H5F_SCOPE_LOCAL_F, hdf5_status )

  endif

  deallocate( rho04disk, rho14disk, rho24disk, &
              ux04disk, ux14disk, ux24disk, &
              uy04disk, uy14disk, uy24disk, &
              uz04disk, uz14disk, uz24disk )

end subroutine write_restart_data

subroutine close_restart_file
! Closes the restart file.  Once closed, no further I/O may be performed on the
! restart file without opening it again.  The timestep dataset identifiers
! (restart_number_steps_dataset_id, restart_wall_time_step_dataset_id,
! and restart_step_numeric_error_dataset_id), the various dataspace identifiers
! (restart_grid_dataspace_id, restart_length1_vector_dataspace_id, and
! restart_scalar_dataspace_id), and the file identifier (restart_file_id) are
! closed.

  use constants, only: rank, root
  use HDF5, only:      hid_t, &
                       h5dclose_f, &
                       h5fclose_f, &
                       h5sclose_f
  use io_restart, only:   restart_file_id, &
                       restart_grid_dataspace_id, &
                       restart_scalar_dataspace_id

  implicit none

  integer :: hdf5_status

  if (rank == root) then

     ! Close the grid, scalar, and vector scalar dataspaces.
     call H5SCLOSE_F( restart_grid_dataspace_id, hdf5_status )
     call H5SCLOSE_F( restart_scalar_dataspace_id, hdf5_status )

     ! Close the file so that all of the pending writes are flushed out.
     call H5FCLOSE_F( restart_file_id, hdf5_status )

  endif
end subroutine close_restart_file
