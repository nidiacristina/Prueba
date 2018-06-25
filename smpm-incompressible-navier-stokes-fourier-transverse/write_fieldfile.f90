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
!       the /field/ group appeared to not contain any sub-groups, nor the
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
!       Currently, the field file is written such that each timestep written
!       lives within its own group.

subroutine open_field_file( field_file_name, success_flag )
! Opens the field file for writing throughout the solver's life time.  If
! the specified field file already exists, it will be overwritten.  Nothing
! is written to the field file.
!
! This routine opens both the scalar and grid data spaces.

  use constants, only: n, nsubx, nsuby, nsubz, rank, root
  use HDF5, only:      H5F_ACC_TRUNC_F, H5S_SCALAR_F, &
                       hid_t, hsize_t, &
                       h5fcreate_f, &
                       h5screate_f, h5screate_simple_f
  use io_field, only:  field_file_id, &
                       field_grid_dataspace_id, field_scalar_dataspace_id

  implicit none

  character(len=*), intent(in)   :: field_file_name
  logical, intent(out)           :: success_flag

  integer                        :: data_rank
  integer(hsize_t), dimension(3) :: data_dimensions
  integer                        :: hdf5_status

  if ( rank == root ) then
     ! Create our output file, overwriting any existing one.
     call H5FCREATE_F( trim ( field_file_name ), &
                       H5F_ACC_TRUNC_F, field_file_id, hdf5_status )

     ! Check if file was created succesfully. If not, stop execution.
     success_flag = (hdf5_status >= 0)
     if (success_flag .eqv. .false. ) then
        return
     end if

     ! Create a scalar dataspace that we'll use to create scalar datatypes
     ! with.  each of the non-grid variables are tied to this dataspace.
     call H5SCREATE_F( H5S_SCALAR_F, field_scalar_dataspace_id, hdf5_status )

     ! Create a 3D dataspace that we'll use for grid-interface variables
     ! (coordinates, velocities, density, etc).
     data_rank          = 3
     data_dimensions(1) = n * nsubz
     data_dimensions(2) = n * nsubx
     data_dimensions(3) = nsuby
     call H5SCREATE_SIMPLE_F( data_rank, data_dimensions, field_grid_dataspace_id, &
                              hdf5_status )

     !GAR (05/04/2018): Our background fields are being set to 3D arrays. This is
     !                  wasteful because the are transverse-independent. Consider
     !                  changing them to 2D and incorporate them into the field
     !                  variables via looping scheme. If such change takes place
     !                  we will need the 2D dataspace back.
!     ! Create a 2D dataspace that we'll use for grid-interface variables that
!     ! are constant in the transverse direction (background currents, etc).
!     data_rank = 2
!     call H5SCREATE_SIMPLE_F( data_rank, data_dimensions, field_grid_2D_dataspace_id, &
!                              hdf5_status )
  end if

end subroutine open_field_file

subroutine write_field_header_grid
! Writes out the solver's grid to the field file's header.

  use constants, only:               n, nsg, nsubx, nsuby, nsubz, rank, root
  use geom, only:                    cx, cz
  use HDF5, only:                    H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER, H5T_STD_I32LE, &
                                     hid_t, hsize_t, &
                                     h5dclose_f, h5dcreate_f, h5dwrite_f, &
                                     h5gclose_f, h5gcreate_f, &
                                     h5screate_f, h5screate_simple_f
  use io_field, only:                field_file_id, &
                                     field_grid_dataspace_id, field_scalar_dataspace_id, &
                                     field_name_grid, field_name_n, field_name_mx, field_name_my, &
                                     field_name_mz, field_name_x, field_name_y, field_name_z
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
     ! full_cz.  Note that these are ignored when specified with our scalar
     ! variables.
     data_dimensions(1) = n * nsubz
     data_dimensions(2) = n * nsubx
     data_dimensions(3) = nsuby

     ! Create the grid group that holds all of the variables we're
     ! writing out.
     call H5GCREATE_F( field_file_id, field_name_grid, group_id, hdf5_status )

     ! Create a scalar for the collocation dimension, n.
     call H5DCREATE_F( field_file_id, field_name_n, H5T_STD_I32LE, &
                       field_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_INTEGER, n, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar for the number of x sub-domains, nsubx.
     call H5DCREATE_F( field_file_id, field_name_mx, H5T_STD_I32LE, &
                       field_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_INTEGER, nsubx, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's x positions,
     ! full_cx.
     call H5DCREATE_F( field_file_id, field_name_x, H5T_IEEE_F64LE, &
                       field_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, reshape( full_cx, (/ n * nsubz, n * nsubx, nsuby /) ), data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar for the number of z sub-domains, nsubz.
     call H5DCREATE_F( field_file_id, field_name_mz, H5T_STD_I32LE, &
                       field_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_INTEGER, nsubz, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's z positions,
     ! full_cz.
     call H5DCREATE_F( field_file_id, field_name_z, H5T_IEEE_F64LE, &
                       field_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, reshape( full_cz, (/ n * nsubz, n * nsubx, nsuby /) ), data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar for the number of y grid points, nsuby.
     call H5DCREATE_F( field_file_id, field_name_my, H5T_STD_I32LE, &
                       field_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_INTEGER, nsuby, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's y positions,
     ! full_cy.
     call H5DCREATE_F( field_file_id, field_name_y, H5T_IEEE_F64LE, &
                       field_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, reshape( full_cy, (/ n * nsubz, n * nsubx, nsuby /) ), data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

  end if

  deallocate( full_cx, full_cy, full_cz )

end subroutine write_field_header_grid

subroutine write_field_variables( timestep_index )
! Writes out the field variables ux, uy, uz, and rho, for a single time.  The
! field is collected onto the root rank and written to the output file.  This
! can only be called after write_field_header() has been called.

  use constants, only:          nsg, nsuby, rank, root
  use field_variables, only:    rho, rho_bar, rho_bar_z, dubcdz, ubc, ux, uy, uz
  use HDF5, only:               H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER, H5T_IEEE_F64LE, H5F_SCOPE_LOCAL_F, &
                                hid_t, hsize_t, H5T_STD_I32LE, &
                                h5dclose_f, h5dcreate_f, h5dwrite_f, &
                                h5fflush_f, &
                                h5gclose_f, h5gcreate_f
  use io_field, only:           field_file_id, &
                                field_grid_dataspace_id, field_name_field_initial_conditions, &
                                field_name_initial_conditions_dubcdz, &
                                field_name_initial_conditions_rho_bar, &
                                field_name_initial_conditions_rho_bar_z, &
                                field_name_initial_conditions_ubc, &
                                field_name_number_steps, field_name_field, &
                                field_name_step_base, field_name_rho, field_name_ux, field_name_uy, field_name_uz, &
                                field_name_time, field_name_time_step, field_scalar_dataspace_id
  use precision, only:          dp
  use timestepping, only:       simulation_time

  implicit none

  integer, intent(in)            :: timestep_index

  ! Define the I/O variables that will hold the full velocity and density.
  real(kind=dp), allocatable, dimension(:, :) :: ux4disk, uy4disk, uz4disk, rho4disk

  ! Define the I/O variables that will hold the background density and its vertical derivative.
  real(kind=dp), allocatable, dimension(:, :) :: rhobar4disk, rhobarz4disk

  ! Define the I/O variables that will hold the background velocity and its vertical derivative.
  real(kind=dp), allocatable, dimension(:, :) :: ubc4disk, dubcdz4disk

  ! Dummy dimension variable for writing datasets.  All of the configuration
  ! values written are scalars, though a vector of dimensions must be supplied
  ! to each h5dwrite_f() call.  Note that this value is ignored by the HDF5
  ! library.
  ! GAR (05/04/2018): A question for Greg, are all calls to h5dwrite_f()
  !                   ignoring data_dimensions or just scalars?
  integer(hsize_t), dimension(1)              :: data_dimensions = (/ 1 /)

  ! HDF5 identifiers and status code.
  integer(hid_t)                              :: dataset_id
  integer(hid_t)                              :: group_id
  integer                                     :: hdf5_status

  ! Allocate arrays for file I/O.
  if ( rank == root ) then
     allocate( rhobar4disk(1:nsg, 1:nsuby), rhobarz4disk(1:nsg, 1:nsuby) )
     allocate( ubc4disk(1:nsg, 1:nsuby), dubcdz4disk(1:nsg, 1:nsuby) )
     allocate( ux4disk(1:nsg, 1:nsuby), uy4disk(1:nsg, 1:nsuby), uz4disk(1:nsg, 1:nsuby), rho4disk(1:nsg, 1:nsuby) )
  else
     allocate( rhobar4disk(1, 1), rhobarz4disk(1, 1) )
     allocate( ubc4disk(1, 1), dubcdz4disk(1, 1) )
     allocate( ux4disk(1, 1), uy4disk(1, 1), uz4disk(1, 1), rho4disk(1, 1) )
  end if

  ! Gather all the field variables onto root.
  call gather_3D_array( ux, ux4disk )
  call gather_3D_array( uy, uy4disk )
  call gather_3D_array( uz, uz4disk )
  call gather_3D_array( rho, rho4disk )

  ! Gather the constant field variables onto root.
  call gather_3D_array( rho_bar, rhobar4disk )
  call gather_3D_array( rho_bar_z, rhobarz4disk )
  call gather_3D_array( ubc, ubc4disk )
  call gather_3D_array( dubcdz, dubcdz4disk )

  if ( rank == root ) then

     ! Write the Field Variables- rho, ux, uy and uz

     ! First, we create the group within the field file called "field". Next, we
     ! proceed to write ux, uy, uz and rho within each "field" group.
     call H5GCREATE_F( field_file_id, field_name_field, group_id, hdf5_status )

     ! Create a scalar grid variable for the grid's horizontal (x) velocity
     ! ux4disk.
     call H5DCREATE_F( group_id, field_name_ux, H5T_IEEE_F64LE, &
                       field_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, ux4disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's transverse (y) velocity
     ! uy4disk.
     call H5DCREATE_F( group_id, field_name_uy, H5T_IEEE_F64LE, &
                       field_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, uy4disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's vertical (z) velocity
     ! uz4disk.
     call H5DCREATE_F( group_id, field_name_uz, H5T_IEEE_F64LE, &
                       field_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, uz4disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's perturbation density (rho)
     ! rho4disk.
     call H5DCREATE_F( group_id, field_name_rho, H5T_IEEE_F64LE, &
                       field_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, rho4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar variable for the current simulation time,
     ! simulation_time.
     call H5DCREATE_F( group_id, field_name_time, H5T_IEEE_F64LE, &
                       field_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, simulation_time, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar variable for the current timestep,
     ! timestep_index.
     call H5DCREATE_F( group_id, field_name_time_step, H5T_STD_I32LE, &
                       field_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_INTEGER, timestep_index, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Flush the step we just constructed to disk.
     call H5FFLUSH_F( group_id, H5F_SCOPE_LOCAL_F, hdf5_status )

     ! Close the "field" group.
     call H5GCLOSE_F( group_id, hdf5_status )


     ! Write the constant field: rho_bar, rho_bar_z, ubc, dubcdz

     ! First, create the group, within the field file, that will hold the contant fields,
     ! step0. The name remains the same as in the previous version of write_fieldfile, just that
     ! we now dump these constant fields at every field file.
     call H5GCREATE_F( field_file_id, field_name_field_initial_conditions, group_id, hdf5_status )

     ! Create a scalar grid variable for the grid's background density profile,
     ! rhobar4disk.
     call H5DCREATE_F( group_id, field_name_initial_conditions_rho_bar, H5T_IEEE_F64LE, &
                       field_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, rhobar4disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's vertical derivative of the background density,
     ! rhobarz4disk.
     call H5DCREATE_F( group_id, field_name_initial_conditions_rho_bar_z, H5T_IEEE_F64LE, &
                       field_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, rhobarz4disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's background velocity,
     ! ubc4disk.
     call H5DCREATE_F( group_id, field_name_initial_conditions_ubc, H5T_IEEE_F64LE, &
                       field_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, ubc4disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's vertical derivative of the background velocity,
     ! dubcdz4disk.
     call H5DCREATE_F( group_id, field_name_initial_conditions_dubcdz, H5T_IEEE_F64LE, &
                       field_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, dubcdz4disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Flush the step we just constructed to disk.
     call H5FFLUSH_F( group_id, H5F_SCOPE_LOCAL_F, hdf5_status )

     ! Close constant fields (step0) group.
     call H5GCLOSE_F( group_id, hdf5_status )

  endif

  deallocate( ux4disk, uy4disk, uz4disk, rho4disk )
  deallocate( rhobar4disk, rhobarz4disk, ubc4disk, dubcdz4disk )

end subroutine write_field_variables

subroutine close_field_file
! Closes the field file.  Once closed, no further I/O may be performed on the
! field file without opening it again.  The timestep dataset identifiers
! (field_number_steps_dataset_id, field_wall_time_step_dataset_id,
! and field_step_numeric_error_dataset_id), the various dataspace identifiers
! (field_grid_dataspace_id, field_length1_vector_dataspace_id, and
! field_scalar_dataspace_id), and the file identifier (field_file_id) are
! closed.
  use constants, only: rank, root
  use HDF5, only:      hid_t, &
                       h5dclose_f, &
                       h5fclose_f, &
                       h5sclose_f
  use io_field, only:  field_file_id, &
                       field_grid_dataspace_id, &
                       field_scalar_dataspace_id

  implicit none

  integer :: hdf5_status

  ! Since the root rank created the field file, only the root rank closes the file
  if (rank == root) then

     ! Close the grid, scalar, and vector scalar dataspaces.
     call H5SCLOSE_F( field_grid_dataspace_id, hdf5_status )

     ! Close the scalar data space (everything that is not an array)
     call H5SCLOSE_F( field_scalar_dataspace_id, hdf5_status )

     ! Close the file so that all of the pending writes are flushed out.
     call H5FCLOSE_F( field_file_id, hdf5_status )

  endif

end subroutine close_field_file

subroutine gather_3D_array( gather_buffer, root_buffer )

  use constants, only:                  nsg, nsuby, rank, root
  use mpi, only:                        MPI_COMM_WORLD, MPI_DOUBLE_PRECISION
  use precision, only:                  dp
  use woodbury_matrices, only:          rpk

  implicit none

  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(inout) :: gather_buffer
  real(kind=dp), dimension(1:nsg, 1:nsuby), intent(out)   :: root_buffer

  integer                                                 :: ierr

  call perfect_shuffle( nsuby, rpk, gather_buffer )
  call MPI_GATHER( gather_buffer, rpk * nsuby, MPI_DOUBLE_PRECISION, &
                     root_buffer, rpk * nsuby, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr )
  call perfect_shuffle( rpk, nsuby, gather_buffer )
  if ( rank == root ) then
     call perfect_shuffle( nsg, nsuby, root_buffer )
  endif

end subroutine gather_3D_array
