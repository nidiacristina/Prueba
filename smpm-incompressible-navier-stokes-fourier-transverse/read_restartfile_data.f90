subroutine read_restartfile_data( restart_file_name, success_flag )
! Reads and parses a restart file, which sets the initial values of all the
! field variables.

  use constants, only:                  n, nsubx, nsuby, nsubz, rank, root
  use field_variables, only:            rho0, rho1, rho2, rho_bar, rho_bar_z, &
                                        ubc, dubcdz, &
                                        ux0, ux1, ux2, &
                                        uy0, uy1, uy2, &
                                        uz0, uz1, uz2, &
                                        ux_b, uy_b, uz_b
  use geom, only:                       cx, cz
  use HDF5, only:                       H5F_ACC_RDWR_F, H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, &
                                        hid_t, hsize_t, &
                                        h5dopen_f, h5dread_f, h5dclose_f, h5fclose_f, h5fopen_f
  use io_restart, only:                 restart_ic_name_rho0, restart_ic_name_rho1, restart_ic_name_rho2, &
                                        restart_ic_name_ux0, restart_ic_name_ux1, restart_ic_name_ux2, &
                                        restart_ic_name_uy0, restart_ic_name_uy1, restart_ic_name_uy2, &
                                        restart_ic_name_uz0, restart_ic_name_uz1, restart_ic_name_uz2, &
                                        restart_ic_name_dt, restart_ic_name_dt1, restart_ic_name_dt2, &
                                        restart_ic_name_time, &
                                        restart_name_dubcdz, &
                                        restart_name_rho_bar, restart_name_rho_bar_z, &
                                        restart_name_ubc, &
                                        restart_name_ux_b, restart_name_uy_b, restart_name_uz_b, &
                                        restart_name_x, restart_name_y, restart_name_z
  use mpi, only:                        MPI_COMM_WORLD, MPI_DOUBLE_PRECISION
  use precision, only:                  dp
  use timestepping, only:               dt, dt1, dt2, simulation_time
  use transverse, only:                 cy

  implicit none

  ! Path to the HDF5 restart file.
  character(len=*), intent(in)                   :: restart_file_name
  logical, intent(out)                           :: success_flag

  ! MPI status variable.
  integer                                        :: ierr

  ! HDF5 file identifier, status variable, and variable dataspace dimensions.
  integer(hid_t)                                 :: file_id = -1
  integer                                        :: hdf5_status
  integer(hsize_t), dimension(3)                 :: data_dimensions

  ! Set up a temporary buffer to read the data into.
  real(kind=dp), allocatable, dimension(:, :, :) :: read_buffer_3D
  real(kind=dp), allocatable, dimension(:, :)    :: read_buffer_2D

  ! Most all of the restart variables are grid-shaped.
  !
  ! NOTE: We rely on the fact that the Y dimension is the slowest as this
  !       vector is passed into routines that read both 3D and 2D buffers.  In
  !       the latter the third dimension (Y) is ignored.
  !
  data_dimensions(1) = n * nsubz
  data_dimensions(2) = n * nsubx
  data_dimensions(3) = nsuby

  ! Allocate buffers for reading and distributing from.
  allocate( read_buffer_2D(1:n*nsubz, 1:n*nsubx)          )
  allocate( read_buffer_3D(1:n*nsubz, 1:n*nsubx, 1:nsuby) )

  ! Root is the only rank to do I/O, so open the file on it.
  if ( rank == root ) then
     call H5FOPEN_F( restart_file_name, H5F_ACC_RDWR_F, file_id, hdf5_status )
  end if

  call sync_flag( hdf5_status )

  ! We assume that since we were able to open the file, everything is fine.
  !
  ! NOTE: Checking each HDF5 operation would *greatly* bloat the code and
  !       isn't likely worth the effort since a non-programmatic failure in
  !       the middle of a series of HDF5 operations is almost certainly
  !       indicative of an environmental problem that is beyond our
  !       control.  Spiraling out of control in a glorious ball of flames
  !       is perfectly acceptable in that situation.
  success_flag = (hdf5_status >= 0)
  if (success_flag .eqv. .false. ) then
     return
  end if

  ! Read amd distribute the scalar fields
  call read_scalar_real_variable( file_id, restart_ic_name_time, simulation_time, rank == root )
  call MPI_BCAST( simulation_time, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr )

  call read_scalar_real_variable( file_id, restart_ic_name_dt, dt, rank == root )
  call MPI_BCAST( dt, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr )

  call read_scalar_real_variable( file_id, restart_ic_name_dt1, dt1, rank == root )
  call MPI_BCAST( dt1, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr )

  call read_scalar_real_variable( file_id, restart_ic_name_dt2, dt2, rank == root )
  call MPI_BCAST( dt2, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr )

  ! Read and distribute the restart conditions from the root rank.
  call read_3D_real_variable( file_id, restart_name_x, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, cx )

  call read_3D_real_variable( file_id, restart_name_y, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, cy )

  call read_3D_real_variable( file_id, restart_name_z, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, cz )

  call read_3D_real_variable( file_id, restart_name_dubcdz, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, dubcdz )

  call read_3D_real_variable( file_id, restart_name_rho_bar, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, rho_bar )

  call read_3D_real_variable( file_id, restart_name_rho_bar_z, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, rho_bar_z )

  call read_3D_real_variable( file_id, restart_ic_name_rho0, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, rho0 )

  call read_3D_real_variable( file_id, restart_ic_name_rho1, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, rho1 )

  call read_3D_real_variable( file_id, restart_ic_name_rho2, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, rho2 )

  call read_3D_real_variable( file_id, restart_name_ubc, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, ubc )

  call read_3D_real_variable( file_id, restart_name_ux_b, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, ux_b )

  call read_3D_real_variable( file_id, restart_ic_name_ux0, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, ux0 )

  call read_3D_real_variable( file_id, restart_ic_name_ux1, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, ux1 )

  call read_3D_real_variable( file_id, restart_ic_name_ux2, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, ux2 )

  call read_3D_real_variable( file_id, restart_name_uy_b, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, uy_b )

  call read_3D_real_variable( file_id, restart_ic_name_uy0, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, uy0 )

  call read_3D_real_variable( file_id, restart_ic_name_uy1, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, uy1 )

  call read_3D_real_variable( file_id, restart_ic_name_uy2, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, uy2 )

  call read_3D_real_variable( file_id, restart_name_uz_b, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, uz_b )

  call read_3D_real_variable( file_id, restart_ic_name_uz0, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, uz0 )

  call read_3D_real_variable( file_id, restart_ic_name_uz1, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, uz1 )

  call read_3D_real_variable( file_id, restart_ic_name_uz2, data_dimensions, read_buffer_3D, rank == root )
  call scatter_3D_array( read_buffer_3D, uz2 )

  ! Close the file.
  if ( rank == root ) then
     call H5FCLOSE_F( file_id, hdf5_status )
  end if

  deallocate( read_buffer_3D, read_buffer_2D )

end subroutine read_restartfile_data

subroutine read_scalar_real_variable( file_id, hdf5_path, scalar_value, do_reading_flag )
! Reads a scalar real-valued variable at the specified path into the provided
! buffer.  Reading only occurs if the do_reading_flag logical is true,
! allowing for rank-specific reading.

  use HDF5, only:                               H5T_IEEE_F64LE, &
                                                hid_t, hsize_t, &
                                                h5dclose_f, h5dopen_f, h5dread_f
  use precision, only:                          dp

  implicit none

  ! Sub-routine parameters.
  integer(hid_t), intent(inout)  :: file_id
  character(len=*), intent(in)   :: hdf5_path
  real(kind=dp), intent(out)     :: scalar_value
  logical, intent(in)            :: do_reading_flag

  ! HDF5 dataset identifier and status variable.
  integer(hid_t)                 :: dataset_id
  integer                        :: hdf5_status

  ! HDF5 "buffer" and dimension to work around the fact that a single
  ! dimensioned variable with length 1 is not compatible with a scalar.
  integer(hsize_t), dimension(1) :: data_dimensions
  real(kind=dp), dimension(1)    :: read_buffer

  data_dimensions(1) = 1

  ! Only read from the HDF5 path if we're requested to.
  if ( do_reading_flag ) then
     call H5DOPEN_F( file_id, hdf5_path, dataset_id, hdf5_status )
     call H5DREAD_F( dataset_id, H5T_IEEE_F64LE, read_buffer, data_dimensions, hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     scalar_value = read_buffer(1)
  end if

end subroutine read_scalar_real_variable

subroutine read_2D_real_variable( file_id, hdf5_path, data_dimensions, read_buffer, do_reading_flag )
! Reads a 2D real-valued variable at the specified path into the provided
! buffer.  Reading only occurs if the do_reading_flag logical is true,
! allowing for rank-specific reading.

  use HDF5, only:                               H5T_IEEE_F64LE, &
                                                hid_t, hsize_t, &
                                                h5dclose_f, h5dopen_f, h5dread_f
  use precision, only:                          dp

  implicit none

  ! Sub-routine parameters.
  integer(hid_t), intent(inout)                               :: file_id
  character(len=*), intent(in)                                :: hdf5_path
  integer(hsize_t), dimension(2), intent(in)                  :: data_dimensions
  real(kind=dp), dimension(1:data_dimensions(1), &
                           1:data_dimensions(2)), intent(out) :: read_buffer
  logical, intent(in)                                         :: do_reading_flag

  ! HDF5 dataset identifier and status variable.
  integer(hid_t)                                              :: dataset_id
  integer                                                     :: hdf5_status

  ! Only read from the HDF5 path if we're requested to.
  if ( do_reading_flag ) then
     call H5DOPEN_F( file_id, hdf5_path, dataset_id, hdf5_status )
     call H5DREAD_F( dataset_id, H5T_IEEE_F64LE, read_buffer, &
                     data_dimensions, hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )
  end if

end subroutine read_2D_real_variable

subroutine read_3D_real_variable( file_id, hdf5_path, data_dimensions, read_buffer, do_reading_flag )
! Reads a 3D real-valued variable at the specified path into the provided
! buffer.  Reading only occurs if the do_reading_flag logical is true,
! allowing for rank-specific reading.

  use HDF5, only:                               H5T_IEEE_F64LE, &
                                                hid_t, hsize_t, &
                                                h5dclose_f, h5dopen_f, h5dread_f
  use precision, only:                          dp

  implicit none

  ! Sub-routine parameters.
  integer(hid_t), intent(inout)                               :: file_id
  character(len=*), intent(in)                                :: hdf5_path
  integer(hsize_t), dimension(3), intent(in)                  :: data_dimensions
  real(kind=dp), dimension(1:data_dimensions(1), &
                           1:data_dimensions(2), &
                           1:data_dimensions(3)), intent(out) :: read_buffer
  logical, intent(in)                                         :: do_reading_flag

  ! HDF5 dataset identifier and status variable.
  integer(hid_t)                                              :: dataset_id
  integer                                                     :: hdf5_status

  ! Only read from the HDF5 path if we're requested to.
  if ( do_reading_flag ) then
     call H5DOPEN_F( file_id, hdf5_path, dataset_id, hdf5_status )
     call H5DREAD_F( dataset_id, H5T_IEEE_F64LE, read_buffer, &
                     data_dimensions, hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )
  end if

end subroutine read_3D_real_variable
