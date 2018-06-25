subroutine read_initfile_data( initial_conditions_file_name, success_flag )
! Reads and parses an initial conditions file, which sets the initial values
! of all the field variables.

  use constants, only:                  n, nsubx, nsuby, nsubz, rank, root
  use field_variables, only:            dubcdz, rho, rho_bar, rho_bar_z, ubc, ux, uy, uz
  use geom, only:                       cx, cz
  use HDF5, only:                       H5F_ACC_RDWR_F, &
                                        hid_t, hsize_t, &
                                        h5fclose_f, h5fopen_f
  use io_field, only:                   ic_name_grid_x, ic_name_grid_y, ic_name_grid_z, &
                                        ic_name_dubcdz, &
                                        ic_name_rho, ic_name_rho_bar, ic_name_rho_bar_z, &
                                        ic_name_ubc, ic_name_ux, ic_name_uy, ic_name_uz
  use precision, only:                  dp
  use transverse, only:                 cy

  implicit none

  ! Path to the HDF5 initial conditions file.
  character(len=*), intent(in)                   :: initial_conditions_file_name
  logical, intent(out)                           :: success_flag

  ! HDF5 file identifier, status variable, and variable dataspace dimensions.
  integer(hid_t)                                 :: file_id = -1
  integer                                        :: hdf5_status
  integer(hsize_t), dimension(3)                 :: data_dimensions

  ! Set up a temporary buffer to read the data into.
  real(kind=dp), allocatable, dimension(:, :, :) :: read_buffer

  ! Most all of the initial conditions are grid-shaped.
  !
  ! NOTE: We rely on the fact that the Y dimension is the slowest as this
  !       vector is passed into routines that read both 3D and 2D buffers.  In
  !       the latter the third dimension (Y) is ignored.
  !
  data_dimensions(1) = n * nsubz
  data_dimensions(2) = n * nsubx
  data_dimensions(3) = nsuby

  ! Allocate a buffer for reading and distributing from.
  allocate( read_buffer(1:n*nsubz, 1:n*nsubx, 1:nsuby) )

  ! Root is the only rank to do I/O, so open the file on it.
  if ( rank == root ) then
     call H5FOPEN_F( initial_conditions_file_name, &
                     H5F_ACC_RDWR_F, file_id, hdf5_status )
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

  ! Read and distribute the initial conditions from the root rank.
  call read_3D_real_variable( file_id, ic_name_grid_x, data_dimensions, read_buffer, rank == root )
  call scatter_3D_array( read_buffer, cx )

  call read_3D_real_variable( file_id, ic_name_grid_y, data_dimensions, read_buffer, rank == root )
  call scatter_3D_array( read_buffer, cy )

  call read_3D_real_variable( file_id, ic_name_grid_z, data_dimensions, read_buffer, rank == root )
  call scatter_3D_array( read_buffer, cz )

  ! Read and distribute the vertical derivative of the background current
  ! XXX: this needs to be extended to 3D - GNT
  call read_3D_real_variable( file_id, ic_name_dubcdz, data_dimensions, read_buffer, rank == root )
  call scatter_3D_array( read_buffer, dubcdz )

  ! Read and distribute the perturbation density field
  call read_3D_real_variable( file_id, ic_name_rho, data_dimensions, read_buffer, rank == root )
  call scatter_3D_array( read_buffer, rho )

  ! Read and distribute the background density profile
  ! XXX: this needs to be extended to 3D - GNT
  call read_3D_real_variable( file_id, ic_name_rho_bar, data_dimensions, read_buffer, rank == root )
  call scatter_3D_array( read_buffer, rho_bar )

  ! Read and distribute the vertical derivative of the background density profile
  ! XXX: this needs to be extended to 3D - GNT
  call read_3D_real_variable( file_id, ic_name_rho_bar_z, data_dimensions, read_buffer, rank == root )
  call scatter_3D_array( read_buffer, rho_bar_z )

  ! Read and distribute the background velocity profile
  ! XXX: this needs to be extended to 3D - GNT
  call read_3D_real_variable( file_id, ic_name_ubc, data_dimensions, read_buffer, rank == root )
  call scatter_3D_array( read_buffer, ubc )

  ! Read and distribute the horizontal perturbation velocity
  call read_3D_real_variable( file_id, ic_name_ux, data_dimensions, read_buffer, rank == root )
  call scatter_3D_array( read_buffer, ux )

  ! Read and distribute the transverse perturbation velocity
  call read_3D_real_variable( file_id, ic_name_uy, data_dimensions, read_buffer, rank == root )
  call scatter_3D_array( read_buffer, uy )

  ! Read and distribute the vertical perturbation velocity
  call read_3D_real_variable( file_id, ic_name_uz, data_dimensions, read_buffer, rank == root )
  call scatter_3D_array( read_buffer, uz )

  ! Close the file.
  if ( rank == root ) then
     call H5FCLOSE_F( file_id, hdf5_status )
  end if

  deallocate( read_buffer )

end subroutine read_initfile_data

subroutine scatter_3D_array( root_buffer, scattered_buffer )

  use constants, only:                  n, nsg, nsubx, nsuby, nsubz, rank, root
  use mpi, only:                        MPI_COMM_WORLD, MPI_DOUBLE_PRECISION
  use precision, only:                  dp
  use woodbury_matrices, only:          rpk

  implicit none

  real(kind=dp), dimension(1:nsg, 1:nsuby), intent(inout) :: root_buffer
  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(out)   :: scattered_buffer
  integer                                                 :: ierr

  if ( rank == root )  call perfect_shuffle( nsuby, n * n * nsubx * nsubz, root_buffer )
  call MPI_SCATTER(      root_buffer, rpk * nsuby, MPI_DOUBLE_PRECISION, &
                    scattered_buffer, rpk * nsuby, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr )
  call perfect_shuffle( rpk, nsuby, scattered_buffer )

end subroutine scatter_3D_array
