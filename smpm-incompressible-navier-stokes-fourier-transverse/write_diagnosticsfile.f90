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

subroutine open_diagnostics_file( diagnostics_file_name, success_flag )
! Opens the diagnostics file for writing throughout the solver's life time.  If
! the specified diagnostics file already exists, it will be overwritten.  Nothing
! is written to the diagnostics file.
!
! This routine opens both the scalar and grid data spaces.

  use constants, only:      rank, root
  use HDF5, only:           H5F_ACC_TRUNC_F, H5S_SCALAR_F, &
                            hid_t, hsize_t, &
                            h5fcreate_f, &
                            h5screate_f, h5screate_simple_f
  use io_diagnostics, only: diagnostics_file_id, &
                            diagnostics_grid_2D_dataspace_id, &
                            diagnostics_scalar_dataspace_id

  implicit none

  character(len=*), intent(in)   :: diagnostics_file_name
  logical, intent(out)           :: success_flag

  integer                        :: data_rank
  integer(hsize_t), dimension(3) :: data_dimensions
  integer                        :: hdf5_status

  if ( rank == root ) then
     ! Create our output file, overwriting any existing one.
     call H5FCREATE_F( trim ( diagnostics_file_name ), &
                       H5F_ACC_TRUNC_F, diagnostics_file_id, hdf5_status )

     ! Check if file was created succesfully. If not, stop execution.
     success_flag = (hdf5_status >= 0)
     if (success_flag .eqv. .false. ) then
        return
     end if

     ! Create a scalar dataspace that we'll use to create scalar datatypes
     ! with.  each of the non-grid variables are tied to this dataspace.
     call H5SCREATE_F( H5S_SCALAR_F, diagnostics_scalar_dataspace_id, hdf5_status )

     ! Create a 2D dataspace that we'll use for grid-interface variables that
     ! pertain to the possion operator (mainly the left and capacitance kernel vectors).
     data_rank = 2
     call H5SCREATE_SIMPLE_F( data_rank, data_dimensions, diagnostics_grid_2D_dataspace_id, &
                              hdf5_status )
  end if

end subroutine open_diagnostics_file

subroutine write_diagnostics_header( solver_start, elapsed_time_null_basis, elapsed_time_null_error, &
                               elapsed_time_setup )
! Writes the diagnostics file's header to disk.  The solver's configuration is
! collected onto the root rank and written to the output file.  This routine
! must be called before write_diagnostics_timestamp() and
! after open_diagnostics_file() is called.  It also must be called after the solver's
! configuration has been finalized since it accesses most of the solver's state
! to write out the diagnostics file's header.
!
! This routine creates the number timesteps dataset identifier.
  use constants, only:               rank, root
  use errors, only:                  l2_poisson_kernel_error, &
                                     l2_poisson_schur_kernel_error
  use HDF5, only:                    H5F_SCOPE_LOCAL_F, H5T_STD_I32LE, &
                                     hid_t, &
                                     h5dcreate_f, &
                                     h5fflush_f, &
                                     h5gclose_f, h5gcreate_f
  use io_diagnostics, only:          diagnostics_file_id
  use precision, only:               dp

  implicit none

  real(kind=dp), intent(in)                :: solver_start
  real(kind=dp), intent(in)                :: elapsed_time_null_basis
  real(kind=dp), intent(in)                :: elapsed_time_null_error
  real(kind=dp), intent(in)                :: elapsed_time_setup

  ! HDF5 identifier and status code.
  integer                                  :: hdf5_status

  ! Write out the configuration.
  call write_diagnostics_header_config

  ! NOTE: The execution configuration is currently the only portion of the
  !       header that the root rank can write without cooperation from other
  !       ranks.
  if ( rank == root ) then
     ! Write out the execution configuration.
     call write_diagnostics_header_execution_config( solver_start, &
                                               elapsed_time_null_basis, &
                                               elapsed_time_null_error, &
                                               elapsed_time_setup, &
                                               l2_poisson_kernel_error, &
                                               l2_poisson_schur_kernel_error )

     ! Flush the contents of the file as they exist so that others see the
     ! current contents.
     call H5FFLUSH_F( diagnostics_file_id, H5F_SCOPE_LOCAL_F, hdf5_status )

  endif

end subroutine write_diagnostics_header

subroutine write_diagnostics_header_config
! Writes out the solver's configuration to the diagnostics file's header.

  use constants, only:         nsg, nu, rank, root
  use HDF5, only:              H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER, H5T_STD_I32LE, &
                               hid_t, hsize_t, &
                               h5dclose_f, h5dcreate_f, h5dwrite_f, &
                               h5gclose_f, h5gcreate_f, &
                               h5sclose_f, h5screate_f, h5screate_simple_f
  use io_diagnostics, only:    diagnostics_file_id, &
                               diagnostics_scalar_dataspace_id, &
                               diagnostics_name_config, &
                               diagnostics_name_uL, diagnostics_name_uC, diagnostics_name_facrobin, &
                               diagnostics_name_facrobin_ppe, diagnostics_name_nu, diagnostics_name_filter_order_xz, &
                               diagnostics_name_filter_order_y, diagnostics_name_dt, diagnostics_name_tend, &
                               diagnostics_name_timesteps_per_write
  use mpi, only:               MPI_COMM_WORLD, MPI_DOUBLE_PRECISION
  use options, only:           facrobin, facrobin_ppe, filter_order_xz, filter_order_y, &
                               timesteps_between_writes
  use precision, only:         dp
  use timestepping, only:      dt, tend
  use woodbury_matrices, only: dimCk, displacements, k, recv_counts, rpk, uC, uL

  implicit none

  integer                                  :: data_rank
  integer(hsize_t), dimension(1)           :: data_dimensions

  ! HDF5 identifiers and status code.
  integer(hid_t)                           :: dataset_id
  integer(hid_t)                           :: group_id
  integer(hid_t)                           :: matrix_dataspace_id
  integer                                  :: hdf5_status

  ! Arrays to gather the kernel vectors.
  real(kind=dp), allocatable, dimension(:) :: uC_all, uL_all

  integer                                  :: ierr

  ! Allocate space to collect uC and uL on the root rank.
  if ( rank == root ) then
     allocate( uC_all(1:k), uL_all(1:nsg) )
  else
     allocate( uC_all(1), uL_all(1) )
  end if

  ! Pull the pieces of uC and uL onto the root rank.
  !
  ! NOTE: uC is distributed unevenly so that the ranks are partitioned into 2
  !       groups:
  !
  !   1) edges    - The 1st and last ranks only have one vertical interface
  !                 within the grid.
  !   2) interior - The other (nprocs - 2) ranks have two vertical interfaces
  !                 within the grid.
  !
  !      with the uneven distribution, the vector version of the gather,
  !      MPI_GATHERV(), needs to be used instead of the normal gather,
  !      MPI_GATHER().
  call MPI_GATHERV( uC, dimCk, MPI_DOUBLE_PRECISION, &
                    uC_all, recv_counts, displacements, MPI_DOUBLE_PRECISION, &
                    root, MPI_COMM_WORLD, ierr )

  call MPI_GATHER( uL, rpk, MPI_DOUBLE_PRECISION, &
                   uL_all, rpk, MPI_DOUBLE_PRECISION, &
                   root, MPI_COMM_WORLD, ierr )

  if ( rank == root ) then

     ! Create the group to hold our configuration variables.
     call H5GCREATE_F( diagnostics_file_id, diagnostics_name_config, group_id, hdf5_status )
     call H5GCLOSE_F( group_id, hdf5_status )

     ! Create a scalar for the time step size, dt.
     call H5DCREATE_F( diagnostics_file_id, diagnostics_name_dt, H5T_IEEE_F64LE, &
                       diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, dt, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar for the total simulation time, tend.
     call H5DCREATE_F( diagnostics_file_id, diagnostics_name_tend, H5T_IEEE_F64LE, &
                       diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, tend, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar for the number of time steps per field update,
     ! timesteps_per_write.
     call H5DCREATE_F( diagnostics_file_id, diagnostics_name_timesteps_per_write, H5T_STD_I32LE, &
                       diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_INTEGER, timesteps_between_writes, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar for the facrobin and facrobin PPE parameters,
     ! facrobin and facrobin_ppe, respectively.
     call H5DCREATE_F( diagnostics_file_id, diagnostics_name_facrobin, H5T_IEEE_F64LE, &
                       diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, facrobin, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )
     call H5DCREATE_F( diagnostics_file_id, diagnostics_name_facrobin_ppe, H5T_IEEE_F64LE, &
                       diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, facrobin_ppe, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar for the viscosity, nu.
     call H5DCREATE_F( diagnostics_file_id, diagnostics_name_nu, H5T_IEEE_F64LE, &
                       diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, nu, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar for the filter order in the x-z plane, filter_order_xz.
     call H5DCREATE_F( diagnostics_file_id, diagnostics_name_filter_order_xz, H5T_IEEE_F64LE, &
                       diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, filter_order_xz, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar for the filter order in the y direction, filter_order_y.
     call H5DCREATE_F( diagnostics_file_id, diagnostics_name_filter_order_y, H5T_IEEE_F64LE, &
                       diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, filter_order_y, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid vector for the Poisson operator's capacitance
     ! matrix's left kernel vector, uC.
     data_rank          = 1
     data_dimensions(1) = k
     call H5SCREATE_SIMPLE_F( data_rank, data_dimensions, matrix_dataspace_id, hdf5_status )
     call H5DCREATE_F( diagnostics_file_id, diagnostics_name_uC, H5T_IEEE_F64LE, &
                       matrix_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, uC_all, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )
     call H5SCLOSE_F( matrix_dataspace_id, hdf5_status )

     ! Create a scalar grid vector for the Poisson operator's left kernel
     ! vector, uL.
     data_rank          = 1
     data_dimensions(1) = nsg
     call H5SCREATE_SIMPLE_F( data_rank, data_dimensions, matrix_dataspace_id, hdf5_status )
     call H5DCREATE_F( diagnostics_file_id, diagnostics_name_uL, H5T_IEEE_F64LE, &
                       matrix_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, uL_all, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )
     call H5SCLOSE_F( matrix_dataspace_id, hdf5_status )

     ! write out the GMRES configuration parameters.
     call write_diagnostics_header_gmres

  end if

  deallocate( uC_all, uL_all )

end subroutine write_diagnostics_header_config

subroutine write_diagnostics_header_gmres
! Writes out the solver's GMRES configuration to the diagnostics file's header.

  use HDF5, only:           H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER, H5T_IEEE_F64LE, H5T_STD_I32LE, &
                            hid_t, hsize_t, &
                            h5gclose_f, h5gcreate_f, &
                            h5dclose_f, h5dcreate_f, h5dwrite_f
  use io_diagnostics, only: diagnostics_file_id, &
                            diagnostics_scalar_dataspace_id, &
                            diagnostics_name_config_gmres, &
                            diagnostics_name_poisson_tolerance, diagnostics_name_poisson_max_iters, &
                            diagnostics_name_poisson_restart, diagnostics_name_viscous_tolerance, &
                            diagnostics_name_viscous_max_iters
  use options, only:        gmres_maxit_poisson, gmres_restart_poisson, gmres_tol_poisson, &
                            gmres_maxit_viscous, gmres_tol_viscous

  implicit none

  ! Dummy dimension variable for writing datasets.  All of the configuration
  ! values written are scalars, though a vector of dimensions must be supplied
  ! to each h5dwrite_f() call.  Note that this value is ignored by the HDF5
  ! library.
  integer(hsize_t), dimension(1) :: data_dimensions = (/ 1 /)

  ! HDF5 identifiers and status code.
  integer(hid_t)                 :: dataset_id
  integer(hid_t)                 :: group_id
  integer                        :: hdf5_status

  ! Create the GMRES configuration group that holds all of the variables we're
  ! writing out.
  call H5GCREATE_F( diagnostics_file_id, diagnostics_name_config_gmres, group_id, hdf5_status )
  call H5GCLOSE_F( group_id, hdf5_status )

  ! Create scalars for the Poisson GMRES parameters, convergence tolerance,
  ! maximum iterations, and restart threshold
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_poisson_tolerance, H5T_IEEE_F64LE, &
                    diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_DOUBLE, gmres_tol_poisson, data_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )

  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_poisson_max_iters, H5T_STD_I32LE, &
                    diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_INTEGER, gmres_maxit_poisson, data_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )

  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_poisson_restart, H5T_STD_I32LE, &
                    diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_INTEGER, gmres_restart_poisson, data_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )

  ! Create scalars for the viscous GMRES parameters, convergence tolerance,
  ! and maximum iterations.
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_viscous_tolerance, H5T_IEEE_F64LE, &
                    diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_DOUBLE, gmres_tol_viscous, data_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )

  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_viscous_max_iters, H5T_STD_I32LE, &
                    diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_INTEGER, gmres_maxit_viscous, data_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )

end subroutine write_diagnostics_header_gmres

subroutine write_diagnostics_header_execution_config( solver_start, elapsed_time_null_basis, &
                                                elapsed_time_null_error, elapsed_time_setup, &
                                                l2_poisson_kernel_error, l2_poisson_schur_kernel_error )
! Writes out the solver's execution configuration to the diagnostics file's header.
!
! The wall time step dataset identifier is created in this routine.

  use constants, only:      nprocs
  use HDF5, only:           H5P_DATASET_CREATE_F, H5S_UNLIMITED_F, &
                            H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER, H5T_STD_I32LE, &
                            hid_t, hsize_t, &
                            h5dclose_f, h5dcreate_f, h5dwrite_f, &
                            h5gclose_f, h5gcreate_f, &
                            h5pclose_f, h5pcreate_f, h5pset_chunk_f, &
                            h5sclose_f, h5screate_f, h5screate_simple_f
  use io_diagnostics, only: diagnostics_file_id, &
                            diagnostics_length1_vector_dataspace_id, diagnostics_scalar_dataspace_id, &
                            diagnostics_name_cfl, &
                            diagnostics_name_cfl_x, &
                            diagnostics_name_cfl_y, &
                            diagnostics_name_cfl_z, &
                            diagnostics_cfl_x_dataset_id, &
                            diagnostics_cfl_y_dataset_id, &
                            diagnostics_cfl_z_dataset_id, &
                            diagnostics_time_step_size_dataset_id, &
                            diagnostics_simulation_time_dataset_id, &
                            diagnostics_name_execution, diagnostics_name_number_ranks, diagnostics_name_number_threads, &
                            diagnostics_name_error, &
                            diagnostics_name_cpu_time, &
                            diagnostics_name_gmres_diffusion_real_iters, diagnostics_name_gmres_diffusion_imag_iters, &
                            diagnostics_name_gmres_poisson_real_iters, diagnostics_name_gmres_poisson_imag_iters, &
                            diagnostics_name_gmres_viscous_real_iters, diagnostics_name_gmres_viscous_imag_iters, &
                            diagnostics_name_l2_poisson_error, &
                            diagnostics_name_l2_poisson_kernel_error, &
                            diagnostics_name_l2_poisson_schur_real_error, &
                            diagnostics_name_l2_poisson_schur_imag_error, &
                            diagnostics_name_l2_poisson_schur_kernel_error, &
                            diagnostics_name_linf_diffusion_real_error, diagnostics_name_linf_diffusion_imag_error, &
                            diagnostics_name_linf_divergence_error, &
                            diagnostics_name_linf_viscous_real_error, diagnostics_name_linf_viscous_imag_error, &
                            diagnostics_name_null_basis_wall_time, &
                            diagnostics_name_null_error_wall_time, &
                            diagnostics_name_step_wall_time, &
                            diagnostics_name_setup_wall_time, &
                            diagnostics_name_simulation_time, &
                            diagnostics_name_start_time, &
                            diagnostics_name_time_step_size, &
                            diagnostics_name_total_wall_time, &
                            diagnostics_name_wall_time, &
                            diagnostics_gmres_diffusion_iterations_real_dataset_id, &
                            diagnostics_gmres_diffusion_iterations_imag_dataset_id, &
                            diagnostics_gmres_poisson_iterations_real_dataset_id, &
                            diagnostics_gmres_poisson_iterations_imag_dataset_id, &
                            diagnostics_gmres_viscous_iterations_real_dataset_id, &
                            diagnostics_gmres_viscous_iterations_imag_dataset_id, &
                            diagnostics_l2_poisson_error_dataset_id, &
                            diagnostics_l2_poisson_schur_real_error_dataset_id, &
                            diagnostics_l2_poisson_schur_imag_error_dataset_id, &
                            diagnostics_linf_diffusion_real_error_dataset_id, &
                            diagnostics_linf_diffusion_imag_error_dataset_id, &
                            diagnostics_linf_divergence_error_dataset_id, &
                            diagnostics_linf_viscous_real_error_dataset_id, &
                            diagnostics_linf_viscous_imag_error_dataset_id, &
                            diagnostics_wall_time_step_dataset_id
!$ use omp_lib, only:       omp_get_max_threads
  use options, only:        timesteps_between_writes
  use precision, only:      dp

  implicit none

  real(kind=dp), intent(in)               :: solver_start
  real(kind=dp), intent(in)               :: elapsed_time_null_basis
  real(kind=dp), intent(in)               :: elapsed_time_null_error
  real(kind=dp), intent(in)               :: elapsed_time_setup
  real(kind=dp), intent(in)               :: l2_poisson_kernel_error
  real(kind=dp), intent(in)               :: l2_poisson_schur_kernel_error

  ! Number of OpenMP threads.
  integer                                 :: number_threads

  ! Rank and dimension vectors for the time step vector.
  integer                                 :: data_rank
  integer(hsize_t), dimension(1)          :: current_dimensions
  integer(hsize_t), dimension(1)          :: maximum_dimensions

  ! HDF5 identifiers and status code.
  integer(hid_t)                          :: dataset_id
  integer(hid_t)                          :: chunked_property_id
  integer(hid_t)                          :: extensible_dataspace_id
  integer(hid_t)                          :: group_id
  integer                                 :: hdf5_status

  ! Identify the maximum number of threads OpenMP can execute with.
  number_threads    = 1
  !$ number_threads = omp_get_max_threads()

  ! Create groups to hold our execution measurements.
  call H5GCREATE_F( diagnostics_file_id, diagnostics_name_execution, group_id, hdf5_status )
  call H5GCLOSE_F( group_id, hdf5_status )
  call H5GCREATE_F( diagnostics_file_id, diagnostics_name_error, group_id, hdf5_status )
  call H5GCLOSE_F( group_id, hdf5_status )
  call H5GCREATE_F( diagnostics_file_id, diagnostics_name_cfl, group_id, hdf5_status )
  call H5GCLOSE_F( group_id, hdf5_status )
  call H5GCREATE_F( diagnostics_file_id, diagnostics_name_wall_time, group_id, hdf5_status )
  call H5GCLOSE_F( group_id, hdf5_status )

  ! Create a scalar for the number of MPI ranks, nprocs.
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_number_ranks, H5T_STD_I32LE, &
                    diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_INTEGER, nprocs, current_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )

  ! Create a scalar for the number of OpenMP threads, number_threads.
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_number_threads, H5T_STD_I32LE, &
                    diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_INTEGER, number_threads, current_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )

  ! Create a 1D dataspace with a fixed length of 1 that can grow to
  ! arbitrary lengths.
  data_rank             = 1
  current_dimensions(1) = 1
  maximum_dimensions(1) = H5S_UNLIMITED_F
  call H5SCREATE_SIMPLE_F( data_rank, current_dimensions, extensible_dataspace_id, &
                           hdf5_status, maximum_dimensions )

  ! Create a property that specifies our chunk size for the 1D vector
  ! dataspace.  We size each chunk to correspond to the number of timesteps
  ! between diagnostics updates (this decision is fairly arbitrary).
  data_rank             = 1
  current_dimensions(1) = timesteps_between_writes
  call H5PCREATE_F( H5P_DATASET_CREATE_F, chunked_property_id, hdf5_status )
  call H5PSET_CHUNK_F( chunked_property_id, data_rank, current_dimensions, hdf5_status )

  ! Create datasets with the 1D dataspace that we will extend in the future
  ! (time, GMRES iterations, and numeric errors per step).  These datasets
  ! have to be chunked to support extension.
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_gmres_diffusion_real_iters, H5T_STD_I32LE, &
                    extensible_dataspace_id, diagnostics_gmres_diffusion_iterations_real_dataset_id, hdf5_status, &
                    chunked_property_id )
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_gmres_diffusion_imag_iters, H5T_STD_I32LE, &
                    extensible_dataspace_id, diagnostics_gmres_diffusion_iterations_imag_dataset_id, hdf5_status, &
                    chunked_property_id )
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_gmres_poisson_real_iters, H5T_STD_I32LE, &
                    extensible_dataspace_id, diagnostics_gmres_poisson_iterations_real_dataset_id, hdf5_status, &
                    chunked_property_id )
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_gmres_poisson_imag_iters, H5T_STD_I32LE, &
                    extensible_dataspace_id, diagnostics_gmres_poisson_iterations_imag_dataset_id, hdf5_status, &
                    chunked_property_id )
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_gmres_viscous_real_iters, H5T_STD_I32LE, &
                    extensible_dataspace_id, diagnostics_gmres_viscous_iterations_real_dataset_id, hdf5_status, &
                    chunked_property_id )
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_gmres_viscous_imag_iters, H5T_STD_I32LE, &
                    extensible_dataspace_id, diagnostics_gmres_viscous_iterations_imag_dataset_id, hdf5_status, &
                    chunked_property_id )
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_l2_poisson_error, H5T_IEEE_F64LE, &
                    extensible_dataspace_id, diagnostics_l2_poisson_error_dataset_id, hdf5_status, &
                    chunked_property_id )
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_l2_poisson_schur_real_error, H5T_IEEE_F64LE, &
                    extensible_dataspace_id, diagnostics_l2_poisson_schur_real_error_dataset_id, hdf5_status, &
                    chunked_property_id )
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_l2_poisson_schur_imag_error, H5T_IEEE_F64LE, &
                    extensible_dataspace_id, diagnostics_l2_poisson_schur_imag_error_dataset_id, hdf5_status, &
                    chunked_property_id )
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_linf_diffusion_real_error, H5T_IEEE_F64LE, &
                    extensible_dataspace_id, diagnostics_linf_diffusion_real_error_dataset_id, hdf5_status, &
                    chunked_property_id )
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_linf_diffusion_imag_error, H5T_IEEE_F64LE, &
                    extensible_dataspace_id, diagnostics_linf_diffusion_imag_error_dataset_id, hdf5_status, &
                    chunked_property_id )
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_linf_divergence_error, H5T_IEEE_F64LE, &
                    extensible_dataspace_id, diagnostics_linf_divergence_error_dataset_id, hdf5_status, &
                    chunked_property_id )
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_linf_viscous_real_error, H5T_IEEE_F64LE, &
                    extensible_dataspace_id, diagnostics_linf_viscous_real_error_dataset_id, hdf5_status, &
                    chunked_property_id )
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_linf_viscous_imag_error, H5T_IEEE_F64LE, &
                    extensible_dataspace_id, diagnostics_linf_viscous_imag_error_dataset_id, hdf5_status, &
                    chunked_property_id )
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_step_wall_time, H5T_IEEE_F64LE, &
                    extensible_dataspace_id, diagnostics_wall_time_step_dataset_id, hdf5_status, &
                    chunked_property_id )
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_cfl_x, H5T_IEEE_F64LE, &
                    extensible_dataspace_id, diagnostics_cfl_x_dataset_id, hdf5_status, &
                    chunked_property_id )
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_cfl_y, H5T_IEEE_F64LE, &
                    extensible_dataspace_id, diagnostics_cfl_y_dataset_id, hdf5_status, &
                    chunked_property_id )
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_cfl_z, H5T_IEEE_F64LE, &
                    extensible_dataspace_id, diagnostics_cfl_z_dataset_id, hdf5_status, &
                    chunked_property_id )
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_time_step_size, H5T_IEEE_F64LE, &
                    extensible_dataspace_id, diagnostics_time_step_size_dataset_id, hdf5_status, &
                    chunked_property_id )
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_simulation_time, H5T_IEEE_F64LE, &
                    extensible_dataspace_id, diagnostics_simulation_time_dataset_id, hdf5_status, &
                    chunked_property_id )

  ! Create a 1x1 memory space used to extend the 1D dataspace.
  data_rank             = 1
  current_dimensions(1) = 1
  call H5SCREATE_SIMPLE_F( data_rank, current_dimensions, diagnostics_length1_vector_dataspace_id, hdf5_status )

  ! Close the extensible dataspace.
  call H5PCLOSE_F( chunked_property_id, hdf5_status )
  call H5SCLOSE_F( extensible_dataspace_id, hdf5_status )

  ! Create a scalar for the simulation's start time, in unix seconds.
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_start_time, H5T_IEEE_F64LE, &
                    diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_DOUBLE, solver_start, current_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )

  ! Create a scalar for the setup wall time elapsed.
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_setup_wall_time, H5T_IEEE_F64LE, &
                    diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_DOUBLE, elapsed_time_setup, current_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )

  ! Create a scalar for the null basis setup wall time elapsed.
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_null_basis_wall_time, H5T_IEEE_F64LE, &
                    diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_DOUBLE, elapsed_time_null_basis, current_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )

  ! Create a scalar for the null error wall time elapsed.
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_null_error_wall_time, H5T_IEEE_F64LE, &
                    diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_DOUBLE, elapsed_time_null_error, current_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )

  ! Create scalars for the Poisson and Poisson Schur kernels error.
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_l2_poisson_kernel_error, H5T_IEEE_F64LE, &
                    diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_DOUBLE, l2_poisson_kernel_error, current_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_l2_poisson_schur_kernel_error, H5T_IEEE_F64LE, &
                    diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_DOUBLE, l2_poisson_schur_kernel_error, current_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )


end subroutine write_diagnostics_header_execution_config

subroutine write_diagnostics_timestep( timestep_index, elapsed_time )
! Writes out timing information for a single timestep.

  use errors, only:         gmres_diffusion_iterations_real, gmres_diffusion_iterations_imag, &
                            gmres_poisson_iterations_real, gmres_poisson_iterations_imag, &
                            gmres_viscous_iterations_real, gmres_viscous_iterations_imag, &
                            linf_diffusion_error, linf_divergence_error, &
                            l2_poisson_error, l2_poisson_schur_error, &
                            linf_viscous_error
  use HDF5, only:           H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER, H5F_SCOPE_LOCAL_F, H5S_SELECT_SET_F, &
                            hid_t, hsize_t, &
                            h5dget_space_f, h5dset_extent_f, h5dwrite_f, &
                            h5fflush_f, &
                            h5sselect_hyperslab_f
  use io_diagnostics, only: diagnostics_cfl_x_dataset_id,&
                            diagnostics_cfl_y_dataset_id,&
                            diagnostics_cfl_z_dataset_id,&
                            diagnostics_gmres_diffusion_iterations_real_dataset_id, &
                            diagnostics_gmres_diffusion_iterations_imag_dataset_id, &
                            diagnostics_gmres_poisson_iterations_real_dataset_id, &
                            diagnostics_gmres_poisson_iterations_imag_dataset_id, &
                            diagnostics_gmres_viscous_iterations_real_dataset_id, &
                            diagnostics_gmres_viscous_iterations_imag_dataset_id, &
                            diagnostics_l2_poisson_error_dataset_id, &
                            diagnostics_l2_poisson_schur_real_error_dataset_id, &
                            diagnostics_l2_poisson_schur_imag_error_dataset_id, &
                            diagnostics_length1_vector_dataspace_id, &
                            diagnostics_linf_diffusion_real_error_dataset_id, &
                            diagnostics_linf_diffusion_imag_error_dataset_id, &
                            diagnostics_linf_divergence_error_dataset_id, &
                            diagnostics_linf_viscous_real_error_dataset_id, &
                            diagnostics_linf_viscous_imag_error_dataset_id, &
                            diagnostics_simulation_time_dataset_id, &
                            diagnostics_time_step_size_dataset_id,&
                            diagnostics_wall_time_step_dataset_id
  use options, only:        timesteps_between_writes
  use precision, only:      dp

  use timestepping, only:   cfl_max_x, cfl_max_y, cfl_max_z, simulation_time, dt

  implicit none

  integer, intent(in)            :: timestep_index
  real(kind=dp), intent(in)      :: elapsed_time

  ! HDF5 identifiers and status code.
  integer(hsize_t), dimension(1) :: data_dimensions
  integer(hsize_t), dimension(1) :: count, offset
  integer(hid_t)                 :: extensible_dataspace_id
  integer                        :: hdf5_status

  ! Extend the datasets' extent by one element to include this timestep.
  data_dimensions(1) = timestep_index
  call H5DSET_EXTENT_F( diagnostics_gmres_diffusion_iterations_real_dataset_id, data_dimensions, hdf5_status )
  call H5DSET_EXTENT_F( diagnostics_gmres_diffusion_iterations_imag_dataset_id, data_dimensions, hdf5_status )
  call H5DSET_EXTENT_F( diagnostics_gmres_poisson_iterations_real_dataset_id, data_dimensions, hdf5_status )
  call H5DSET_EXTENT_F( diagnostics_gmres_poisson_iterations_imag_dataset_id, data_dimensions, hdf5_status )
  call H5DSET_EXTENT_F( diagnostics_gmres_viscous_iterations_real_dataset_id, data_dimensions, hdf5_status )
  call H5DSET_EXTENT_F( diagnostics_gmres_viscous_iterations_imag_dataset_id, data_dimensions, hdf5_status )
  call H5DSET_EXTENT_F( diagnostics_l2_poisson_error_dataset_id, data_dimensions, hdf5_status )
  call H5DSET_EXTENT_F( diagnostics_l2_poisson_schur_real_error_dataset_id, data_dimensions, hdf5_status )
  call H5DSET_EXTENT_F( diagnostics_l2_poisson_schur_imag_error_dataset_id, data_dimensions, hdf5_status )
  call H5DSET_EXTENT_F( diagnostics_linf_divergence_error_dataset_id, data_dimensions, hdf5_status )
  call H5DSET_EXTENT_F( diagnostics_linf_diffusion_real_error_dataset_id, data_dimensions, hdf5_status )
  call H5DSET_EXTENT_F( diagnostics_linf_diffusion_imag_error_dataset_id, data_dimensions, hdf5_status )
  call H5DSET_EXTENT_F( diagnostics_linf_viscous_real_error_dataset_id, data_dimensions, hdf5_status )
  call H5DSET_EXTENT_F( diagnostics_linf_viscous_imag_error_dataset_id, data_dimensions, hdf5_status )
  call H5DSET_EXTENT_F( diagnostics_wall_time_step_dataset_id, data_dimensions, hdf5_status )
  call H5DSET_EXTENT_F( diagnostics_cfl_x_dataset_id, data_dimensions, hdf5_status )
  call H5DSET_EXTENT_F( diagnostics_cfl_y_dataset_id, data_dimensions, hdf5_status )
  call H5DSET_EXTENT_F( diagnostics_cfl_z_dataset_id, data_dimensions, hdf5_status )
  call H5DSET_EXTENT_F( diagnostics_time_step_size_dataset_id, data_dimensions, hdf5_status )
  call H5DSET_EXTENT_F( diagnostics_simulation_time_dataset_id, data_dimensions, hdf5_status )

  ! Set the dataspace to a single scalar after the last timestep, and write
  ! out the new values.
  offset(1) = timestep_index - 1
  count(1)  = 1

  ! GMRES variables.
  call H5DGET_SPACE_F( diagnostics_gmres_diffusion_iterations_real_dataset_id, extensible_dataspace_id, hdf5_status )
  call H5SSELECT_HYPERSLAB_F( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call H5DWRITE_F( diagnostics_gmres_diffusion_iterations_real_dataset_id, H5T_NATIVE_INTEGER, gmres_diffusion_iterations_real, &
                   count, hdf5_status, diagnostics_length1_vector_dataspace_id, extensible_dataspace_id )
  call H5DGET_SPACE_F( diagnostics_gmres_diffusion_iterations_imag_dataset_id, extensible_dataspace_id, hdf5_status )
  call H5SSELECT_HYPERSLAB_F( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call H5DWRITE_F( diagnostics_gmres_diffusion_iterations_imag_dataset_id, H5T_NATIVE_INTEGER, gmres_diffusion_iterations_imag, &
                   count, hdf5_status, diagnostics_length1_vector_dataspace_id, extensible_dataspace_id )
  call H5DGET_SPACE_F( diagnostics_gmres_poisson_iterations_real_dataset_id, extensible_dataspace_id, hdf5_status )
  call H5SSELECT_HYPERSLAB_F( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call H5DWRITE_F( diagnostics_gmres_poisson_iterations_real_dataset_id, H5T_NATIVE_INTEGER, gmres_poisson_iterations_real, &
                   count, hdf5_status, diagnostics_length1_vector_dataspace_id, extensible_dataspace_id )
  call H5DGET_SPACE_F( diagnostics_gmres_poisson_iterations_imag_dataset_id, extensible_dataspace_id, hdf5_status )
  call H5SSELECT_HYPERSLAB_F( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call H5DWRITE_F( diagnostics_gmres_poisson_iterations_imag_dataset_id, H5T_NATIVE_INTEGER, gmres_poisson_iterations_imag, &
                   count, hdf5_status, diagnostics_length1_vector_dataspace_id, extensible_dataspace_id )
  call H5DGET_SPACE_F( diagnostics_gmres_viscous_iterations_real_dataset_id, extensible_dataspace_id, hdf5_status )
  call H5SSELECT_HYPERSLAB_F( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call H5DWRITE_F( diagnostics_gmres_viscous_iterations_real_dataset_id, H5T_NATIVE_INTEGER, gmres_viscous_iterations_real, &
                   count, hdf5_status, diagnostics_length1_vector_dataspace_id, extensible_dataspace_id )
  call H5DGET_SPACE_F( diagnostics_gmres_viscous_iterations_imag_dataset_id, extensible_dataspace_id, hdf5_status )
  call H5SSELECT_HYPERSLAB_F( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call H5DWRITE_F( diagnostics_gmres_viscous_iterations_imag_dataset_id, H5T_NATIVE_INTEGER, gmres_viscous_iterations_imag, &
                   count, hdf5_status, diagnostics_length1_vector_dataspace_id, extensible_dataspace_id )

  ! L2 Poisson errors.
  call H5DGET_SPACE_F( diagnostics_l2_poisson_error_dataset_id, extensible_dataspace_id, hdf5_status )
  call H5SSELECT_HYPERSLAB_F( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call H5DWRITE_F( diagnostics_l2_poisson_error_dataset_id, H5T_NATIVE_DOUBLE, l2_poisson_error, &
                   count, hdf5_status, diagnostics_length1_vector_dataspace_id, extensible_dataspace_id )
  call H5DGET_SPACE_F( diagnostics_l2_poisson_schur_real_error_dataset_id, extensible_dataspace_id, hdf5_status )
  call H5SSELECT_HYPERSLAB_F( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call H5DWRITE_F( diagnostics_l2_poisson_schur_real_error_dataset_id, H5T_NATIVE_DOUBLE, real( l2_poisson_schur_error, kind=dp ), &
                   count, hdf5_status, diagnostics_length1_vector_dataspace_id, extensible_dataspace_id )
  call H5DGET_SPACE_F( diagnostics_l2_poisson_schur_imag_error_dataset_id, extensible_dataspace_id, hdf5_status )
  call H5SSELECT_HYPERSLAB_F( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call H5DWRITE_F( diagnostics_l2_poisson_schur_imag_error_dataset_id, H5T_NATIVE_DOUBLE, aimag( l2_poisson_schur_error ), &
                   count, hdf5_status, diagnostics_length1_vector_dataspace_id, extensible_dataspace_id )

  ! Linf errors for divergence, viscosity, and diffusion.
  call H5DGET_SPACE_F( diagnostics_linf_divergence_error_dataset_id, extensible_dataspace_id, hdf5_status )
  call H5SSELECT_HYPERSLAB_F( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call H5DWRITE_F( diagnostics_linf_divergence_error_dataset_id, H5T_NATIVE_DOUBLE, linf_divergence_error, &
                   count, hdf5_status, diagnostics_length1_vector_dataspace_id, extensible_dataspace_id )
  call H5DGET_SPACE_F( diagnostics_linf_viscous_real_error_dataset_id, extensible_dataspace_id, hdf5_status )
  call H5SSELECT_HYPERSLAB_F( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call H5DWRITE_F( diagnostics_linf_viscous_real_error_dataset_id, H5T_NATIVE_DOUBLE, real( linf_viscous_error, kind=dp ), &
                   count, hdf5_status, diagnostics_length1_vector_dataspace_id, extensible_dataspace_id )
  call H5DGET_SPACE_F( diagnostics_linf_viscous_imag_error_dataset_id, extensible_dataspace_id, hdf5_status )
  call H5SSELECT_HYPERSLAB_F( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call H5DWRITE_F( diagnostics_linf_viscous_imag_error_dataset_id, H5T_NATIVE_DOUBLE, aimag( linf_viscous_error ) , &
                   count, hdf5_status, diagnostics_length1_vector_dataspace_id, extensible_dataspace_id )
  call H5DGET_SPACE_F( diagnostics_linf_diffusion_real_error_dataset_id, extensible_dataspace_id, hdf5_status )
  call H5SSELECT_HYPERSLAB_F( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call H5DWRITE_F( diagnostics_linf_diffusion_real_error_dataset_id, H5T_NATIVE_DOUBLE, real( linf_diffusion_error, kind=dp ), &
                   count, hdf5_status, diagnostics_length1_vector_dataspace_id, extensible_dataspace_id )
  call H5DGET_SPACE_F( diagnostics_linf_diffusion_imag_error_dataset_id, extensible_dataspace_id, hdf5_status )
  call H5SSELECT_HYPERSLAB_F( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call H5DWRITE_F( diagnostics_linf_diffusion_imag_error_dataset_id, H5T_NATIVE_DOUBLE, aimag( linf_diffusion_error ), &
                   count, hdf5_status, diagnostics_length1_vector_dataspace_id, extensible_dataspace_id )

  ! Timestep wall time.
  call H5DGET_SPACE_F( diagnostics_wall_time_step_dataset_id, extensible_dataspace_id, hdf5_status )
  call H5SSELECT_HYPERSLAB_F( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call H5DWRITE_F( diagnostics_wall_time_step_dataset_id, H5T_NATIVE_DOUBLE, elapsed_time, &
                   count, hdf5_status, diagnostics_length1_vector_dataspace_id, extensible_dataspace_id )

  ! CFL information in x
  call H5DGET_SPACE_F( diagnostics_cfl_x_dataset_id, extensible_dataspace_id, hdf5_status )
  call H5SSELECT_HYPERSLAB_F( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call H5DWRITE_F( diagnostics_cfl_x_dataset_id, H5T_NATIVE_DOUBLE, cfl_max_x, &
                   count, hdf5_status, diagnostics_length1_vector_dataspace_id, extensible_dataspace_id )

  ! CFL information in y.
  call H5DGET_SPACE_F( diagnostics_cfl_y_dataset_id, extensible_dataspace_id, hdf5_status )
  call H5SSELECT_HYPERSLAB_F( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call H5DWRITE_F( diagnostics_cfl_y_dataset_id, H5T_NATIVE_DOUBLE, cfl_max_y, &
                   count, hdf5_status, diagnostics_length1_vector_dataspace_id, extensible_dataspace_id )

  ! CFL information in z.
  call H5DGET_SPACE_F( diagnostics_cfl_z_dataset_id, extensible_dataspace_id, hdf5_status )
  call H5SSELECT_HYPERSLAB_F( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call H5DWRITE_F( diagnostics_cfl_z_dataset_id, H5T_NATIVE_DOUBLE, cfl_max_z, &
                   count, hdf5_status, diagnostics_length1_vector_dataspace_id, extensible_dataspace_id )

  ! Timestep size
  call H5DGET_SPACE_F( diagnostics_time_step_size_dataset_id, extensible_dataspace_id, hdf5_status )
  call H5SSELECT_HYPERSLAB_F( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call H5DWRITE_F( diagnostics_time_step_size_dataset_id, H5T_NATIVE_DOUBLE, dt, &
                   count, hdf5_status, diagnostics_length1_vector_dataspace_id, extensible_dataspace_id )

  ! Simulation Time
  call H5DGET_SPACE_F( diagnostics_simulation_time_dataset_id, extensible_dataspace_id, hdf5_status )
  call H5SSELECT_HYPERSLAB_F( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call H5DWRITE_F( diagnostics_simulation_time_dataset_id, H5T_NATIVE_DOUBLE, simulation_time, &
                   count, hdf5_status, diagnostics_length1_vector_dataspace_id, extensible_dataspace_id )

  ! Periodically, flush the updates to disk.  We only do this when we're going
  ! to write field variables so we don't unnecessarily introduce overhead
  ! without much gain.  While individual step timings/error are useful, they
  ! don't have to be available immediately.
  if (mod( timestep_index, timesteps_between_writes ) == 0) then
     call H5FFLUSH_F( diagnostics_gmres_diffusion_iterations_real_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
     call H5FFLUSH_F( diagnostics_gmres_diffusion_iterations_imag_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
     call H5FFLUSH_F( diagnostics_gmres_poisson_iterations_real_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
     call H5FFLUSH_F( diagnostics_gmres_poisson_iterations_imag_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
     call H5FFLUSH_F( diagnostics_gmres_viscous_iterations_real_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
     call H5FFLUSH_F( diagnostics_gmres_viscous_iterations_imag_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
     call H5FFLUSH_F( diagnostics_l2_poisson_error_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
     call H5FFLUSH_F( diagnostics_l2_poisson_schur_real_error_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
     call H5FFLUSH_F( diagnostics_l2_poisson_schur_imag_error_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
     call H5FFLUSH_F( diagnostics_linf_divergence_error_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
     call H5FFLUSH_F( diagnostics_linf_viscous_real_error_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
     call H5FFLUSH_F( diagnostics_linf_viscous_imag_error_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
     call H5FFLUSH_F( diagnostics_cfl_x_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
     call H5FFLUSH_F( diagnostics_cfl_y_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
     call H5FFLUSH_F( diagnostics_cfl_z_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
     call H5FFLUSH_F( diagnostics_time_step_size_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
     call H5FFLUSH_F( diagnostics_simulation_time_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
     call H5FFLUSH_F( diagnostics_linf_diffusion_imag_error_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
     call H5FFLUSH_F( diagnostics_linf_diffusion_real_error_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
     call H5FFLUSH_F( diagnostics_wall_time_step_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
  end if

end subroutine write_diagnostics_timestep

subroutine write_diagnostics_execution_stats( elapsed_cpu_time, elapsed_field_io_time, &
                                        elapsed_restart_io_time, &
                                        elapsed_diffusion_time, &
                                        elapsed_poisson_time, &
                                        elapsed_viscous_time, &
                                        elapsed_total_time, &
                                        total_norm_numeric_error )
! Writes out the solver's execution statistics.
!
! XXX: pull the gather of CPU time into this routine when parallel writing
!      is implemented.

  use constants, only:      nprocs
  use HDF5, only:           H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, &
                            hid_t, hsize_t, &
                            h5dclose_f, h5dcreate_f, h5dwrite_f, &
                            h5sclose_f, h5screate_simple_f
  use io_diagnostics, only: diagnostics_file_id, &
                            diagnostics_scalar_dataspace_id, &
                            diagnostics_name_cpu_time, diagnostics_name_total_numeric_error, &
                            diagnostics_name_null_basis_wall_time, diagnostics_name_setup_wall_time, &
                            diagnostics_name_field_io_time, diagnostics_name_restart_io_time, &
                            diagnostics_name_diffusion_solve_time, &
                            diagnostics_name_poisson_solve_time, &
                            diagnostics_name_viscous_solve_time, &
                            diagnostics_name_total_wall_time
  use precision, only:      dp

  implicit none

  real(kind=dp), intent(in)      :: elapsed_cpu_time
  real(kind=dp), intent(in)      :: elapsed_field_io_time
  real(kind=dp), intent(in)      :: elapsed_restart_io_time
  real(kind=dp), intent(in)      :: elapsed_diffusion_time
  real(kind=dp), intent(in)      :: elapsed_poisson_time
  real(kind=dp), intent(in)      :: elapsed_viscous_time
  real(kind=dp), intent(in)      :: elapsed_total_time
  real(kind=dp), intent(in)      :: total_norm_numeric_error

  integer                        :: data_rank
  integer(hsize_t), dimension(2) :: data_dimensions

  ! HDF5 identifiers and status code.
  integer(hid_t)                 :: dataset_id
  integer(hid_t)                 :: vector_dataspace_id
  integer                        :: hdf5_status

  ! Create a scalar for the numeric error.
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_total_numeric_error, H5T_IEEE_F64LE, &
                    diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_DOUBLE, total_norm_numeric_error, data_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )

  ! Create a 1D dataspace to hold timings from each of the MPI ranks.
  data_rank          = 1
  data_dimensions(1) = nprocs
  data_dimensions(2) = 1
  call H5SCREATE_SIMPLE_F( data_rank, data_dimensions, vector_dataspace_id, hdf5_status )

  ! Create a vector of scalars for the elapsed CPU time for each of the MPI
  ! ranks, log_cpu_time.
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_cpu_time, H5T_IEEE_F64LE, &
                    vector_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_DOUBLE, elapsed_cpu_time, data_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )

  ! Close the 1D dataspace.
  call H5SCLOSE_F( vector_dataspace_id, hdf5_status )

  ! Create a scalar for the elapsed diagnostics I/O time.
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_field_io_time, H5T_IEEE_F64LE, &
                    diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_DOUBLE, elapsed_field_io_time, data_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )

  ! Create a scalar for the elapsed restart I/O time.
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_restart_io_time, H5T_IEEE_F64LE, &
                    diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_DOUBLE, elapsed_restart_io_time, data_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )

  ! Create a scalar for the diffusion solve time.
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_diffusion_solve_time, H5T_IEEE_F64LE, &
                    diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_DOUBLE, elapsed_diffusion_time, data_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )

  ! Create a scalar for the Poisson solve time.
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_poisson_solve_time, H5T_IEEE_F64LE, &
                    diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_DOUBLE, elapsed_poisson_time, data_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )

  ! Create a scalar for the viscous solve time.
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_viscous_solve_time, H5T_IEEE_F64LE, &
                    diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_DOUBLE, elapsed_viscous_time, data_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )

  ! Create a scalar for the total wall time elapsed.
  call H5DCREATE_F( diagnostics_file_id, diagnostics_name_total_wall_time, H5T_IEEE_F64LE, &
                    diagnostics_scalar_dataspace_id, dataset_id, hdf5_status )
  call H5DWRITE_F( dataset_id, &
                   H5T_NATIVE_DOUBLE, elapsed_total_time, data_dimensions, &
                   hdf5_status )
  call H5DCLOSE_F( dataset_id, hdf5_status )

end subroutine write_diagnostics_execution_stats

subroutine close_diagnostics_file
! Closes the diagnostics file.  Once closed, no further I/O may be performed on the
! diagnostics file without opening it again.  The timestep dataset identifiers
! (diagnostics_number_steps_dataset_id, diagnostics_wall_time_step_dataset_id,
! and diagnostics_step_numeric_error_dataset_id), the various dataspace identifiers
! (diagnostics_grid_dataspace_id, diagnostics_length1_vector_dataspace_id, and
! diagnostics_scalar_dataspace_id), and the file identifier (diagnostics_file_id) are
! closed.

  use constants, only:       rank, root
  use HDF5, only:            hid_t, &
                             h5dclose_f, &
                             h5fclose_f, &
                             h5sclose_f
  use io_diagnostics, only:  diagnostics_file_id, &
                             diagnostics_grid_2D_dataspace_id, &
                             diagnostics_length1_vector_dataspace_id, &
                             diagnostics_scalar_dataspace_id, &
                             diagnostics_gmres_diffusion_iterations_real_dataset_id, &
                             diagnostics_gmres_diffusion_iterations_imag_dataset_id, &
                             diagnostics_gmres_poisson_iterations_real_dataset_id, &
                             diagnostics_gmres_poisson_iterations_imag_dataset_id, &
                             diagnostics_gmres_viscous_iterations_real_dataset_id, &
                             diagnostics_gmres_viscous_iterations_imag_dataset_id, &
                             diagnostics_l2_poisson_error_dataset_id, &
                             diagnostics_l2_poisson_schur_real_error_dataset_id, &
                             diagnostics_l2_poisson_schur_imag_error_dataset_id, &
                             diagnostics_linf_divergence_error_dataset_id, &
                             diagnostics_linf_viscous_real_error_dataset_id, &
                             diagnostics_linf_viscous_imag_error_dataset_id, &
                             diagnostics_linf_diffusion_real_error_dataset_id, &
                             diagnostics_linf_diffusion_imag_error_dataset_id, &
                             diagnostics_wall_time_step_dataset_id, &
                             diagnostics_cfl_x_dataset_id, &
                             diagnostics_cfl_y_dataset_id, &
                             diagnostics_cfl_z_dataset_id, &
                             diagnostics_time_step_size_dataset_id, &
                             diagnostics_simulation_time_dataset_id

  implicit none

  integer :: hdf5_status

  if (rank == root) then

     ! Close the datasets used to maintain the timing, errors, etc...
     call H5DCLOSE_F( diagnostics_gmres_diffusion_iterations_real_dataset_id, hdf5_status )
     call H5DCLOSE_F( diagnostics_gmres_diffusion_iterations_imag_dataset_id, hdf5_status )
     call H5DCLOSE_F( diagnostics_gmres_poisson_iterations_real_dataset_id, hdf5_status )
     call H5DCLOSE_F( diagnostics_gmres_poisson_iterations_imag_dataset_id, hdf5_status )
     call H5DCLOSE_F( diagnostics_gmres_viscous_iterations_real_dataset_id, hdf5_status )
     call H5DCLOSE_F( diagnostics_gmres_viscous_iterations_imag_dataset_id, hdf5_status )
     call H5DCLOSE_F( diagnostics_l2_poisson_error_dataset_id, hdf5_status )
     call H5DCLOSE_F( diagnostics_l2_poisson_schur_real_error_dataset_id, hdf5_status )
     call H5DCLOSE_F( diagnostics_l2_poisson_schur_imag_error_dataset_id, hdf5_status )
     call H5DCLOSE_F( diagnostics_linf_divergence_error_dataset_id, hdf5_status )
     call H5DCLOSE_F( diagnostics_linf_viscous_real_error_dataset_id, hdf5_status )
     call H5DCLOSE_F( diagnostics_linf_viscous_imag_error_dataset_id, hdf5_status )
     call H5DCLOSE_F( diagnostics_linf_diffusion_real_error_dataset_id, hdf5_status )
     call H5DCLOSE_F( diagnostics_linf_diffusion_imag_error_dataset_id, hdf5_status )
     call H5DCLOSE_F( diagnostics_wall_time_step_dataset_id, hdf5_status )
     call H5DCLOSE_F( diagnostics_cfl_x_dataset_id, hdf5_status )
     call H5DCLOSE_F( diagnostics_cfl_y_dataset_id, hdf5_status )
     call H5DCLOSE_F( diagnostics_cfl_z_dataset_id, hdf5_status )
     call H5DCLOSE_F( diagnostics_time_step_size_dataset_id, hdf5_status )
     call H5DCLOSE_F( diagnostics_simulation_time_dataset_id, hdf5_status )

     ! Close the grid, scalar, and vector scalar dataspaces.
     call H5SCLOSE_F( diagnostics_grid_2D_dataspace_id, hdf5_status )
     call H5SCLOSE_F( diagnostics_length1_vector_dataspace_id, hdf5_status )
     call H5SCLOSE_F( diagnostics_scalar_dataspace_id, hdf5_status )

     ! Close the file so that all of the pending writes are flushed out.
     call H5FCLOSE_F( diagnostics_file_id, hdf5_status )

  endif

end subroutine close_diagnostics_file
