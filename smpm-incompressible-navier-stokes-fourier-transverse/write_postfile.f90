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
!       the /post/ group appeared to not contain any sub-groups, nor the
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
!       Currently, the post file is written such that each timestep written
!       lives within its own group.

subroutine open_post_file( post_file_name )
! Opens the post file for writing throughout the solver's life time.  If
! the specified post file already exists, it will be overwritten.  Nothing
! is written to the post file.
!
! This routine opens both the scalar and grid data spaces.

  use constants, only: n, nsubx, nsuby, nsubz
  use HDF5, only:      H5F_ACC_TRUNC_F, H5S_SCALAR_F, &
                       hid_t, hsize_t, &
                       h5fcreate_f, &
                       h5screate_f, h5screate_simple_f
  use io_post, only:   post_file_id, &
                       post_grid_dataspace_id, post_scalar_dataspace_id

  implicit none

  character(len=*), intent(in)   :: post_file_name

  integer                        :: data_rank
  integer(hsize_t), dimension(3) :: data_dimensions
  integer                        :: hdf5_status

  ! Create our output file, overwriting any existing one.
  call H5FCREATE_F( trim (post_file_name ), &
                    H5F_ACC_TRUNC_F, post_file_id, hdf5_status )

  ! Create a scalar dataspace that we'll use to create scalar datatypes
  ! with.  each of the non-grid variables are tied to this dataspace.
  call H5SCREATE_F( H5S_SCALAR_F, post_scalar_dataspace_id, hdf5_status )

  ! Create a 2D dataspace that we'll use for grid-interface variables
  ! (coordinates, velocities, density, etc).
  data_rank          = 3
  data_dimensions(1) = n * nsubz
  data_dimensions(2) = n * nsubx
  data_dimensions(3) = nsuby
  call H5SCREATE_SIMPLE_F( data_rank, data_dimensions, post_grid_dataspace_id, &
                           hdf5_status )

end subroutine open_post_file

subroutine write_post_grid
! Writes out the solver's grid to the post file's header.

  use constants, only:               n, nsg, nsubx, nsuby, nsubz, rank, root
  use geom, only:                    cx, cz
  use HDF5, only:                    H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER, H5T_STD_I32LE, &
                                     hid_t, hsize_t, &
                                     h5dclose_f, h5dcreate_f, h5dwrite_f, &
                                     h5gclose_f, h5gcreate_f, &
                                     h5screate_f, h5screate_simple_f
  use io_post, only:                 post_file_id, &
                                     post_grid_dataspace_id, post_scalar_dataspace_id, &
                                     post_name_grid, post_name_n, post_name_mx, &
                                     post_name_my, post_name_mz, post_name_x, post_name_y, post_name_z
  use precision, only:               dp
  use transverse, only:              cy

  implicit none

  real(kind=dp), allocatable, dimension(:,:)    :: full_cx, full_cy, full_cz

  integer(hsize_t), dimension(3)                :: data_dimensions

  ! HDF5 identifiers and status code.
  integer(hid_t)                                :: dataset_id
  integer(hid_t)                                :: group_id
  integer                                       :: hdf5_status

  ! Allocate space for a consolidated grid on the root.
  if (rank == root )then
     allocate( full_cx(1:nsg,1:nsuby), full_cy(1:nsg,1:nsuby), full_cz(1:nsg,1:nsuby) )
  else
     ! Avoid warnings about using uninitialized memory.  Non-root ranks wont
     ! touch the full grid during the collective operation.
     allocate( full_cx(1,1), full_cy(1,1), full_cz(1,1) )
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
     call H5GCREATE_F( post_file_id, post_name_grid, group_id, hdf5_status )
     call H5GCLOSE_F( group_id, hdf5_status )

     ! Create a scalar for the collocation dimension, n.
     call H5DCREATE_F( post_file_id, post_name_n, H5T_STD_I32LE, &
                       post_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_INTEGER, n, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar for the number of x sub-domains, nsubx.
     call H5DCREATE_F( post_file_id, post_name_mx, H5T_STD_I32LE, &
                       post_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_INTEGER, nsubx, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's x positions,
     ! cx4disk.
     call H5DCREATE_F( post_file_id, post_name_x, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, reshape( full_cx, (/ n * nsubz, n * nsubx, nsuby /) ), data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar for the number of z sub-domains, nsubz.
     call H5DCREATE_F( post_file_id, post_name_mz, H5T_STD_I32LE, &
                       post_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_INTEGER, nsubz, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's z positions,
     ! cz4disk.
     call H5DCREATE_F( post_file_id, post_name_z, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, reshape( full_cz, (/ n * nsubz, n * nsubx, nsuby /) ), data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar for the number of y grid points, nsuby.
     call H5DCREATE_F( post_file_id, post_name_my, H5T_STD_I32LE, &
                       post_scalar_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_INTEGER, int( nsuby, 4 ), data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's y positions,
     ! full_cy.
     call H5DCREATE_F( post_file_id, post_name_y, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, reshape( full_cy, (/ n * nsubz, n * nsubx, nsuby /) ), data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )



     deallocate( full_cx, full_cy, full_cz )
  end if

end subroutine write_post_grid

subroutine write_post_data
! Writes out the solver's initial conditions as a timestep #0 entry.

  use constants, only:               n, nsg, nsubx, nsuby, nsubz, rank, root
  use field_variables, only:         divergence, &
                                     drhodx, drhody, drhodz, &
                                     dudx, dudy, dudz, &
                                     dvdx, dvdy, dvdz, &
                                     dwdx, dwdy, dwdz, &
                                     enstrophy, kinetic_energy, &
                                     omega_x, omega_y, omega_z, &
                                     p, rho, rho_bar, &
                                     streamfunction, ux, uy, uz
  use HDF5, only:                    H5T_IEEE_F64LE, H5F_SCOPE_LOCAL_F, &
                                     H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER, H5T_STD_I32LE, &
                                     hid_t, hsize_t, &
                                     h5dclose_f, h5dcreate_f, h5dwrite_f, &
                                     h5fflush_f, &
                                     h5gclose_f, h5gcreate_f, &
                                     h5screate_f, h5screate_simple_f
  use io_post, only:                 post_file_id, &
                                     post_grid_dataspace_id, &
                                     post_name_div, &
                                     post_name_drhodx, post_name_drhody, post_name_drhodz, &
                                     post_name_dudx, post_name_dudy, post_name_dudz, &
                                     post_name_dvdx, post_name_dvdy, post_name_dvdz, &
                                     post_name_dwdx, post_name_dwdy, post_name_dwdz, &
                                     post_name_enstrophy, post_name_field_time, &
                                     post_name_kinetic_energy, &
                                     post_name_omega_x, post_name_omega_y, post_name_omega_z, &
                                     post_name_p, post_name_post, &
                                     post_name_rho, post_name_rho_bar, &
                                     post_name_stream, &
                                     post_name_ux, post_name_uy, post_name_uz
  use postprocessor,only:            field_time
  use precision, only:               dp

  implicit none

  ! Vectors for the full field arrays.
  real(kind=dp), allocatable, dimension(:,:) :: div4disk
  real(kind=dp), allocatable, dimension(:,:) :: drhodx4disk, drhody4disk, drhodz4disk
  real(kind=dp), allocatable, dimension(:,:) :: dudx4disk,   dudy4disk,   dudz4disk
  real(kind=dp), allocatable, dimension(:,:) :: dvdx4disk,   dvdy4disk,   dvdz4disk
  real(kind=dp), allocatable, dimension(:,:) :: dwdx4disk,   dwdy4disk,   dwdz4disk
  real(kind=dp), allocatable, dimension(:,:) :: enstrophy4disk, kineticenergy4disk
  real(kind=dp), allocatable, dimension(:,:) :: omegax4disk, omegay4disk, omegaz4disk, p4disk
  real(kind=dp), allocatable, dimension(:,:) :: rhobar4disk, ux4disk, uy4disk, uz4disk, rho4disk
  real(kind=dp), allocatable, dimension(:,:) :: stream4disk

  ! Dummy dimension variable for writing datasets.  All of the configuration
  ! values written are scalars, though a vector of dimensions must be supplied
  ! to each h5dwrite_f() call.  Note that this value is ignored by the HDF5
  ! library.
  integer(hsize_t), dimension(3)             :: data_dimensions

  ! HDF5 identifiers and status code.
  integer(hid_t)                             :: dataset_id
  integer(hid_t)                             :: group_id
  integer                                    :: hdf5_status

  ! Allocate space for the post processing variables.
  if ( rank == root ) then
     allocate( div4disk(1:nsg,1:nsuby) )
     allocate( drhodx4disk(1:nsg,1:nsuby), drhody4disk(1:nsg,1:nsuby), drhodz4disk(1:nsg,1:nsuby) )
     allocate( dudx4disk(1:nsg,1:nsuby), dudy4disk(1:nsg,1:nsuby), dudz4disk(1:nsg,1:nsuby) )
     allocate( dvdx4disk(1:nsg,1:nsuby), dvdy4disk(1:nsg,1:nsuby), dvdz4disk(1:nsg,1:nsuby) )
     allocate( dwdx4disk(1:nsg,1:nsuby), dwdy4disk(1:nsg,1:nsuby), dwdz4disk(1:nsg,1:nsuby) )
     allocate( kineticenergy4disk(1:nsg,1:nsuby), enstrophy4disk(1:nsg,1:nsuby) )
     allocate( omegax4disk(1:nsg,1:nsuby), omegay4disk(1:nsg,1:nsuby), omegaz4disk(1:nsg,1:nsuby) )
     allocate( p4disk(1:nsg,1:nsuby), rho4disk(1:nsg,1:nsuby), rhobar4disk(1:nsg,1:nsuby) )
     allocate( stream4disk(1:nsg,1:nsuby) )
     allocate( ux4disk(1:nsg,1:nsuby),uy4disk(1:nsg,1:nsuby), uz4disk(1:nsg,1:nsuby) )
  else
     allocate( div4disk(1:1,1:1) )
     allocate( drhodx4disk(1:1,1:1), drhody4disk(1:1,1:1), drhodz4disk(1:1,1:1) )
     allocate( dudx4disk(1:1,1:1), dudy4disk(1:1,1:1), dudz4disk(1:1,1:1) )
     allocate( dvdx4disk(1:1,1:1), dvdy4disk(1:1,1:1), dvdz4disk(1:1,1:1) )
     allocate( dwdx4disk(1:1,1:1), dwdy4disk(1:1,1:1), dwdz4disk(1:1,1:1) )
     allocate( kineticenergy4disk(1:1,1:1), enstrophy4disk(1:1,1:1) )
     allocate( omegax4disk(1:1,1:1), omegay4disk(1:1,1:1), omegaz4disk(1:1,1:1) )
     allocate( p4disk(1:1,1:1), rho4disk(1:1,1:1), rhobar4disk(1:1,1:1) )
     allocate( stream4disk(1:1,1:1) )
     allocate( ux4disk(1:1,1:1),uy4disk(1:1,1:1), uz4disk(1:1,1:1) )
  endif

  ! Gather all the data onto root.
  call gather_3D_array( divergence, div4disk )
  call gather_3D_array( drhodx, drhodx4disk )
  call gather_3D_array( drhody, drhody4disk )
  call gather_3D_array( drhodz, drhodz4disk )
  call gather_3D_array( dudx, dudx4disk )
  call gather_3D_array( dudy, dudy4disk )
  call gather_3D_array( dudz, dudz4disk )
  call gather_3D_array( dvdx, dvdx4disk )
  call gather_3D_array( dvdy, dvdy4disk )
  call gather_3D_array( dvdz, dvdz4disk )
  call gather_3D_array( dwdx, dwdx4disk )
  call gather_3D_array( dwdy, dwdy4disk )
  call gather_3D_array( dwdz, dwdz4disk )
  call gather_3D_array( enstrophy, enstrophy4disk )
  call gather_3D_array( kinetic_energy, kineticenergy4disk )
  call gather_3D_array( p, p4disk )
  call gather_3D_array( omega_x, omegax4disk )
  call gather_3D_array( omega_y, omegay4disk )
  call gather_3D_array( omega_z, omegaz4disk )
  call gather_3D_array( rho, rho4disk )
  call gather_3D_array( rho_bar, rhobar4disk )
  call gather_3D_array( ux, ux4disk )
  call gather_3D_array( uy, uy4disk )
  call gather_3D_array( uz, uz4disk )
  call gather_3D_array( streamfunction, stream4disk )


  ! Write out the post-processing solution to disk.
  if ( rank == root ) then

     ! Specify the dimensions of our grid variables, full_cx and full_cz.
     ! Note that these are ignored when specified with our scalar variables.
     data_dimensions(1) = n * nsubz
     data_dimensions(2) = n * nsubx
     data_dimensions(3) = nsuby

     ! Create the post group that holds all of the variables we're
     ! writing out.
     call H5GCREATE_F( post_file_id, post_name_post, group_id, hdf5_status )
     call H5GCLOSE_F( group_id, hdf5_status )

     ! Create a scalar grid variable for the grid's velocity divergence
     ! div4disk.
     call H5DCREATE_F( post_file_id, post_name_div, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, div4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's density derivatives
     ! drhodx4disk.
     call H5DCREATE_F( post_file_id, post_name_drhodx, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, drhodx4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! drhody4disk.
     call H5DCREATE_F( post_file_id, post_name_drhody, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, drhody4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! drhodz4disk.
     call H5DCREATE_F( post_file_id, post_name_drhodz, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, drhodz4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's u-velocity derivatives
     ! dudx4disk.
     call H5DCREATE_F( post_file_id, post_name_dudx, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, dudx4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! dudy4disk.
     call H5DCREATE_F( post_file_id, post_name_dudy, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, dudy4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! dudz4disk.
     call H5DCREATE_F( post_file_id, post_name_dudz, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, dudz4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's v-velocity derivatives
     ! dvdx4disk.
     call H5DCREATE_F( post_file_id, post_name_dvdx, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, dvdx4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! dvdy4disk.
     call H5DCREATE_F( post_file_id, post_name_dvdy, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, dvdy4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! dvdz4disk.
     call H5DCREATE_F( post_file_id, post_name_dvdz, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, dvdz4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's w-velocity derivatives
     ! dwdx4disk.
     call H5DCREATE_F( post_file_id, post_name_dwdx, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, dwdx4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! dwdy4disk.
     call H5DCREATE_F( post_file_id, post_name_dwdy, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, dwdy4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! dwdz4disk.
     call H5DCREATE_F( post_file_id, post_name_dwdz, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, dwdz4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the grid's enstrophy,
     ! enstrophy4disk
     call H5DCREATE_F( post_file_id, post_name_enstrophy, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, enstrophy4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar variable for the field time, field_time.
     call H5DCREATE_F( post_file_id, post_name_field_time, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, field_time, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for the kinetic energy,
     ! kineticenergy4disk
     call H5DCREATE_F( post_file_id, post_name_kinetic_energy, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, kineticenergy4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's vorticity values
     ! omegax4disk.
     call H5DCREATE_F( post_file_id, post_name_omega_x, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, omegax4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's vorticity values
     ! omegay4disk.
     call H5DCREATE_F( post_file_id, post_name_omega_y, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, omegay4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's vorticity values
     ! omegaz4disk.
     call H5DCREATE_F( post_file_id, post_name_omega_z, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, omegaz4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's pressure
     ! p4disk.
     call H5DCREATE_F( post_file_id, post_name_p, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, p4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's densities
     ! rho4disk.
     call H5DCREATE_F( post_file_id, post_name_rho, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, rho4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's densities
     ! rhobar4disk.
     call H5DCREATE_F( post_file_id, post_name_rho_bar, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, rhobar4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's streamfunction values
     ! stream4disk.
     call H5DCREATE_F( post_file_id, post_name_stream, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                       H5T_NATIVE_DOUBLE, stream4disk, data_dimensions, &
                       hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's x
     ! velocities, ux4disk.
     call H5DCREATE_F( post_file_id, post_name_ux, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, ux4disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's y
     ! velocities, uy4disk.
     call H5DCREATE_F( post_file_id, post_name_uy, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, uy4disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Create a scalar grid variable for each of the grid's z
     ! velocities, uz4disk.
     call H5DCREATE_F( post_file_id, post_name_uz, H5T_IEEE_F64LE, &
                       post_grid_dataspace_id, dataset_id, hdf5_status )
     call H5DWRITE_F( dataset_id, &
                      H5T_NATIVE_DOUBLE, uz4disk, data_dimensions, &
                      hdf5_status )
     call H5DCLOSE_F( dataset_id, hdf5_status )

     ! Flush the step we just constructed to disk.
     call H5FFLUSH_F( post_file_id, H5F_SCOPE_LOCAL_F, hdf5_status )

  endif

  deallocate( ux4disk, uy4disk, uz4disk, rho4disk, rhobar4disk, div4disk, stream4disk )
  deallocate( kineticenergy4disk, enstrophy4disk, omegax4disk, omegay4disk, omegaz4disk, p4disk)

end subroutine write_post_data

subroutine close_post_file
! Closes the post file.  Once closed, no further I/O may be performed on the
! post file without opening it again.  The timestep dataset identifiers
! (post_number_steps_dataset_id, post_wall_time_step_dataset_id,
! and post_step_numeric_error_dataset_id), the various dataspace identifiers
! (post_grid_dataspace_id, post_length1_vector_dataspace_id, and
! post_scalar_dataspace_id), and the file identifier (post_file_id) are
! closed.

  use HDF5, only:      hid_t, &
                       h5dclose_f, &
                       h5fclose_f, &
                       h5sclose_f
  use io_post, only:   post_file_id, &
                       post_grid_dataspace_id, &
                       post_scalar_dataspace_id

  implicit none

  integer :: hdf5_status

  ! Close the grid, scalar, and vector scalar dataspaces.
  call H5SCLOSE_F( post_grid_dataspace_id, hdf5_status )
  call H5SCLOSE_F( post_scalar_dataspace_id, hdf5_status )

  ! Close the file so that all of the pending writes are flushed out.
  call H5FCLOSE_F( post_file_id, hdf5_status )

end subroutine close_post_file

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
