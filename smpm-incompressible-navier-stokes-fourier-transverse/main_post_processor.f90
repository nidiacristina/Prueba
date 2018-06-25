program smpm_post_processor

  ! Post Processor for SMPM. Flow properties computed include:
  ! - derivatives for density and velocity field
  !   -> curl
  !   -> divergence
  ! - enstrophy integral kernel
  ! - kinetic energy integral kernel
  ! - instantaneous pressure
  ! - streamfunction
  ! - pressure
  !
  ! Note: There streamfunction has been commented out but not
  !       removed. We are still not sure how useful it is to
  !       keep it since not all runs will require such an
  !       expensive calculation.
  !
  ! 9 August 2018
  ! Gustavo Rivera

  ! Declare modules.
  use, intrinsic :: iso_c_binding
  use constants, only:               g, nprocs, nsg, nsuby, rank, root, rho_0
  use field_variables, only:         Cux0, Cuy0, Cuz0, divergence, &
                                     drhodx, drhody, drhodz, dubcdz, &
                                     dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, &
                                     enstrophy, kinetic_energy, Nux0, Nuy0, Nuz0, omega_x, &
                                     omega_y, omega_z, p, rho, streamfunction, ubc, ux, uy, uz
  use HDF5, only:                    h5close_f, h5open_f
  use legendre, only:                filter_xz
  use mpi, only:                     MPI_COMM_WORLD
  use options, only:                 fname_init, fname_setup, &
                                     read_from_setupfile, use_gravity_force, write_to_setupfile
  use parallel_linear_algebra, only: pdot_product, pmaxval, pnorm2
  use postprocessor, only:           compute_streamfunction, compute_pressure
  use precision, only:               dp
  use timestepping, only:            dt, logflag
  use woodbury_matrices, only:       rpk

  implicit none
  include 'fftw3.f03'

  ! Allocatable arrays.
  real(kind=dp), allocatable, dimension(:,:)        :: stream_rhs

  ! General use variables.
  real(kind=dp)                                     :: pi

  ! MPI and HDF5 status variables
  integer                                           :: ierr
  integer                                           :: hdf5_status

  ! Variables for GMRES computation associated with streamfunction
  character(len=64)                                 :: caststr
  integer                                           :: gmres_iterations, gmres_restart
  real(kind=dp)                                     :: gmres_tol, gmres_tol_viscous

  ! Declare variables for reading command line inputs.
  integer                                           :: iargc
  integer                                           :: n_arg
  character(len=100)                                :: infile_name
  character(len=100)                                :: outfile_name
  character(len=5)                                  :: procstr
  logical                                           :: io_flag

  ! Declare the functions used with GMRES to solve the viscous problems.
  external                                          :: apply_streamfunction_matrix

  ! Declare variables associated with the pressure solver
  real(kind=dp)                                     :: elapsed_time_null_basis

  ! Initialize MPI.
  call MPI_INIT( ierr )

  ! Initialize the HDF5 library.
  call H5OPEN_F( hdf5_status )

  ! Get my rank and total number of processes.
  call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

  ! Get pi.
  pi = 2.0_dp * acos( 0.0_dp )

  ! Write the process number to a string for display.
  write( procstr, '(I3)' ) rank
  procstr = trim(procstr)

  ! Read command line inputs for input file.
  n_arg = iargc()
  if (n_arg < 2) then

     if ( rank == root ) then
        write(*,*) 'Error: too few input arguments:', n_arg
        write(*,'(a)') '  Usage: smpm_post_processor runname_in fieldname_out'
     endif

     call MPI_FINALIZE( ierr )

     stop
  else

     ! Get the input and output file names.
     call GETARG( 1, infile_name )
     call GETARG( 2, outfile_name )

  endif

  ! Present welcome message.
  call notify( '==============================================================================================' )
  call notify( ' ' )
  call notify( '  SMPM post-processing tool                                                                   ' )
  call notify( ' ' )
  call notify( '==============================================================================================' )
  call notify( ' ' )

  ! Read the input file.
  ! Initialize our configuration from the input file.  This reads parameters,
  ! defaults optional values, and initializes parameters derived from the
  ! inputs.
  call notify( 'Reading configuration from ' // infile_name )
  call initialize_configuration( infile_name, io_flag )
  call smpm_assert( io_flag, 'Failed to read a valid configuration.' )

  ! Allocate space for our internal state.
  call allocate_solver_variables()

  ! Display input information.
  call display_inputs

  ! Set the initial conditions.
  call notify( 'Reading the initialization file.' )
  call read_initfile_data( fname_init, io_flag )
  call smpm_assert( io_flag, 'Failed to open ' // fname_init )


  ! Compute the deformation maps.
  call notify( 'Computing the deformation maps.' )
  call setup_deformation_maps

  ! Setup the transverse grid.
  call notify( 'Setting up the transverse grid.' )
  call setup_transverse()
  call setup_filter_y()

  ! Setup the Poisson capacitance and deflation operators, along with the null
  ! vectors.  These will be read from disk if they were created in a previous
  ! run, otherwise computed.
  if ( compute_pressure ) then

     call setup_poisson_operators( read_from_setupfile, write_to_setupfile, &
                                   fname_setup, io_flag, elapsed_time_null_basis )
     call smpm_assert( io_flag, 'Failed to read or write ' // fname_setup )
  else

     call notify( 'No pressure computation requested. Skipping Poisson operator setup.' )

  endif

  ! Compute the divergence of the flow.
  allocate( divergence( 1:rpk, 1:nsuby ) )
  call notify( 'Computing the divergence.' )
  divergence = 0.0_dp
  call compute_3D_divergence( ux, uy, uz, divergence )

  ! Compute the field derivatives
  call notify( 'Computing the field derivatives.')
  dudx = 0.0_dp
  dudy = 0.0_dp
  dudz = 0.0_dp
  call compute_3D_gradient( ux, dudx, dudy, dudz)
  dvdx = 0.0_dp
  dvdy = 0.0_dp
  dvdz = 0.0_dp
  call compute_3D_gradient( uy, dvdx, dvdy, dvdz)
  dwdx = 0.0_dp
  dwdy = 0.0_dp
  dwdz = 0.0_dp
  call compute_3D_gradient( uz, dwdx, dwdy, dwdz)
  drhodx = 0.0_dp
  drhody = 0.0_dp
  drhodz = 0.0_dp
  call compute_3D_gradient( rho, drhodx, drhody, drhodz)

  ! Compute the vorticity.
  allocate( omega_x(1:rpk ,1:nsuby) )
  allocate( omega_y(1:rpk ,1:nsuby) )
  allocate( omega_z(1:rpk ,1:nsuby) )
  call notify( 'Computing the vorticity.' )
  omega_x = 0.0_dp
  omega_y = 0.0_dp
  omega_z = 0.0_dp
  call compute_3D_curl( ux, uy, uz, omega_x, omega_y, omega_z )

  ! Compute the enstrophy integral kernel
  allocate( enstrophy(1:rpk,1:nsuby) )
  enstrophy = 0.0_dp
  call notify( 'Computing the enstrophy.' )
  enstrophy = omega_x * omega_x + omega_y * omega_y + omega_z * omega_z

  ! Compute the streamfunction of the flow.
  allocate( streamfunction(1:rpk,1:nsuby), stream_rhs(1:rpk,1:nsuby) )
  streamfunction = 0.0_dp
  stream_rhs     = 0.0_dp
  if ( compute_streamfunction ) then
     if (nsuby > 1 ) then
        call notify( '3D run... Skipping streamfunction.' )
     else
        call notify( 'Computing the streamfunction.' )
        call compute_divergence( ux, uz, streamfunction )
        gmres_tol                  = gmres_tol_viscous
        gmres_iterations           = nsg
        gmres_restart              = nsg
        streamfunction             = 0.0_dp
        stream_rhs                 = -1.0_dp * omega_y
        call compute_gmres_householder( streamfunction, stream_rhs, nsuby * rpk, gmres_tol, gmres_iterations, &
                                        gmres_restart, apply_streamfunction_matrix )
        write( caststr, '(D17.10)' ) gmres_tol
        call notify('   Streamfunction computation converged with residual: ' // adjustl( caststr ) )
        write( caststr, '(I10)' ) gmres_iterations
        call notify('   Streamfunction computation converged on iteration : ' // adjustl( caststr ) )
     endif
  endif

  ! Compute the Kinetic Energy Integral Kernel
  allocate( kinetic_energy(1:rpk,1:nsuby) )
  kinetic_energy = 0.0_dp
  call notify( 'Computing the kinetic energy.' )
  kinetic_energy = ux * ux + uy * uy + uz * uz

  ! Compute the instantaneous pressure
  if ( compute_pressure ) then
     call notify( 'Computing the pressure. ' )

     ! Step 1: Compute the Nonlinear terms

     ! Step 1.1: Treat the horizontal (x) momentum equation

     ! Step 1.1.1: Compute nonlinear terms
     call apply_smpm_advection( Nux0, ux, ux, uy, uz)

     ! Step 1.1.2: Add contribution of background current
     Nux0 = Nux0 - ( ubc * dudx + uz + dubcdz)

     ! Step 1.1.3: Filter the horizontal momentum
     call apply_filter_xz( Nux0, filter_xz )
     call apply_filter_y(  Nux0, Nux0 )

     ! Step 1.2: Treat the transverse (y) momentum equation
     ! Step 1.2.1: Compute nonlinear terms
     call apply_smpm_advection( Nuy0, uy, ux, uy, uz )

     ! Step 1.2.2: Add contribution of background current
     Nuy0 = Nuy0 - (ubc * dvdx)

     ! Step 1.2.3: Filter the transverse momentum
     call apply_filter_xz( Nuy0, filter_xz )
     call apply_filter_y(  Nuy0, Nuy0 )

     ! Step 1.3: Treat the vertical (z) momentum equation
     ! Step 1.3.1: Compute nonlinear terms
     call apply_smpm_advection( Nuz0, uz, ux, uz, uz )

     ! Step 1.3.1.1: Include gravity force
     if ( use_gravity_force .eqv. .true.) then
        Nuz0 = Nuz0 - g * rho / rho_0
     endif

     ! Step 1.3.2: Add contribution of background current
     Nuz0 = Nuz0 - ubc * dwdx

     ! Step 1.3.3: Filter the vertical momentum equation
     call apply_filter_xz( Nuz0, filter_xz )
     call apply_filter_y(  Nuz0, Nuz0 )

     ! Step 2: Setup the RHS of the Pressure Poisson Equation
     ! Step 2.1: Compute divergence of velocity field
     call compute_3D_divergence( Nux0, Nuy0, Nuz0, p)

     ! Step 2.1.1: Compute the double curl terms
     call compute_double_curl( ux, uy, uz, Cux0, Cuy0, Cuz0 )

     ! Step 2.1.2: Reset the timestepping coefficients
     call set_timestepping_coefficients( 1, dt, 0.0_dp, 0.0_dp )

     ! Step 3: Solve the Pressure Poisson Equation
     logflag = .true. ! Activate log for Pressure-GMRES information
     call solve_3D_pressure( p )

  endif ! End of pressure

  ! Open the post processor output file.
  if ( rank == root ) then
     call notify( 'Writing to post file: ' // trim( outfile_name ) // '.h5' )
     call open_post_file( trim( outfile_name ) // '.h5' )
  endif

  ! Write the header and the data (this subroutine ought to be renamed).
  call notify( 'Writing grid to post file. ' )
  call write_post_grid( )

  call notify( 'Writing data to post file. ' )
  call write_post_data

  ! Close the post file.
  if ( rank == root ) call close_post_file( )

  ! Create XDMF File describing the file
  call write_xdmf_postprocessor_file( trim( outfile_name) // '.xmf', &
                                      trim( outfile_name) // '.h5' )

  ! shutdown the HDF5 library before exiting.
  call H5CLOSE_F( hdf5_status )

  call notify( 'Postprocessing Complete. ' )

  ! Finalize MPI.
  call MPI_FINALIZE( ierr )

end program smpm_post_processor
