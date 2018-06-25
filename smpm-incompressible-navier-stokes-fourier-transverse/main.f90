program smpm_incompressible_navier_stokes

  ! Declare modules.

  use, intrinsic :: iso_c_binding
  use constants, only:               g, nprocs, nsuby, rank, rho_0, root
  use errors, only:                  compute_numeric_error, &
                                     linf_divergence_error
  use field_variables, only:         Cux0, Cuy0, Cuz0, Cux1, Cuy1, Cuz1, &
                                     Cux2, Cuy2, Cuz2, div_u_int, &
                                     drhodx, drhody, drhodz, &
                                     dubcdz, dudx,dudy, dudz,  &
                                     dvdx, dvdy, dvdz, &
                                     dwdx, dwdy, dwdz, &
                                     Nrho0, Nrho1, Nrho2, &
                                     Nux0, Nux1, Nux2, Nuy0, Nuy1, Nuy2, Nuz0, Nuz1, Nuz2, p, px, py, pz, &
                                     rho, rho0, rho1, rho2, rho_bar_z, rho_bar_zz, rho_int, &
                                     ubc, &
                                     ux, ux0, ux1, ux2, ux_int, &
                                     uy, uy0, uy1, uy2, uy_int, &
                                     uz, uz0, uz1, uz2, uz_int
  use HDF5, only:                    h5close_f, h5open_f
  use legendre, only:                filter_xz
  use mpi, only:                     MPI_COMM_WORLD, MPI_DOUBLE_PRECISION
  use options, only:                 apply_restart, check_null_error, check_numerical_error, &
                                     do_interfacial_averaging, enforce_strong_velocity_bc, &
                                     exact_nullspace_projection, fname_setup, &
                                     read_from_setupfile, setup_and_stop, &
                                     timesteps_between_logs, timesteps_between_writes, &
                                     solve_momentum_equation, use_gravity_force, &
                                     use_scalar_transport, use_scalar_diffusion, &
                                     write_to_setupfile, apply_sponge_layer, timesteps_between_restarts
  use parallel_linear_algebra, only: pdot_product, pmaxval, pnorm2, znorm2
  use precision, only:               dp, int64
  use profiling, only:               time_cpu, time_setup, &
                                     time_setup_null_basis, &
                                     time_setup_null_error, &
                                     time_field_io, &
                                     time_restart_io, &
                                     time_timesteps, &
                                     time_compute_poisson, &
                                     time_compute_diffusion, &
                                     time_compute_viscous
  use sponge, only:                  raycoeff
  use stopwatch, only:               stopwatch_initialize, stopwatch_tick, &
                                     stopwatch_get_elapsed
  use timestepping, only:            a0, a1, a2, b0, b1, b2, dt, dt1, dt2, g0, logflag, simulation_time, &
                                     time_ndx, tend
  use transverse, only:              p_real_buff_FFT, p_complex_buff_FFT
  use woodbury_matrices, only:       rpk, uL

  implicit none
  include 'fftw3.f03'

  ! Allocatable arrays.
  real(kind=dp), allocatable, dimension(:, :)   :: Ax
  real(kind=dp), allocatable, dimension(:, :)   :: junk_in_x, junk_in_y

  ! Variable to store field index file.
  character(len=6)                              :: field_index_string

  ! Variable to store restart index file.
  character(len=6)                              :: restart_index_string

  ! Variables for simulation time advancement.
  integer                                       :: Nt

  ! Elapsed CPU seconds for the beginning and end of the solver's run to
  ! determine how much time was spent computing rather than waiting.
  real(kind=dp)                                 :: CPUtimeStart, CPUtimeEnd, CPUtimeElapsed

  ! Time references for the major components of the solver: solver start
  ! (beginning of setup), timesteps start (end of setup), and I/O start
  ! (per-timestep).  Also a smaller, section-specific reference that covers
  ! everything else.
  !
  ! NOTE: We explicitly request 8-byte integers to ensure the highest
  !       precision timer available is used by system_clock().
  integer(kind=int64)                           :: solver_start_time, &
                                                   timesteps_start_time, &
                                                   compute_start_time, &
                                                   io_start_time
  integer(kind=int64)                           :: clock_start

  ! Start of the solver in seconds since the epoch.
  real(kind=dp)                                 :: solver_start

  ! Time, in seconds, the current timestep took to execute.
  real(kind=dp)                                 :: elapsed_time_timestep

  ! General use variables.
  integer                                       :: ii
  real(kind=dp)                                 :: pi

  ! MPI and HDF5 status variables
  integer                                       :: ierr
  integer                                       :: hdf5_status

  ! Error variables.
  real(kind=dp)                                 :: thiserror

  ! Total norm'd numeric error.
  real(kind=dp)                                 :: total_norm_numeric_error

  ! Declare variables for reading command line inputs.
  integer                                       :: iargc
  integer                                       :: n_arg
  character(len=100)                            :: diagfile_name
  character(len=100)                            :: infile_name
  character(len=100)                            :: outfile_name
  character(len=100)                            :: restart_name
  character(len=64)                             :: caststr
  character(len=5)                              :: procstr
  logical                                       :: io_flag

  ! Initialize MPI.
  call MPI_INIT( ierr )

  ! Initialize the HDF5 library.
  call H5OPEN_F( hdf5_status )

  ! Initialize the timing module.
  call stopwatch_initialize()

  ! Get my rank and total number of processes.
  call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

  ! Record that our universe started at a specific time.
  call seconds_since_epoch( solver_start )

  ! Get pi.
  pi = 2.0_dp * acos( 0.0_dp )

  ! Set default execution options.
  check_null_error = .false.

  ! Commence timing.  Get references for elapsed CPU time and clock ticks
  ! so we can compute deltas throughout the setup and solve loop.
  call CPU_TIME( CPUtimeStart )
  call stopwatch_tick( solver_start_time )

  ! Read command line inputs for input file.
  n_arg = iargc()
  if (n_arg .ne. 4) then

     if ( rank == root ) then
        write(*,'(A,I1,A)') 'ERROR: Expected 3 input arguments though received ', n_arg, '.'
        write(*,*)          ''
        write(*,'(A)')      '  Usage: smpm_incompressible_navier_stokes <runname_in> <fieldname_out> <restart_file> <diagfile_name>'
        write(*,*)          ''
     endif

     call MPI_FINALIZE( ierr )

     stop
  else

     ! Get the input, output and restart file names.
     call GETARG( 1, infile_name )
     call GETARG( 2, outfile_name )
     call GETARG( 3, restart_name )
     call GETARG( 4, diagfile_name )

  endif

  ! Present welcome message.
  call notify( '==============================================================================================' )
  call notify( ' ' )
  call notify( '   Spectral Multidomain Penalty Method model for the Incompressible Navier Stokes Equations   ' )
  call notify( ' ' )
  call notify( '==============================================================================================' )
  call notify( ' ' )

  ! Initialize our configuration from the input file.  This reads parameters,
  ! defaults optional values, and initializes parameters derived from the
  ! inputs.
  call notify( 'Reading configuration from ' // infile_name )
  call initialize_configuration( infile_name, io_flag )
  call smpm_assert( io_flag, 'Failed to read a valid configuration.' )

  ! Allocate space for our internal state.
  call allocate_solver_variables()

  ! Initialize the sponge layer.
  raycoeff = 0.0_dp

  ! Initialize the simulation time.
  simulation_time = 0.0_dp

  ! Display input information.
  call display_inputs()

  ! Acquire our initial conditions.
  call setup_initial_conditions()

  ! Configure the deformation maps.
  call setup_deformation_maps()

  ! Setup the sponge layer if needed.
  if ( apply_sponge_layer ) then
     call notify( 'Generating Sponge Layer.' )
     call setup_sponge_layer()
  end if

  ! Setup the transverse grid.
  call notify( 'Setting up the transverse grid.' )
  call setup_transverse()
  call setup_filter_y()

  ! Compute the second derivative of the background density profile
  allocate( junk_in_x(1:rpk, 1:nsuby) )
  allocate( junk_in_y(1:rpk, 1:nsuby) )
  call compute_3D_gradient( rho_bar_z, junk_in_x, junk_in_y, rho_bar_zz )
  deallocate( junk_in_x )
  deallocate( junk_in_y )

  ! Check if we are solving the full NS Equations
  ! If we are, then we proceed to build the Poisson Operator. If not, we skip
  if ( solve_momentum_equation .eqv. .true. ) then
     ! Setup the Poisson capacitance and deflation operators, along with the null
     ! vectors.  These will be read from disk if they were created in a previous
     ! run, otherwise computed.
     call setup_poisson_operators( read_from_setupfile, write_to_setupfile, &
                                   fname_setup, io_flag, time_setup_null_basis )
     call smpm_assert( io_flag, 'Failed to read or write ' // fname_setup )

     ! Keep track of how long was spent in Woodbury setup.
     time_setup = stopwatch_get_elapsed( solver_start_time )

     ! If asked, compute the nullspace basis.
     if ( exact_nullspace_projection ) then
        call notify( 'Computing exact divergence null-space basis.' )
        call setup_nullspace_projection()
     endif

     ! If we want to just run setup and stop, stop the run.
     if ( setup_and_stop ) then
        call notify( 'Setup complete.  Halting execution.' )
        call MPI_FINALIZE( ierr )
        stop
     endif

     ! Check the error in the null space computation.
     call stopwatch_tick( clock_start )
     call compute_null_space_error()
     time_setup_null_error = stopwatch_get_elapsed( clock_start )
  else
        call notify( 'Skipping Poisson Operator Setup.' )
  endif

  ! Create our field file and write out the header.
  call notify( 'Writing to diagnostics file: ' // trim( diagfile_name ) // '.h5' )

  call stopwatch_tick( clock_start )
  call open_diagnostics_file( trim( diagfile_name ) // '.h5' , io_flag)
  call smpm_assert( io_flag, 'Failed to open the diagnostics file.' )

  call write_diagnostics_header( solver_start, time_setup_null_basis, &
                           time_setup_null_error, time_setup )
  time_field_io = stopwatch_get_elapsed( clock_start )

  ! Prepare for error computation.
  allocate( Ax(1:rpk, 1:nsuby) )
  total_norm_numeric_error = 0.0_dp

  ! Generate nonlinear and pressure bc's if simulation is a restart.
  ! Note: We emulate Step 1 of the main solver loop in generating the
  !       nonlinear terms. The user should remember that any changes
  !       to the computation of the nonlinear terms in the main solver
  !       loop needs to be reflected in the following section as well.
  if ( apply_restart ) then

     ! Notify the user, again, that this is a restart.
     call notify( ' ' )
     write( caststr, '(D17.10)' ) simulation_time
     call notify( 'Restarting Simulation at: ' // caststr )

     ! Compute some time-stepping parameters.
     Nt = ceiling( tend / dt )
     write( caststr, '(I8.0)' ) Nt
     call notify( 'Total time steps: ' // caststr )
     call notify( ' ' )

     ! Part 1: Advection/Transport Treatment of the Momentum and
     !         Scalar Equation

     ! Step 1.1: Advection treatment in u_x at time n-1 (ux1) and n-2 (ux2)

     ! Step 1.1.1: Advection solve in u_x
     call apply_smpm_advection( Nux1, ux1, ux1, uy1, uz1)
     call apply_smpm_advection( Nux2, ux2, ux2, uy2, uz2)

     ! Step 1.1.2: Add the background current's contribution
     call compute_3D_gradient( ux1, dudx, dudy, dudz )
     Nux1 = Nux1 - (ubc * dudx + uz1 * dubcdz)
     dudx = 0.0_dp
     dudy = 0.0_dp
     dudz = 0.0_dp
     call compute_3D_gradient( ux2, dudx, dudy, dudz )
     Nux2 = Nux2 - (ubc * dudx + uz2 * dubcdz)
     dudx = 0.0_dp
     dudy = 0.0_dp
     dudz = 0.0_dp

     ! Step 1.1.3: Apply the Sponge Layer
     Nux1 = Nux1 - (ux1 * raycoeff)
     Nux2 = Nux2 - (ux2 * raycoeff)

     ! Step 1.2: Advection treatment in u_y at time n-1 (uy1) and n-2 (uy2)

     ! Step 1.2.1: Advection solve in u_y
     call apply_smpm_advection( Nuy1, uy1, ux1, uy1, uz1)
     call apply_smpm_advection( Nuy2, uy2, ux2, uy2, uz2)

     ! Step 1.2.2: Add the background current's contribution
     call compute_3D_gradient( uy1, dvdx, dvdy, dvdz )
     Nuy1 = Nuy1 - (ubc * dvdx)
     dvdx = 0.0_dp
     dvdy = 0.0_dp
     dvdz = 0.0_dp
     call compute_3D_gradient( uy2, dvdx, dvdy, dvdz )
     Nuy2 = Nuy2 - (ubc * dvdx)
     dvdx = 0.0_dp
     dvdy = 0.0_dp
     dvdz = 0.0_dp

     ! Step 1.2.3: Apply the Sponge Layer
     Nuy1 = Nuy1 - (uy1 * raycoeff)
     Nuy2 = Nuy2 - (uy2 * raycoeff)

     ! Step 1.3: Advection treatment in u_z at time n-1 (uz1) and n-2 (uz2)

     ! Step 1.3.1: Advection solve in u_z
     call apply_smpm_advection( Nuz1, uz1, ux1, uy1, uz1)
     call apply_smpm_advection( Nuz2, uz2, ux2, uy2, uz2)

     ! Step 1.3.1.2: Include gravity force if requested
     if ( use_gravity_force .eqv. .true.) then
        Nuz1 = Nuz1 - g * rho1 / (rho_0)
        Nuz2 = Nuz2 - g * rho2 / (rho_0)
     endif

     ! Step 1.3.2: Add the background current's contribution.
     call compute_3D_gradient( uz1, dwdx, dwdy, dwdz )
     Nuz1 = Nuz1 - (ubc * dwdx)
     dwdx = 0.0_dp
     dwdy = 0.0_dp
     dwdz = 0.0_dp
     call compute_3D_gradient( uz2, dwdx, dwdy, dwdz )
     Nuz2 = Nuz2 - (ubc * dwdx)
     dwdx = 0.0_dp
     dwdy = 0.0_dp
     dwdz = 0.0_dp

     ! Step 1.3.3: Apply the Sponge Layer
     Nuz1 = Nuz1 - (uz1 * raycoeff)
     Nuz2 = Nuz2 - (uz2 * raycoeff)

     ! Step 1.4: Transport treatment of rho at time n-1 (uz1) and n-2 (uz2)

     ! Step 1.4.1: Nonlinear advection solve in rho + stratification terms
     call apply_smpm_transport( Nrho1, rho1, ux1, uy1, uz1)
     call apply_smpm_transport( Nrho2, rho2, ux2, uy2, uz2)

     ! Step 1.4.2: Add Contribution of background current and density and diffusion
     !             of background density.
     call compute_3D_gradient( rho1, drhodx, drhody, drhodz )
     ! XXX: these straficiation terms are multiplied by -1.
     !Nrho1 = Nrho1 - (uz1 * rho_bar_z + ubc * drhodx - nu_d * rho_bar_zz) ! W/ bg density diffusion
     ! XXX: this needs to be extended to 3D - GNT
     Nrho1 = Nrho1 - (uz1 * rho_bar_z + ubc * drhodx)   ! W/o bg density diffusion
     drhodx = 0.0_dp
     drhody = 0.0_dp
     drhodz = 0.0_dp
     call compute_3D_gradient( rho2, drhodx, drhody, drhodz )
     ! XXX: these straficiation terms are multiplied by -1.
     !Nrho2 = Nrho2 - (uz2 * rho_bar_z + ubc * drhodx - nu_d * rho_bar_zz) ! W/ bg density diffusion
     ! XXX: this needs to be extended to 3D - GNT
     Nrho2 = Nrho2 - (uz2 * rho_bar_z + ubc * drhodx)   ! W/o bg density diffusion
     drhodx = 0.0_dp
     drhody = 0.0_dp
     drhodz = 0.0_dp

     ! Step 1.4.3: Apply the sponge layer
     Nrho1 = Nrho1 - (rho1 * raycoeff)
     Nrho2 = Nrho2 - (rho2 * raycoeff)

     ! Part 4: Compute Pressure Boundary Conditions

     ! Step 4.1: Compute the double-curl for the pressure boundary conditions
     call compute_double_curl( ux0, uy0, uz0, Cux0, Cuy0, Cuz0 )
     call compute_double_curl( ux1, uy1, uz1, Cux1, Cuy1, Cuz1 )
     call compute_double_curl( ux2, uy2, uz2, Cux2, Cuy2, Cuz2 )

  else

     ! Compute some time-stepping parameters.
     Nt = ceiling( tend / dt )
     write( caststr, '(I8.0)' ) Nt
     call notify( 'Total time steps: ' // caststr )

  endif

  ! Track the beginning of all timesteps.
  call stopwatch_tick( timesteps_start_time )

  ! Commence time-stepping.
  call notify( ' ' )
  call notify( 'Commencing time-stepping.' )
  call notify( ' ' )

  ! Beginning of time loop.
  do time_ndx = 1, Nt

     ! Track the start of this timestep.
     call stopwatch_tick( clock_start )

     ! Check to see if logging is turned on for this time-step.
     if ( mod( time_ndx, timesteps_between_logs ) == 0 ) then
        logflag = .true.
     else
        logflag = .false.
     endif

     ! Solve this time-step.
     write( caststr, '(I8.0)' ) time_ndx
     call notify_cond( logflag, ' ' )
     call notify_cond( logflag, '   Solving time-step ' // caststr )
     call notify_cond( logflag, ' ' )

     ! Note the current simulation time.
     write( caststr, '(D17.10)' ) simulation_time
     call notify_cond( logflag, ' ' )
     call notify_cond( logflag, '   Simulation Time (s) ' // caststr )
     call notify_cond( logflag, ' ' )

     ! Step 0.0: If it's time, write out the field to disk.
     if ( mod( time_ndx, timesteps_between_writes ) == 0 ) then
        call notify_cond( logflag, ' ' )
        call notify_cond( logflag, '   Writing field file. ' )
        call notify_cond( logflag, ' ' )

        call stopwatch_tick( io_start_time )


        ! Create the field file by concatenating the underscore and write index
        write( field_index_string, '(I6.6)' ) time_ndx / timesteps_between_writes
        call open_field_file( trim( outfile_name ) // '_' // field_index_string // '.h5', io_flag )
        call smpm_assert( io_flag, 'Failed to open ' // trim( outfile_name ) // '_' // field_index_string // '.h5' )

        ! Write out the grid.
        call write_field_header_grid

        ! Write out the field variables rho, ux, uy, uz
        call write_field_variables( time_ndx )

        ! Close the field file.
        call close_field_file

        ! Create an XDMF file describing the file.
        call write_xdmf_field_file( trim( outfile_name ) // '_' // field_index_string //  '.xmf', &
                                      trim( outfile_name ) // '_' // field_index_string //  '.h5' )

        time_field_io = (time_field_io + &
                         stopwatch_get_elapsed( io_start_time ))
     endif

     ! Step 0.1: Check the CFL number and adjust time-step accordingly.
     call notify_cond( logflag, '      substep 0: checking CFL condition. ' )
     call check_cfl_number( dt, dt )

     ! Step 0.2: Compute our timestep coefficients.
     call stopwatch_tick( compute_start_time )
     call set_timestepping_coefficients( time_ndx, dt, dt1, dt2 )

     ! Part 1: Advection/Transport Treatment of Momentum and Scalar Equation

     ! Advection Treatment of Momentum Equation
     if ( solve_momentum_equation .eqv. .true. ) then

        call notify_cond( logflag, ' ' )
        call notify_cond( logflag, '      substep 1.1: applying advection. ')

        ! Step 1: Advection Treatment in u_x

        ! Step 1.1.1: Advection solve in u_x.
        call apply_smpm_advection( Nux0, ux0, ux0, uy0, uz0 )

        ! Step 1.1.2: Add the background current's contribution.
        call compute_3D_gradient( ux0, dudx, dudy, dudz )
        Nux0 = Nux0 - (ubc * dudx + uz0 * dubcdz)

        ! Step 1.1.3: Reset derivatives to zero.
        dudx = 0.0_dp
        dudy = 0.0_dp
        dudz = 0.0_dp

        ! Step 1.1.4: Apply the sponge layer.
        Nux0 = Nux0 - (ux0 * raycoeff)

        ! Step 1.1.5: Apply explicit AB3.
        ux_int   =       a0 * ux0  + a1 * ux1  + a2 * ux2 + &
                   dt * (b0 * Nux0 + b1 * Nux1 + b2 * Nux2)

        ! Step 1.1.6 Filter Velocity
        call apply_filter_xz( ux_int, filter_xz )
        call apply_filter_y( ux_int,ux_int )

        ! Step 1.2.0: Advection Treament in u_y
        if ( nsuby > 1 ) then

           ! Step 1.2.1: Nonlinear advection solve in u_y.
           call apply_smpm_advection( Nuy0, uy0, ux0, uy0, uz0 )

           ! Step 1.2.2: Add the background current's contribution.
           call compute_3D_gradient( uy0, dvdx, dvdy, dvdz )
           Nuy0 = Nuy0 - (ubc * dvdx)

           ! Step 1.2.3: Reset derivatives to zero.
           dvdx = 0.0_dp
           dvdy = 0.0_dp
           dvdz = 0.0_dp

           ! Step 1.2.4: Apply the sponge layer.
           Nuy0 = Nuy0 - (uy0 * raycoeff)

           ! Step 1.2.5: Apply explicit AB3.
           uy_int =       a0 * uy0  + a1 * uy1  + a2 * uy2 + &
                    dt * (b0 * Nuy0 + b1 * Nuy1 + b2 * Nuy2)

           ! Step 1.2.6: Filter Transverse Velocity
           call apply_filter_xz( uy_int, filter_xz )
           call apply_filter_y( uy_int,uy_int )

        else

           ! Initialize the y-component of the AB3 Scheme to 0
           uy_int = 0.0_dp

        endif

        ! Step 1.3.0: Advection Treatment in u_z/ with gravity

        ! Step 1.3.1: Nonlinear advection solve in u_z + gravity.
        call apply_smpm_advection( Nuz0, uz0, ux0, uy0, uz0 )

        ! Step 1.3.1.2: Include gravity force if requested
        if ( use_gravity_force .eqv. .true.) then
           Nuz0 = Nuz0 - g * rho0 / (rho_0)
        endif

        ! Step 1.3.2: Add the background current's contribution.
        call compute_3D_gradient( uz0, dwdx, dwdy, dwdz )
        Nuz0 = Nuz0 - (ubc * dwdx)

        ! Step 1.3.3 Reset derivatives to zero.
        dwdx = 0.0_dp
        dwdy = 0.0_dp
        dwdz = 0.0_dp

        ! Step 1.3.4: Apply the sponge layer.
        Nuz0 = Nuz0 - (uz0 * raycoeff)

        ! Step 1.3.5: Apply explicit AB3.
        uz_int =       a0 * uz0  + a1 * uz1  + a2 * uz2 + &
                 dt * (b0 * Nuz0 + b1 * Nuz1 + b2 * Nuz2)

        ! Step 1.3.6: Filter Vertical Velocity
        call apply_filter_xz( uz_int, filter_xz )
        call apply_filter_y( uz_int,uz_int )

     endif ! Advection Treatment of Momentum Equation

     ! Transport treatment for solve in rho + stratification terms
     if ( use_scalar_transport .eqv. .true. ) then

        call notify_cond( logflag, ' ' )
        call notify_cond( logflag, '      substep 1.1: applying transport. ')

        ! Step 1.4.1: Nonlinear advection solve in rho + stratification terms.
        call apply_smpm_transport( Nrho0, rho0, ux0, uy0, uz0 )

        ! Step 1.4.2: Add contribution of background current, background density
        !             and viscous diffusion of background density.
        call compute_3D_gradient( rho0, drhodx, drhody, drhodz )
        ! XXX: these stratifcation terms are multiplied by -1.
        !Nrho0 = Nrho0 - (uz0 * rho_bar_z + ubc * drhodx - nu_d * rho_bar_zz) ! W/ bg density diffusion

        ! XXX: this needs to be extended to 3D - GNT
        Nrho0 = Nrho0 - (uz0 * rho_bar_z + ubc * drhodx) ! W/o bg density diffusion

        ! Step 1.4.3: Reset derivatives to zero.
        drhodx = 0.0_dp
        drhody = 0.0_dp
        drhodz = 0.0_dp

        ! Step 1.4.4: Apply the sponge layer.
        Nrho0 = Nrho0 - (rho0 * raycoeff)

        ! Step 1.4.5: Apply explicit AB3.
        rho_int =       a0 * rho0  + a1 * rho1  + a2 * rho2 + &
                  dt * (b0 * Nrho0 + b1 * Nrho1 + b2 * Nrho2)

        ! Step 1.4.6: Filter perturbation density
        call apply_filter_xz( rho_int, filter_xz )
        call apply_filter_y( rho_int,rho_int )

     else ! Unlike the momentum equation, we include the explicit AB3 for the scalar
          ! equation because we are also interested in solving scalar advection-diffusion
          ! problems.

        ! Step 1.4.5: Apply explicit AB3
        rho_int =       a0 * rho0 + a1 * rho1 + a2 * rho2

     endif ! End transport treament for solve in rho +  stratification terms

     ! XXXXXXXX: what stopwatch needs to be ticked here?

     ! Part 2: Pressure Solve
     if ( solve_momentum_equation .eqv. .true. ) then

        call notify_cond( logflag, ' ' )
        call notify_cond( logflag, '      substep 2: solving Poisson problem. ')
        call notify_cond( logflag, ' ' )

        ! Step 2.0: Compute divergence of velocity field
        call compute_3D_divergence( ux_int, uy_int, uz_int, div_u_int )

        ! Step 2.1: Assemble the Right Hand Side
        p = - div_u_int / dt

        ! Step 2.2: set up the boundary conditions for the pressure solve.
        call solve_3D_pressure( p )

        ! Step 2.3: Compute the gradient of the computed pressure.
        call compute_3D_gradient( p, px, py, pz )

        ! Step 2.4: Filter the pressure gradients.
        call apply_filter_xz( px, filter_xz )
        call apply_filter_xz( py, filter_xz )
        call apply_filter_xz( pz, filter_xz )

        ! Step 2.5: Zero the gradient of pressure on the no-slip boundaries.
        ! XXX: This is wrong for deformed external boundaries.
        !      For deformed domains the proper logic is that if a boundary is
        !      free-slip then make <grad(p),n> = 0 on that boundary.  If the
        !      boundary is no-slip then make grap(p) = 0 on that boundary.
        if ( enforce_strong_velocity_bc ) then
           call notify_cond( logflag, ' ' )
           call notify_cond( logflag, '         Zeroing pressure gradient at boundary. ')
           call notify_cond( logflag, ' ' )

           call enforce_velocity_bc( px, py, pz )
        endif

        ! Step 2.6: Update the velocities with the pressure solution.
        ux_int = ux_int + dt * px
        uy_int = uy_int + dt * py
        uz_int = uz_int + dt * pz

        ! Step 2.7: Apply the exact null-space projection if requested.
        if ( exact_nullspace_projection .and. nsuby == 1 ) then
           call notify_cond( logflag, ' ' )
           call notify_cond( logflag, '         Applying Nullspace Projection to Velocity Field. ')
           call notify_cond( logflag, ' ' )

           call apply_nullspace_projection( ux_int, uz_int )
        endif

        ! Step 2.8: If asked, report divergence after velocity update.
        if ( check_numerical_error ) then

           ! Report div(u) prior to pressure update.
           div_u_int = 0.0_dp
           ux_int    = ux_int - dt * px
           uy_int    = uy_int - dt * py
           uz_int    = uz_int - dt * pz
           call compute_3D_divergence( ux_int, uy_int, uz_int, div_u_int )
           thiserror = pmaxval( reshape( abs( div_u_int ), (/ rpk * nsuby /) ) )
           write( caststr, '(D17.10)' ) thiserror
           call notify_cond( logflag, '         Maximum |div(u)| prior to update     : ' // caststr )

           ! Report the interaction between the original divergence and kernel
           ! vectors.
           div_u_int(:, 1) = uL * pdot_product( uL, div_u_int(:, 1) )
           thiserror = pmaxval( reshape( abs( div_u_int(:, 1) ), (/ rpk /) ) )
           write( caststr, '(D17.10)' ) thiserror
           call notify_cond( logflag, "               Maximum |u0*u0'*div(u)|        : " // caststr )

           ! Report div(u) after pressure update.
           ux_int = ux_int + dt * px
           uy_int = uy_int + dt * py
           uz_int = uz_int + dt * pz
           call compute_3D_divergence( ux_int, uy_int, uz_int, div_u_int )
           linf_divergence_error = pmaxval( reshape( abs( div_u_int ), (/ nsuby * rpk /) ) )
           write( caststr, '(D17.10)' ) linf_divergence_error
           call notify_cond( logflag, '         Maximum |div(u)| after update        : ' // caststr )

        endif ! End report of divergence after velocity update

        ! Step 2.9: Filter the intermediate velocities.
        call apply_filter_xz( ux_int, filter_xz )
        call apply_filter_xz( uy_int, filter_xz )
        call apply_filter_xz( uz_int, filter_xz )

     endif ! End of Pressure Step

     time_compute_poisson = (time_compute_poisson + &
                             stopwatch_get_elapsed( compute_start_time ))

     ! Step 3.0: Viscous solves.
     call stopwatch_tick( compute_start_time )

     ! Part 3: Viscous and Diffusion solves

     ! Viscous solve
     if ( solve_momentum_equation .eqv. .true. ) then
        call notify_cond( logflag, ' ' )
        call notify_cond( logflag, '      substep 3.1: solving viscous problem. ' )
        call notify_cond( logflag, ' ' )

        ! Step 3.1.1: Implicit Viscous solve.
        call solve_3D_viscous( ux_int, uy_int, uz_int )

        ! Step 3.1.2: Store velocity
        ux = ux_int
        uy = uy_int
        uz = uz_int

        ! Step 3.1.3: Filter the Solution
        call apply_filter_xz( ux, filter_xz  )
        call apply_filter_xz( uy, filter_xz  )
        call apply_filter_xz( uz, filter_xz  )

     endif ! End Viscous solves

     call stopwatch_tick( compute_start_time )
     time_compute_viscous = (time_compute_viscous + &
                             stopwatch_get_elapsed( compute_start_time ))

     ! Diffusion solve
     if ( use_scalar_diffusion .eqv. .true. ) then
        call notify_cond( logflag, ' ' )
        call notify_cond( logflag, '      substep 3.2: solving diffusive problem. ' )
        call notify_cond( logflag, ' ' )

        ! Step 3.2.1: Implicit diffusive solve
        call solve_3D_diffusion( rho_int )

        ! Step 3.2.2: Store density
        rho = rho_int

        ! Step 3.2.3: Filter the solution
        call apply_filter_xz( rho, filter_xz )

     else ! Unlike the momentum equation, we include the explicit AB3 for the scalar
          ! equation because we are also interested in solve scalar advection-diffusion
          ! problems.

        call notify_cond( logflag, ' ' )
        call notify_cond( logflag, '      substep 3: normalizing by time step coefficient ' )
        call notify_cond( logflag, '                 in diffusive solve                   ' )
        call notify_cond( logflag, ' ' )

       ! Step 3.2.1: Apply AB3 coefficient
       rho_int = rho_int / g0

       ! Step 3.2.2: Store density
       rho = rho_int

     endif ! End Diffusive solve

     ! XXXXXXXXX: need to tick a stopwatch for diffusion

     ! Step 3.3: Apply interfacial averaging.
     if ( do_interfacial_averaging) then
        call apply_interfacial_averaging( ux )
        call apply_interfacial_averaging( uy )
        call apply_interfacial_averaging( uz )
        call apply_interfacial_averaging( rho )
     endif

     ! Part 4: Advance Time, Compute Pressure Boundary Conditions

     ! Step 4.0: Time advance.
     call advance_time()

     ! Step 4.1: Compute the double-curl for the pressure boundary conditions.
     call compute_double_curl( ux0, uy0, uz0, Cux0, Cuy0, Cuz0 )


     ! Part 5: Cleanup and report errors

     ! Step 5.1: See how long it took for this timestep.
     elapsed_time_timestep = stopwatch_get_elapsed( clock_start )

     ! Step 5.2: Write errors and stats to file
     if ( rank == root ) then
        call write_diagnostics_timestep( time_ndx, elapsed_time_timestep )

        ! Update our running error summation.
        total_norm_numeric_error = total_norm_numeric_error + compute_numeric_error()

     endif

     ! Part 6: Update time, counters & write restart file

     ! Step 6.1: Update the simulation time.
     simulation_time = simulation_time + dt

     ! Step 6.2: Check if it is time to write out restart file to disk.
     if ( mod( time_ndx, timesteps_between_restarts ) == 0 ) then
        call notify_cond( logflag, ' ' )
        call notify_cond( logflag, '   Writing restart file ' )
        call notify_cond( logflag, ' ' )

        call stopwatch_tick( io_start_time )

        ! Create the restart file by concatenating the underscore and restart index
        write( restart_index_string, '(I6.6)' ) time_ndx / timesteps_between_restarts
        call open_restart_file( trim( restart_name ) // '_' // restart_index_string // '.h5', io_flag )
        call smpm_assert( io_flag, 'Failed to open ' // trim( restart_name ) // '_' // restart_index_string // '.h5' )

        ! Write the grid and data.
        call write_restart_grid
        call write_restart_constant_fields
        call write_restart_data

        ! Close the restart file.
        call close_restart_file

        ! Create an XDMF file describing the file.
        call write_xdmf_restart_file( trim( restart_name ) // '_' // restart_index_string //  '.xmf', &
                                      trim( restart_name ) // '_' // restart_index_string //  '.h5' )
        time_restart_io = (time_restart_io + &
                           stopwatch_get_elapsed( io_start_time ))
     endif

  enddo ! Time Loop

  ! Stop timing and display elapsed time.
  call CPU_TIME( CPUtimeEnd )
  CPUtimeElapsed = CPUtimeEnd - CPUtimeStart

  time_timesteps = stopwatch_get_elapsed( timesteps_start_time )

  ! Transfer the CPU times from each of the proceses back to the master to
  ! record how well each core was utilized.
  allocate( time_cpu(nprocs) )
  call MPI_GATHER( CPUtimeElapsed,   1, MPI_DOUBLE_PRECISION, &
                   time_cpu, 1, MPI_DOUBLE_PRECISION, &
                   root, MPI_COMM_WORLD, ierr )

  do ii = 0, nprocs - 1
     write( caststr, '(f12.2)' ) time_cpu(ii+1)
     write( procstr, '(I0)' ) ii
     call notify( 'Process ' // procstr // ' CPU time elapsed : '// trim( adjustl( caststr ) ) // ' seconds' )
  end do

  ! Compute the total, and per-timestep, elapsed wall-clock time.
  write( caststr, '(f12.2)' ) time_setup + time_timesteps
  call notify_cond( logflag, ' ' )
  call notify_cond( logflag, 'Wall time elapsed (est.): '// trim( adjustl( caststr ) ) // ' seconds' )
  call notify_cond( logflag, ' ' )

  write( caststr, '(f12.5) ' ) time_timesteps / Nt
  call notify_cond( logflag, 'Wall time per time-step : '// trim( adjustl( caststr ) ) // ' seconds' )
  call notify_cond( logflag, ' ' )

  if ( rank == root ) then
     ! Write out the execution statististcs and close the field file.
     ! XXX: should we write out all ranks' CPU time?
     ! XXX: how do we handle io vs restart io?  total run time?
     call write_diagnostics_execution_stats( time_cpu(1), time_field_io, &
                                       time_restart_io, &
                                       time_compute_diffusion, &
                                       time_compute_poisson, &
                                       time_compute_viscous, &
                                       time_setup + time_timesteps, &
                                       total_norm_numeric_error )
     call close_diagnostics_file

  endif

  deallocate( Ax )
  deallocate( time_cpu )

  call FFTW_FREE( p_real_buff_FFT )
  call FFTW_FREE( p_complex_buff_FFT )

  ! Shutdown the HDF5 library before exiting.
  call H5CLOSE_F( hdf5_status )

  ! Finalize MPI.
  call MPI_FINALIZE( ierr )

end program smpm_incompressible_navier_stokes
