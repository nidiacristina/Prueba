subroutine initialize_configuration( infile_name, status_flag )
! Initializes the solver's configuration from the contents of infile_name and
! returns a flag indicating whether configuration was successful or not.
! Parameters are read from the input file, optional parameters are defaulted,
! and internal parameters are derived.  After this routine returns all of the
! parameters necessary for setting up the solver's internal state are available.

  use constants, only: nsubx, nprocs

  implicit none

  character(len=*), intent(in) :: infile_name
  logical, intent(out)         :: status_flag

  character(len=256)           :: caststr

  status_flag = .true.

  ! Parse our configuration from the input file and initialize any optional
  ! parameters that were not provided.
  call read_input_file( infile_name, status_flag )
  call initialize_optional_inputs()

  ! Derive the remainder of our configuration that is not specified in the
  ! input file.
  call derive_configuration()

  ! Ensure that the grid and our MPI rank count are compatible.
  if ( (nprocs > 1) .and. (mod( nsubx, nprocs ) .ne. 0) ) then
     write( caststr, '(A,I0,A,I0,A)' ) &
          'ERROR: Incompatible execution parameters.  nprocs (', nprocs, &
          ') does not divide the number of sub-domains (', nsubx, ').'
     call notify( caststr )

     status_flag = .false.
     return
  endif

end subroutine initialize_configuration

subroutine read_input_file( fname, success_flag )
! Reads an _in file and parses it.

  use constants, only:    n, nsubx, nsuby, nsubz, nu, nu_d, rho_0
  use options, only:      adaptive_timestep, apply_restart, apply_sponge_layer, bc_flag_diffusion, bc_flag_lhsgmres, &
                          bc_flag_viscous, check_null_error, check_numerical_error, do_interfacial_averaging, &
                          exact_nullspace_projection, enforce_strong_velocity_bc, facrobin, facrobin_ppe, &
                          filter_order_xz, filter_order_y, fname_init, fname_restart, fname_runname, fname_setup, &
                          gmres_maxit_poisson, gmres_maxit_viscous, gmres_restart_poisson, gmres_restart_viscous, &
                          gmres_tol_poisson, gmres_tol_viscous, read_bcs_from_initfile, &
                          read_from_setupfile, setup_and_stop, solve_momentum_equation, &
                          timesteps_between_logs, timesteps_between_restarts, &
                          timesteps_between_writes, &
                          use_gravity_force, &
                          use_scalar_transport, &
                          use_scalar_diffusion, &
                          use_capacitance_preconditioner, use_deflation, &
                          write_to_setupfile
  use postprocessor,only: compute_pressure, compute_streamfunction, field_time
  use sponge, only:       left_fraction, right_fraction, sponge_layer_location, time_scale_in_x
  use timestepping, only: dt, tend

  implicit none

  character(len=*), intent(in) :: fname
  logical, intent(out)         :: success_flag

  character(len=132)           :: line_buffer
  character(len=132)           :: line_buffer2
  character(len=10)            :: line_number_buffer
  integer                      :: ii, readstat, pound_ndx, ndx

  ! Flags indicating whether we've read each of our mandatory parameters.
  ! None of these have sensible defaults and the user needs to be alerted if
  ! we don't read one of them.
  logical                      :: read_n, read_nsubx, read_nsuby, read_nsubz, read_nu, read_nu_d, read_rho_0, &
                                  read_bc_flag_diffusion, read_bc_flag_viscous, read_facrobin, read_facrobin_ppe, &
                                  read_filter_order_xz, read_filter_order_y, read_fname_runname, & 
                                  read_dt, read_tend

  read_n                 = .false.
  read_nsubx             = .false.
  read_nsuby             = .false.
  read_nsubz             = .false.
  read_nu                = .false.
  read_nu_d              = .false.
  read_rho_0             = .false.
  read_bc_flag_diffusion = .false.
  read_bc_flag_viscous   = .false.
  read_facrobin          = .false.
  read_facrobin_ppe      = .false.
  read_filter_order_xz   = .false.
  read_filter_order_y    = .false.
  read_fname_runname     = .false.
  read_dt                = .false.
  read_tend              = .false.

  ! Open the file for reading.
  open( 200, file=trim( fname ), action='read', iostat=readstat )

  call sync_flag( readstat )

  ! We assume that since we were able to open the file, everything is fine.
  ! Elsewhere we'll validate our input parameters and make a determination
  ! if we read enough to start the solver or not.
  success_flag = (readstat == 0)
  if ( success_flag .eqv. .false. ) then
     return
  end if

  ! Read the file line-by-line.
  do ii = 1, 1000

     ! Read next line.
     read( 200, '(A)', iostat=readstat ) line_buffer

     ! Left adjust and trim the buffer.
     line_buffer = trim( adjustl( line_buffer ) )

     ! If we're at the end of the file, break out of our file reading loop.
     if ( readstat < 0 ) then
        close( 200 )
        exit
     endif

     ! Ignore comment lines and blank lines.
     if ( index( line_buffer, "#" ) /= 1 .and. len( trim( adjustl( line_buffer ) ) ) > 0 ) then

        if ( len( line_buffer ) == len( trim( line_buffer ) ) ) then
           write( line_number_buffer, '(I10)' ) ii
           call notify( '  WARNING: Configuration likely truncated.  Please review line ' // adjustl( line_number_buffer ) )
        end if

        ! Strip out any whitespace.
        call strip_all_whitespace( line_buffer )

        ! Strip out any trailing comments.
        pound_ndx = index( line_buffer, "#" )
        if (pound_ndx /= 0 ) line_buffer = line_buffer(1:pound_ndx-1)

        ! To upper case.
        line_buffer2 = line_buffer
        call to_upper( line_buffer )

        ! Parse the string for data.
        !
        ! Note: We list the data in alphabetical order starting with
        !       the first module, constants, then proceding to options,
        !       sponge and ending with timestepping. This makes the
        !       list easier to scan by the user.
        !
        !      20 June 2017
        !      Gustavo Rivera w/ Greg Thomsen

           ! MODULE CONSTANTS

           ! This is the number of points per subdomain.
           if ( index( line_buffer, 'N=') ==  1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) n
              read_n = .true.
              cycle
           endif

           ! This is the number of subdomains in the x-coordinate.
           if ( index( line_buffer, 'NSUBX=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) nsubx
              read_nsubx = .true.
              cycle
           endif

           ! This is the number of points in the y-coordinate.
           if ( index( line_buffer, 'NSUBY=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) nsuby
              read_nsuby = .true.
              cycle
           endif

           ! This is the number of subdomains in the z-coordinate.
           if ( index( line_buffer, 'NSUBZ=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) nsubz
              read_nsubz = .true.
              cycle
           endif

           ! This is the kinematic viscosity. (units of Length^2/Time)
           if ( index( line_buffer, 'NU=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) nu
              read_nu = .true.
              cycle
           endif

           ! This is the diffusivity. (units of Length^2/Time)
           if ( index( line_buffer, 'NU_D=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) nu_d
              read_nu_d = .true.
              cycle
           endif

           ! This is the reference density. (units of Mass/Length^3)
           if ( index( line_buffer, 'RHO_0=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) rho_0
              read_rho_0 = .true.
              cycle
           endif


           ! MODULE OPTIONS

           ! This parameter turns on/off the adaptive timestepping.
           ! (i.e. on == 1, off == 0)
           if ( index( line_buffer, 'ADAPTIVE_TIMESTEP=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) adaptive_timestep
              cycle
           endif

           ! This parameter turns on/off the restart capability.
           ! (i.e. on == 1, off == 0)
           if ( index( line_buffer, 'APPLY_RESTART=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) apply_restart
              cycle
           endif

           ! This parameter turns on/off the sponge layer.
           ! (i.e. on == 1, off == 0)
           if ( index( line_buffer, 'APPLY_SPONGE_LAYER=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) apply_sponge_layer
              cycle
           endif

           ! This parameter specifies the boundary conditions at each
           ! boundary domain (in x-z) for the density equation.
           ! (i.e. 1 == Dirichlet, 2 == Neumann)
           if ( index( line_buffer, 'BC_DIFFUSION=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) bc_flag_diffusion
              read_bc_flag_diffusion = .true.
              cycle
           endif

           ! This parameter specifies the boundary conditions at each
           ! boundary domain (x-z) for the Pressure Poisson equation.
           ! (i.e. 1 == Dirichlet, 2 == Neumann)
           if ( index( line_buffer, 'BC_LHSGMRES=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) bc_flag_lhsgmres
              cycle
           endif

           ! This parameter specifies the boundary conditions at each
           ! boundary domain (x-z) for the viscous solve.
           ! (i.e. 1 == Dirichlet, 2 == Neumann)
           if ( index( line_buffer, 'BC_VISCOUS=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) bc_flag_viscous
              read_bc_flag_viscous = .true.
              cycle
           endif

           ! This parameter is for computing the error ||Ax - b||
           ! for the nullspace computation.
           ! (i.e. do == 1, don't == 0)
           if ( index( line_buffer, 'CHECK_NULL_ERROR=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) check_null_error
              cycle
           endif

           ! This parameter is for computing the error ||Ax - b||,
           ! for both viscous and pressure solve, at every timestep.
           ! (i.e. do == 1, don't == 0)
           if ( index( line_buffer, 'CHECK_NUMERICAL_ERROR=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) check_numerical_error
              cycle
           endif

           ! This parameter is for performing the interfacial averaging
           ! across element boundaries.
           ! (i.e. do == 1, don't == 0)
           if ( index( line_buffer, 'DO_INTERFACIAL_AVERAGING=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) do_interfacial_averaging
              cycle
           endif

           ! This parameter applies the exact nullspace projection from
           ! Joshi et al. 2016 in the x-z domain.
           ! (i.e. do == 1, don't == 0)
           ! Note: only works for 2D runs (i.e. nsuby = 1)
           if ( index( line_buffer, 'EXACT_NULLSPACE_PROJECTION=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) exact_nullspace_projection
              cycle
           endif

           ! This parameter is for strongly enforcing the boundary condition
           ! after the pressure solve.
           ! (i.e. do == 1, don't == 0)
           if ( index( line_buffer, 'ENFORCE_STRONG_VELOCITY_BC=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) enforce_strong_velocity_bc
              cycle
           endif

           ! The value of the penalty factor for the viscous solve.
           if ( index( line_buffer, 'FACROBIN=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) facrobin
              read_facrobin = .true.
              cycle
           endif

           ! The value of the penalty factor for the pressure solve.
           if ( index( line_buffer, 'FACROBIN_PPE=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) facrobin_ppe
              read_facrobin_ppe = .true.
              cycle
           endif

           ! The order of the filtering matrix in the x-z
           if ( index( line_buffer, 'FILTER_ORDER_XZ=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) filter_order_xz
              read_filter_order_xz = .true.
              cycle
           endif

           ! The order of the exponential filter in the y-coordinate.
           if ( index( line_buffer, 'FILTER_ORDER_Y=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) filter_order_y
              read_filter_order_y = .true.
              cycle
           endif

           ! The name of the HDF5 that contains the initial conditions.
           if ( index( line_buffer, 'FNAME_INIT=' ) == 1 ) then
              ndx        = index( line_buffer, '=' )
              fname_init = trim( adjustl( line_buffer2(ndx+1:) ) )
              cycle
           endif

           ! The name of the restart file that we will be reading the
           ! the restart field from.
           if ( index( line_buffer, 'FNAME_RESTART=' ) == 1 ) then
              ndx           = index( line_buffer, '=' )
              fname_restart = trim( adjustl( line_buffer2(ndx+1:) ) )
              cycle
           endif

           ! The name of the output file on which we write the field data.
           if ( index( line_buffer, 'FNAME_RUNNAME=' ) == 1 ) then
              ndx           = index( line_buffer, '=' )
              fname_runname = trim( adjustl( line_buffer2(ndx+1:) ) )
              read_fname_runname = .true.
              cycle
           endif

           ! The name of the setup file on which we store the poisson operators.
           if ( index( line_buffer, 'FNAME_SETUP=' ) == 1 ) then
              ndx         = index( line_buffer, '=' )
              fname_setup = trim( adjustl( line_buffer2(ndx+1:) ) )
              cycle
           endif

           ! The maximum number of iterations for the GMRES of the poisson solver.
           if ( index( line_buffer, 'GMRES_MAXIT_POISSON=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) gmres_maxit_poisson
              cycle
           endif

           ! The maximum number of iterations for the GMRES of the viscous &
           ! diffusive solver.
           if ( index( line_buffer, 'GMRES_MAXIT_VISCOUS=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) gmres_maxit_viscous
              cycle
           endif

           ! The maximum number of iterations for the GMRES of the poisson
           ! solver at restart.
           if ( index( line_buffer, 'GMRES_RESTART_POISSON=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) gmres_restart_poisson
              cycle
           endif

           ! The maximum number of iterations for the GMRES of the viscous
           ! and diffusive solver at restart.
           if ( index( line_buffer, 'GMRES_RESTART_VISCOUS=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) gmres_restart_viscous
              cycle
           endif

           ! The minimum tolerance to exit the GMRES of the poisson
           ! solver.
           if ( index( line_buffer, 'GMRES_TOL_POISSON=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) gmres_tol_poisson
              cycle
           endif

           ! The minimum tolerance to exit the GMRES of the viscous and
           ! diffusive solver.
           if ( index( line_buffer, 'GMRES_TOL_VISCOUS=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) gmres_tol_viscous
              cycle
           endif

           ! This parameters specifies whether we will read the boundary
           ! conditions from the init file or not.
           ! (i.e. do == 1, don't == 0)
           if ( index( line_buffer, 'READ_BCS_FROM_INITFILE=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) read_bcs_from_initfile
              cycle
           endif

           ! This parameter specifies whether we read the poisson operators
           ! from a setupfile or we build them from scratch.
           ! (i.e. do == 1, don't == 0)
           if ( index( line_buffer, 'READ_FROM_SETUPFILE=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) read_from_setupfile
              cycle
           endif

          ! This parameters specifies if we just want to setup the poisson
          ! operators, not setup and solve.
          ! (i.e. do == 1, don't == 0)
           if ( index( line_buffer, 'SETUP_AND_STOP=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) setup_and_stop
              cycle
           endif

           ! This parameter specifies whether we want to solve the momentum
           ! Equation or not.
           ! (i.e. do == 1, don't == 0)
           if ( index( line_buffer, 'SOLVE_MOMENTUM_EQUATION=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) solve_momentum_equation
              cycle
           endif

           ! This parameter specifies the number of timesteps needed to write
           ! the field data to the HDF5 *_out.h5
           if ( index( line_buffer, 'TIMESTEPS_BETWEEN_WRITES=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) timesteps_between_writes
              cycle
           endif

           ! This parameter specifies the number of timesteps needed to write
           ! the solver log. This is not the same as writing the statistics
           ! of the solver. The latter are written to the *_out.h5 file regardless.
           ! of the timestep.
           if ( index( line_buffer, 'TIMESTEPS_BETWEEN_LOGS=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) timesteps_between_logs
              cycle
           endif

           ! This parameter specifies the number of timesteps needed to write
           ! the restart field, in case we want to restart the run.
           if ( index( line_buffer, 'TIMESTEPS_BETWEEN_RESTARTS=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) timesteps_between_restarts
              cycle
           endif

           ! This parameter specifies whether we want to use the capacitance
           ! preconditioner or not during the poisson solve.
           ! (i.e. yes == 1, no == 0)
           if ( index( line_buffer, 'USE_CAPACITANCE_PRECONDITIONER=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) use_capacitance_preconditioner
              cycle
           endif

           ! This parameter specifies whether we want to use the deflation
           ! technique or not for the poisson solve.
           ! (i.e. yes == 1, no == 0)
           if ( index( line_buffer, 'USE_DEFLATION=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) use_deflation
              cycle
           endif

           ! This parameter specifies whether we want to have the density
           ! equation be a passive or active scalar. If we want an active
           ! scalar, then the density is fed to the momentum equation
           ! via the gravitational force.
           ! (i.e. for active scalar == 1, for passive scalar == 0)
           if ( index( line_buffer, 'USE_GRAVITY_FORCE=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) use_gravity_force
              cycle
           endif

           ! This parameter specifies whether we want to solve the
           ! scalar (density) diffusion equation or not.
           ! (i.e. do == 1, don't == 0)
           if ( index( line_buffer, 'USE_SCALAR_DIFFUSION=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) use_scalar_diffusion
              cycle
           endif

           ! This parameter specifies whether we want to solve the scalar
           ! advection or not.
           ! (i.e. do == 1, don't == 0).
           if ( index( line_buffer, 'USE_SCALAR_TRANSPORT=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) use_scalar_transport
              cycle
           endif

           ! This parameter specifies whether we want to write a setupfile or not.
           ! (i.e. do == 1, don't == 0)
           if ( index( line_buffer, 'WRITE_TO_SETUPFILE=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) write_to_setupfile
              cycle
           endif

           ! MODULE POSTPROCCESOR

           ! This parameter specifies whether we want to compute the streamfunction
           ! or not.
           ! (i.e. do == 1, don't == 0)
           if ( index( line_buffer, 'COMPUTE_STREAMFUNCTION=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) compute_streamfunction
              cycle
           endif

           ! This parameter specifies whether we want to compute the pressure field
           ! or not from the post processor.
           ! (i.e. do == 1, don't == 0)
           if ( index( line_buffer, 'COMPUTE_PRESSURE=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) compute_pressure
              cycle
           endif

           ! This is the simulation time during postprocessing.
           if ( index( line_buffer, 'FIELD_TIME=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) field_time
              cycle
           endif

           ! MODULE SPONGE

           ! This parameter corresponds to the fraction of the total domain length
           ! that will be occupied by the sponge, at the left boundary.
           if ( index( line_buffer, 'LEFT_FRACTION=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), *) left_fraction
              cycle
           endif

           ! This parameter corresponds to the fraction of the total domain length
           ! that will be occupied by the sponge, at the right boundary.
           if ( index( line_buffer, 'RIGHT_FRACTION=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), *) right_fraction
              cycle
           endif

           ! This parameter corresponds to the location of the sponge layer. It only
           ! works for either the left or right boundary.
           if ( index( line_buffer, 'SPONGE_LAYER_LOCATION=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) sponge_layer_location
              cycle
           endif

           ! This parameter corresponds to the time scale of the sponge layer.
           ! Further information can be found in Ammar's Thesis.
           if ( index( line_buffer, 'TIME_SCALE_SPONGE_IN_X=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) time_scale_in_x
              cycle
           endif

           ! MODULE TIMESTEPPING

           ! This is the timestep size.
           if ( index( line_buffer, 'DT=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) dt
              read_dt = .true.
              cycle
           endif

           ! This is the final time of the simulation, in seconds.
           if ( index( line_buffer, 'TEND=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) tend
              read_tend = .true.
              cycle
           endif

           ! If we get here, throw an error.
           write( line_number_buffer, '(I10)' ) ii
           call notify( 'Error: Unclassifiable statement at line ' // adjustl( line_number_buffer ) )
           call notify( line_buffer )
           success_flag = .false.

           close( 200 )
           return

     else

        ! Do nothing.  This line is either blank or #-commented.
        !
        ! XXX: This is broken.  If a user adds comments at the end of a line
        !      of meaningful input, this breaks.  I should figure out a way to
        !      test for this.

     endif

  enddo

  ! Check that we read each of the mandatory parameters and notify someone
  ! if we didn't.  We check every parameter so the caller has a complete
  ! idea of what was missing rather than the first thing we find.
  if ( read_n .eqv. .false. ) then
     call notify( 'Error: Failed to read ''N''.' )
     success_flag = .false.
  end if
  if ( read_nsubx .eqv. .false. ) then
     call notify( 'Error: Failed to read ''NSUBX''.' )
     success_flag = .false.
  end if
  if ( read_nsuby .eqv. .false. ) then
     call notify( 'Error: Failed to read ''NSUBY''.' )
     success_flag = .false.
  end if
  if ( read_nsubz .eqv. .false. ) then
     call notify( 'Error: Failed to read ''NSUBZ''.' )
     success_flag = .false.
  end if
  if ( read_nu .eqv. .false. ) then
     call notify( 'Error: Failed to read ''NU''.' )
     success_flag = .false.
  end if
  if ( read_nu_d .eqv. .false. ) then
     call notify( 'Error: Failed to read ''NU_D''.' )
     success_flag = .false.
  end if
  if ( read_rho_0 .eqv. .false. ) then
     call notify( 'Error: Failed to read ''RHO_0''.' )
     success_flag = .false.
  end if
  if ( read_bc_flag_diffusion .eqv. .false. ) then
     call notify( 'Error: Failed to read ''BC_DIFFUSION''.' )
     success_flag = .false.
  end if
  if ( read_bc_flag_viscous .eqv. .false. ) then
     call notify( 'Error: Failed to read ''BC_VISCOUS''.' )
     success_flag = .false.
  end if
  if ( read_facrobin .eqv. .false. ) then
     call notify( 'Error: Failed to read ''FACROBIN''.' )
     success_flag = .false.
  end if
  if ( read_facrobin_ppe .eqv. .false. ) then
     call notify( 'Error: Failed to read ''FACROBIN_PPE''.' )
     success_flag = .false.
  end if
  if ( read_filter_order_xz .eqv. .false. ) then
     call notify( 'Error: Failed to read ''FILTER_ORDER_XZ''.' )
     success_flag = .false.
  end if
  if ( read_filter_order_y .eqv. .false. ) then
     call notify( 'Error: Failed to read ''FILTER_ORDER_Y''.' )
     success_flag = .false.
  end if
  if ( read_fname_runname .eqv. .false. ) then
     call notify( 'Error: Failed to read ''FNAME_RUNNAME''.' )
     success_flag = .false.
  end if
  if ( read_dt .eqv. .false. ) then
     call notify( 'Error: Failed to read ''DT''.' )
     success_flag = .false.
  end if
  if ( read_tend .eqv. .false. ) then
     call notify( 'Error: Failed to read ''TEND''.' )
     success_flag = .false.
  end if

  close( 200 )

end subroutine read_input_file

subroutine initialize_optional_inputs
! Initializes input parameters that are considered optional and need a default
! value based on the input file read.

  use options, only: fname_init, fname_restart, fname_runname, fname_setup, &
                     gmres_restart_poisson, gmres_restart_viscous, &
                     gmres_maxit_poisson, gmres_maxit_viscous

  implicit none

  ! If the restart parameters were unset, set them to defaults.
  if ( gmres_restart_poisson == -1 ) then
     gmres_restart_poisson = gmres_maxit_poisson
  endif
  if ( gmres_restart_viscous == -1 ) then
     gmres_restart_viscous = gmres_maxit_viscous
  endif

  ! Derive file names from the specified run name if they weren't specified.
  if ( len( trim( fname_init ) ) == 0 ) then
     fname_init = trim( fname_runname ) // '_init.h5'
  endif
  if ( len( trim( fname_restart ) ) == 0 ) then
     fname_restart = trim( fname_runname ) // '_restart.h5'
  endif
  if ( len( trim( fname_setup ) ) == 0 ) then
     fname_setup = trim( fname_runname ) // '_setup'
  endif

end subroutine initialize_optional_inputs

subroutine derive_configuration()
! Derives solver parameters from those specified by the user via the command
! line and input file.

  use constants, only:         n, nky, nprocs, nsg, nsubx, nsuby, nsubz, r
  use options, only:           bc_flag_viscous, bc_flag_viscous_x, &
                               bc_flag_viscous_y, bc_flag_viscous_z
  use woodbury_matrices, only: dimA, dimC, dimblock, k, numA_per_rank, &
                               numblocks, rpk, s

 implicit none

  ! Configure the boundary conditions for the vector viscous equation.
  !   The following is the bc configuration:
  !   [ 1, 2, 3, 4 ] = [bottom wall, right wall, top wall, left wall]
  !           w/
  !   dirichlet = 1
  !   neumann   = 2
  ! Note: Now that we have incorporated a transverse velocity, we must
  ! remember that this velocity runs parallel to the walls and is
  ! periodic in the transverse. Therefore, it will obey the same
  ! boundary condition as that of the corresponding tangential
  ! velocity in a 2D case.

  bc_flag_viscous_x = 1
  bc_flag_viscous_y = 1
  bc_flag_viscous_z = 1
  if ( bc_flag_viscous(1) == 2 ) then
     bc_flag_viscous_x(1) = 2
     bc_flag_viscous_y(1) = 2
  endif
  if ( bc_flag_viscous(2) == 2 ) then
     bc_flag_viscous_y(2) = 2
     bc_flag_viscous_z(2) = 2
  endif
  if ( bc_flag_viscous(3) == 2 ) then
     bc_flag_viscous_x(3) = 2
     bc_flag_viscous_y(3) = 2
  endif
  if ( bc_flag_viscous(4) == 2 ) then
     bc_flag_viscous_y(4) = 2
     bc_flag_viscous_z(4) = 2
  endif

  ! Compute some secondary input variables.
  nsg = n * n * nsubx * nsubz          ! global number of GLL points
  r   = n * n * nsubx * nsubz          ! total number of grid points.
  rpk = n * n * nsubx * nsubz / nprocs ! number of grid points per rank.

  ! Compute some variables for the Woodbury domain decomposition.
  numA_per_rank = nsubx / nprocs
  numblocks     = nprocs
  dimblock      = nsubz * n * n * numA_per_rank
  s             = nsubz * n
  k             = 2 * s * (nsubx - 1)

  ! Compute the number of transverse wavenumbers we're working with.
  if ( nsuby .gt. 1 ) then
     nky           = nsuby / 2 + 1
  else
     nky           = 1
  endif

  ! Set some domain decomposition constants.
  dimA          = n * n * nsubz
  dimC          = 2 * s * (nsubx - 1)

  ! Calculate the indices delineating the ranges of various matrices owned by
  ! this rank.
  call set_ownership_ranges()

end subroutine derive_configuration

subroutine strip_all_whitespace( str )

  implicit none

  character(len=*), intent(inout) :: str
  character(len=len( str ))       :: str2
  integer                         :: ii, counter

  counter = 1
  do ii = 1, len( str )
     if ( str(ii:ii) .ne. ' ' ) then
        str2(counter:counter) = str(ii:ii)
        counter               = counter + 1
     endif
  enddo

  str2(counter:) = ' '
  str            = str2

end subroutine strip_all_whitespace

subroutine to_upper( strIn )
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)

     implicit none

     character(len=*), intent(inout) :: strIn
     character(len=len(strIn))       :: strOut
     integer                         :: i,j

     do i = 1, len( strIn )
          j = iachar( strIn(i:i) )
          if (j >= iachar( "a" ) .and. j <= iachar( "z" ) ) then
               strOut(i:i) = achar( iachar( strIn(i:i) ) - 32 )
          else
               strOut(i:i) = strIn(i:i)
          end if
     end do

     strIn = strOut

end subroutine to_upper
