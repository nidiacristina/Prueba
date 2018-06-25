program smpm_validator
! Validation tool for the SMPM code base.  One or more analytic validation
! routines are executed after initializing the SMPM solver's state from an
! input file.  On success, the validator returns 0.  On failure, the non-zero
! value is indicates what went awry: either an error in setting up 
!
! NOTE: This program is F2008 due to its use of features from the F2003 and
!       F2008 standards (e.g. command line parsing and return codes).

  use command_line, only:      get_command_line_double, get_command_line_string, &
                               get_command_line_strings
  use constants, only:         nprocs, rank, root
  use HDF5, only:              h5close_f, h5open_f
  use mpi, only:               MPI_COMM_WORLD, MPI_DOUBLE_PRECISION
  use precision, only:         dp
  use validation

  implicit none

  ! Command line parsing.
  integer                                     :: n_arg
  character(len=:), allocatable               :: infile_name
  real(kind=dp)                               :: validation_tolerance
  character(len=:), allocatable               :: validation_tolerance_string
  character(len=:), allocatable, dimension(:) :: validator_names

  ! Status variables for MPI and HDF5 calls.
  integer                                     :: ierr
  integer                                     :: hdf5_status

  ! Flag indicating whether the SMPM code base has been validated against the
  ! supplied tolerance.
  logical                                     :: success_flag = .true.

  ! Index for iterating through requested tests.
  integer                                     :: test_index

  interface
     subroutine validate_test_names( names, success_flag )
       character(len=:), allocatable, dimension(:), intent(in) :: names
       logical, intent(out) :: success_flag
     end subroutine validate_test_names

     subroutine log_validation_state( tolerance_string, validator_names )
       character(len=:), allocatable, intent(in)               :: tolerance_string
       character(len=:), allocatable, dimension(:), intent(in) :: validator_names
     end subroutine log_validation_state
  end interface

  ! Initialize MPI.
  call MPI_INIT( ierr )

  ! Get my rank and total number of processes.
  call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

  n_arg = COMMAND_ARGUMENT_COUNT()

  if (n_arg < 3) then
     ! No arguments means the user asked for an enumeration of the available
     ! tests.
     if (n_arg == 0) then
        if (rank == root) then
           write(*,'(A,I1,A)') 'Available tests:'
           write(*,*)          ''
           do test_index = 1, size( validation_routines )
              write(*,*)       '  ', validation_routines(test_index)
           end do
           write(*,*)          ''
        end if
        ! Asking for help is not an error.
        test_index = 0
     else
        if (rank == root) then
           write(*,'(A,I1,A)') 'ERROR: Expected at least 3 input arguments though received ', n_arg, '.'
           write(*,*)          ''
           write(*,'(A)')      '  Usage: smpm_vaidator <runname_in> <tolerance> <validator> [...]'
           write(*,*)          ''
        end if

        ! Failure to setup error.
        test_index = -1
     end if

     call MPI_FINALIZE( ierr )
     if (test_index /= 0) then
        call stop_wrapper( test_index )
     else
        stop
     end if

  else
     ! Get the input file, tolerance, and a list of validators to run.
     call get_command_line_string( 1, infile_name )
     call get_command_line_double( 2, validation_tolerance, success_flag )
     call get_command_line_string( 2, validation_tolerance_string )
     call get_command_line_strings( 3, n_arg, validator_names )

     call smpm_assert_code( success_flag, -1, 'Failed to acquire a valid comparison tolerance ("' // &
                                               validation_tolerance_string // ')".' )
  endif

  ! Initialize the HDF5 library.
  call H5OPEN_F( hdf5_status )

  ! See if we can take an early out if an unknown test case was requested.
  ! This skips a potentially lengthy setup time if we already know we can't
  ! succeed.
  call validate_test_names( validator_names, success_flag )
  call smpm_assert_code( success_flag, -1, 'Unknown validation routine requested.' )

  ! Keep track of what we've been requested to do.  We do this before
  ! initializing anything in case there is a problem and to keep it from being
  ! lost in the volume of log messages.
  call log_validation_state( validation_tolerance_string, validator_names )

  ! Initialize the solver state.
  call setup_smpm_state( infile_name, success_flag )
  call smpm_assert_code( success_flag, -1, 'Failed to setup the SMPM state.' )

  ! Walk through each of the tests requested by the user and see if we know
  ! what they're asking about.
  do test_index = 1, size( validator_names )
     if     ( validator_names(test_index) == validation_routines(ADVECTION_TERMS) ) then
        call validate_advection_terms( validation_tolerance, success_flag )

     elseif ( validator_names(test_index) == validation_routines(DIVERGENCE) ) then
        call validate_divergence( validation_tolerance, success_flag )

     elseif ( validator_names(test_index) == validation_routines(DOUBLE_CURL) ) then
        call validate_double_curl( validation_tolerance, success_flag )

     elseif ( validator_names(test_index) == validation_routines(FILTER_Y) ) then
        call validate_filter_y( validation_tolerance, success_flag )

     elseif ( validator_names(test_index) == validation_routines(FFTW) ) then
        call validate_fftw( validation_tolerance, success_flag )

     elseif ( validator_names(test_index) == validation_routines(GRADIENTS) ) then
        call validate_gradients( validation_tolerance, success_flag )

     elseif ( validator_names(test_index) == validation_routines(HELMHOLTZ_2D_SCALAR_DIRICHLET) ) then
        call validate_helmholtz_2D_scalar_dirichlet( validation_tolerance, success_flag )

     elseif ( validator_names(test_index) == validation_routines(HELMHOLTZ_3D_SCALAR_DIRICHLET) ) then
        call validate_helmholtz_3D_scalar_dirichlet( validation_tolerance, success_flag )

     elseif ( validator_names(test_index) == validation_routines(HELMHOLTZ_3D_SCALAR_NEUMANN) ) then
        call validate_helmholtz_3D_scalar_neumann( validation_tolerance, success_flag )

     elseif ( validator_names(test_index) == validation_routines(HELMHOLTZ_3D_VECTOR_DIRICHLET) ) then
        call validate_helmholtz_3D_vector_dirichlet( validation_tolerance, success_flag )

     elseif ( validator_names(test_index) == validation_routines(HELMHOLTZ_3D_VECTOR_NEUMANN) ) then
        call validate_helmholtz_3D_vector_neumann( validation_tolerance, success_flag )

     elseif ( validator_names(test_index) == validation_routines(LAPLACIAN) ) then
        call validate_laplacian( validation_tolerance, success_flag )

     elseif ( validator_names(test_index) == validation_routines(POISSON_3D_ITERATIVELY) ) then
        call validate_poisson_3D_iteratively( validation_tolerance, success_flag )

     elseif ( validator_names(test_index) == validation_routines(TRANSPORT_TERM) ) then
        call validate_transport_term( validation_tolerance, success_flag )

     end if

     if ( success_flag ) then
        call notify( 'The "' // trim( validator_names(test_index) ) // &
                     '" test achieved a tolerance of ' // validation_tolerance_string // &
                     '.' )
     else
        call notify( 'The "' // trim( validator_names(test_index) ) // &
                     '" test failed to meet a tolerance of ' // validation_tolerance_string // &
                     '.' )
        exit
     end if
  end do

  ! Shutdown the HDF5 library before exiting.
  call H5CLOSE_F( hdf5_status )

  ! Finalize MPI.
  call MPI_FINALIZE( ierr )

  ! Indicate what test case failed to validate.
  if ( success_flag ) then
     test_index = 0
  end if

  ! Most Fortran run-times will print out "STOP N" from each process that issues
  ! a stop, so we only do it from the root rank when it would be an abnormal
  ! termination.  This reduces the overall clutter in our output and avoids a
  ! "STOP 0" when all of the validation tests succeed.
  if (( rank == root ) .and. (test_index > 0)) then
     call stop_wrapper( test_index )
  end if

end program smpm_validator

subroutine setup_smpm_state( infile_name, success_flag )
! Sets up the SMPM solver's state from the supplied input file.  Returns its
! success or failure to the caller rather than exiting.

  use field_variables, only:   rho_bar, rho_bar_z, rho_bar_zz
  use options, only:           apply_sponge_layer, exact_nullspace_projection, &
                               fname_setup, read_from_setupfile, &
                               write_to_setupfile
  use precision, only:         dp
  use woodbury_matrices, only: rpk

  implicit none

  character(*), intent(in)                 :: infile_name
  logical, intent(out)                     :: success_flag

  ! Flag indicating success of file I/O.
  real(kind=dp)                            :: elapsed_time_null_basis

  ! Used to compute the background stratification.
  real(kind=dp), allocatable, dimension(:) :: junk

  call initialize_configuration( infile_name, success_flag )
  if ( .not. success_flag ) then
     call notify( "Failed to read the input file." )
     return
  end if

  ! Log the input parameters used.
  call display_inputs()

  ! Allocate space for our internal state.
  call allocate_solver_variables()

  ! Acquire our initial conditions.
  call setup_initial_conditions()

  ! Configure the deformation maps.
  call setup_deformation_maps()

  ! Setup the sponge layer if needed.
  if ( apply_sponge_layer ) then
     call setup_sponge_layer
  end if

  call setup_transverse()

  ! Compute the derivatives of the background stratification.
  allocate( junk(1:rpk) )
  call compute_gradient( rho_bar,   junk, rho_bar_z )
  call compute_gradient( rho_bar_z, junk, rho_bar_zz )
  deallocate( junk )

  ! Setup the Poisson capacitance and deflation operators, along with the null
  ! vectors.  These will be read from disk if they were created in a previous
  ! run, otherwise computed.
  call setup_poisson_operators( read_from_setupfile, write_to_setupfile, &
                                fname_setup, success_flag, elapsed_time_null_basis )
  call smpm_assert( success_flag, 'Failed to read or write ' // fname_setup )

  ! If asked, compute the nullspace basis.
  if ( exact_nullspace_projection ) then
     call setup_nullspace_projection()
  endif

end subroutine setup_smpm_state

subroutine validate_test_names( names, success_flag )
! Takes an array of test names and indicates whether all of them are known
! validation routines or not.

  use validation, only: validation_routines

  implicit none

  character(len=:), allocatable, dimension(:), intent(in) :: names
  logical, intent(out)                                    :: success_flag

  ! Indices in the names and validation_routines vectors, respectively.
  integer                                                 :: test_index
  integer                                                 :: validator_index

  ! Flag indicating if a given test name is known.
  logical                                                 :: known_test

  ! Ensure the tests requested are ones we know about.  We check this up front
  ! rather than below so we can exit early if any unknown tests are requested.
  ! That is particularly helpful when the solver configuration takes a while
  ! to setup.
  do test_index = 1, size( names )
     known_test = .false.
     do validator_index = 1, size( validation_routines )
        known_test = known_test .or. (names(test_index) == &
                                      validation_routines(validator_index))
     end do

     if ( .not. known_test ) then
        success_flag = .false.
        call notify( '"' // trim( names(test_index) ) // '" is not a valid routine.' )
        return
     end if
  end do

end subroutine validate_test_names

subroutine log_validation_state( tolerance_string, validator_names )
! Takes a tolerance string and an array of test names and logs information
! relevant to the validator's configuration.

  implicit none

  character(len=:), allocatable, intent(in)               :: tolerance_string
  character(len=:), allocatable, dimension(:), intent(in) :: validator_names

  integer                                                 :: case_index
  character(len=8)                                        :: case_index_string

  call notify( "Executing the following tests with a tolerance of " // tolerance_string // ":" )
  call notify( "" )

  do case_index = 1, size( validator_names )
     write( case_index_string, "(I0)" ) case_index
     call notify( "  " // trim( case_index_string ) // " - " // validator_names(case_index) )
  end do

  call notify( "" )

end subroutine log_validation_state
