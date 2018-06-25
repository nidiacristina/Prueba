subroutine setup_poisson_operators( read_setup_file, write_setup_file, fname_setup, success_flag, elapsed_time_null_basis )
! Sets up the Poisson capacitance and deflation operators, as well as computes
! the kernel vectors.  Since this is an expensive operation, this routine can
! optionally read previously-assembled operators from disk and/or write them to
! disk for future runs, which is useful for debugging, benchmarking,
! optimization, etc.

  use precision, only:         dp, int64
  use options, only:           use_capacitance_preconditioner, use_deflation
  use woodbury_matrices, only: C_poisson_block, S_poisson

  implicit none

  logical, intent(in)        :: read_setup_file
  logical, intent(in)        :: write_setup_file
  character(*), intent(in)   :: fname_setup
  logical, intent(out)       :: success_flag
  real(kind=dp), intent(out) :: elapsed_time_null_basis

  ! NOTE: We explicitly request 8-byte integers to ensure the highest
  !       precision timer available is used by system_clock().
  integer(kind=int64)        :: clock_start, clock_stop, clock_rate

  if ( read_setup_file ) then
     call notify( 'Reading from setup file: ' // fname_setup )

     ! Read from the setup file.
     call read_setupfile( trim( fname_setup ), success_flag )
     if ( success_flag .eqv. .false. ) then
        call notify( 'Failed to read the setup file.' )
        return
     end if

     ! Setup the Schur matrix for the zero wavenumber.
     C_poisson_block = S_poisson(:, :, :, 1)

     call compute_schur_right_nullity()

     ! NOTE: This isn't technically correct since we compute the right kernel
     !       vector from what we've read from disk, though this is only a tiny
     !       faction of what compute_poisson_kernel() computes.  Rather than
     !       provide meaningless information, we use zero as an indicator that
     !       we loaded it from disk.
     elapsed_time_null_basis = 0.0_dp

     ! XXX: Check that the coarse matrix and its singular vectors are correctly
     !      being read.
  else

     ! Start the fast assembly routine.
     call notify( 'Beginning Poisson matrix assembly.' )
     call assemble_laplacian

     ! Set up and factor the necessary matrices for the Woodbury matrix
     ! decomposition.
     call notify( 'Beginning solver setup. ' )
     call setup_3D_poisson_solver

     ! Note: The following two processes apply to the computation of the
     !       left kernel vector. For the analogous pertaining to the
     !       setup of the poisson matrix, refer to
     !       setup_3D_schur_preconditioner and setup_3D_deflation.
     ! Setup the capacitance preconditioner.
     if ( use_capacitance_preconditioner .or. use_deflation ) then
        use_capacitance_preconditioner = .true.
        call notify( 'Setting up capacitance preconditioner.' )
        call setup_capacitance_preconditioner
     endif

     ! Setup the deflation operator.
     if ( use_deflation ) then
        call notify( 'Setting up the deflation preconditioner.' )
        call setup_deflation
     endif

     call SYSTEM_CLOCK( COUNT_RATE=clock_rate )
     call SYSTEM_CLOCK( COUNT=clock_start )

     ! Compute the poisson kernel vectors with inverse iteration, unless we
     ! already read it from the setup file.
     call compute_poisson_kernel

     call SYSTEM_CLOCK( COUNT=clock_stop )
     elapsed_time_null_basis = real( (clock_stop - clock_start) , kind=int64 ) / clock_rate

  end if

  ! If asked to write to a setup file, write out the assembled matrices and
  ! their factorizations to disk.
  if ( write_setup_file ) then
     call notify( 'Writing to setup file ' // fname_setup )
     call write_setupfile( trim( fname_setup ), success_flag )
     call notify_cond( .not. success_flag, 'Failed to write the setup file.' )
  endif

end subroutine setup_poisson_operators
