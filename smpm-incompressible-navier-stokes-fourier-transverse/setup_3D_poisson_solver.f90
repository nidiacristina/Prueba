subroutine setup_3D_poisson_solver
! This routine drives the assembly of A and then builds the Schur blocks for
! ever wavenumber.

  use constants, only:             nky
  use options, only:               use_capacitance_preconditioner, use_deflation
  use precision, only:             dp
  use transverse, only:            qy
  use woodbury_matrices, only:     C_poisson_block, S_poisson

  implicit none

  ! Internal variables.
  integer                            :: ii
  character(len=32)                  :: caststr
  real(kind=dp)                      :: percent

  ! Notify.
  call notify('Assembling/factoring block-diagonal part of the Poisson operator.')

  ! Build the block-diagonal component of the Poisson operator.
  call setup_A_matrix()

  ! Build the Schur matrices for all the wavenumbers.
  call notify('Assembling the Poisson-Schur matrix for each wavenumber.')
  do ii = 1, nky

     ! Assemble this Schur matrix.
     call setup_2D_schur( S_poisson(:, :, :, ii), cmplx( qy(ii)**2, kind=dp ) )

     ! Notify the user of progress.
     percent = 100.0_dp * real( ii, kind=dp ) / real( nky, kind=dp )
     write( caststr, '(f6.1)' ) percent
     call notify('   Completed ' // trim(caststr) // ' percent of Schur assembly.' )

  enddo

  ! Notify.
  call notify( 'Poisson solver setup complete.' )

  ! For now, store the zero wavenumber's Schur matrix into another array for
  ! poisson kernel computation.  This is just so that we can keep using the
  ! same framework to solve the kernel problem.
  C_poisson_block = S_poisson(:, :, :, 1)

  ! Setup the Schur block-Jacobi preconditioner.
  if ( use_capacitance_preconditioner .or. use_deflation ) then
     use_capacitance_preconditioner = .true.
     call notify( 'Setting up 3D Block-Jacobi-Schur preconditioner.' )
     call setup_3D_schur_preconditioner()
  endif

  ! Setup the deflation matrices for the Schur deflation method.
  if ( use_deflation ) then
     call notify( 'Setting up the 3D deflation-Schur preconditioner.' )
     call setup_3D_deflation()
  endif

end subroutine setup_3D_poisson_solver
