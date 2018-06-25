subroutine compute_poisson_kernel
! Computes the left kernel vectors of the Poisson operator and the capacitance
! matrix using shifted inverse iteration.
!
! This subroutine will automatically adjust the spectral shift as well as
! determine the stopping condition from properties of the right kernel space
! which is assumed to be the constant vector.

  use constants, only:               nsg, sigma
  use options, only:                 check_null_error, gmres_maxit_poisson, gmres_restart_poisson, &
                                     gmres_tol_poisson, use_capacitance_preconditioner
  use parallel_linear_algebra, only: pdot_product, pnorm2
  use precision, only:               dp
  use woodbury_matrices, only:       dimCk, k, rpk, uC, uC_right, uL

  implicit none

  real(kind=dp)                            :: gmres_tol, tol, err, errC
  integer                                  :: maxit, gmres_iterations, gmres_restart, ii
  real(kind=dp), allocatable, dimension(:) :: uii, Ax, vC
  character(len=32)                        :: caststr
  real(kind=dp), allocatable, dimension(:) :: tmp
  real(kind=dp), allocatable, dimension(:) :: ones
  real(kind=dp), allocatable, dimension(:) :: BTuC
  real(kind=dp), allocatable, dimension(:) :: errC_array

  external                                 :: apply_poisson_capacitance_shift
  external                                 :: apply_preconditioned_capacitance_shift

  ! Set some constants.
  sigma            = 1.0e-7_dp
  maxit            = 20
  tol              = 1.0e-12_dp
  gmres_iterations = 2 * gmres_maxit_poisson

  ! Allocate the array to store the error.
  allocate( errC_array(1:maxit) )

  ! Notify.
  call notify( 'Computing the left kernel vector of the poisson capacitance problem. ' )
  call notify( ' ')

  ! Allocate scratch space for computing the kernels.
  !
  ! NOTE: These are the sized large enough to hold the biggest computations
  !       encountered in this routine.  When a smaller computation is
  !       required, that does not implicitly reference a subset of the vector,
  !       explicit indexing is used.
  allocate( ones(1:rpk), Ax(1:rpk), tmp(1:k) )

  ! Get the right nullity of the full operator.
  ones = 1.0_dp / sqrt( real( nsg, kind=dp ) )
  call apply_smpm_poisson( Ax, ones )
  write( caststr, '(D17.10)' ) sqrt( pdot_product( Ax, Ax ) )
  call notify( '      Right null singular value of L is : ' // caststr )

  ! Set the tolerance of the LNSV computation to be equal to the computed
  ! singular value as determined by computing ||L * 1 ||.
  tol   = pnorm2( Ax )
  sigma = tol
  deallocate( ones )

  ! Get the right nullity of the capacitance operator.
  call compute_schur_right_nullity()

  deallocate(tmp)
  allocate( tmp(1:dimCk) )
  tmp = uC_right

  deallocate( Ax )
  allocate( Ax(1:dimCk) )
  Ax = 0.0_dp

  call apply_poisson_capacitance_sparse( Ax, tmp )

  allocate( vC(1:dimCk) )
  vC = tmp / pnorm2( tmp )

  write( caststr, '(D17.10)' ) pnorm2( Ax )
  call notify( '      Right null singular value of C is : ' // caststr )
  call notify( ' ' )

  ! Set the tolerance to the right nullity of the capacitance problem.
  tol   = pnorm2( Ax )
  sigma = max( 10.0_dp * tol, 1.0d-6 )

  ! Initialize.
  call random_number( uC )
  uC = 1.0_dp
  uC = uC / pnorm2( uC )

  ! Ensure we have space to copy uC into.
  allocate( uii(1:dimCk) )
  deallocate( tmp )
  allocate( tmp(1:rpk) )

  ! Allocate an intermediate array.
  allocate( BTuC(1:rpk) )

  ! Start looping the shift and invert.
  call notify('Starting inverse iteration.' )
  call notify(' ' )
  do ii = 1, maxit

     ! Dynamically adjust the shift if we've stagnanted (i.e. converged to the
     ! wrong singular value).
     if ( ii > 2 ) then
        if ( errC_array(ii-2) / errC_array(ii-1) < 10.0_dp ) then
           sigma = sigma / 10.0_dp
           call notify( '   Inverse iteration stagnated.  Reducing the magnitude of the shift. ' )
           call notify( ' ' )
        endif
     endif

     ! Notify user that we're starting this iteration.
     write( caststr, '(I10)' ) ii
     call notify( '      Starting inverse iteration number ' // caststr )
     write( caststr, '(D17.10)' ) sigma
     call notify( '         shift/sigma: ' // caststr )

     ! Solve the shifted system.
     uii                   = uC
     gmres_iterations = min( 2 * gmres_maxit_poisson, k - 1 )
     gmres_tol        = min( tol / 10.0_dp, gmres_tol_poisson )
     gmres_restart    = min( 2 * gmres_restart_poisson, gmres_iterations )
     if ( use_capacitance_preconditioner ) then
        call compute_gmres_householder( uii, uC, dimCk, gmres_tol, gmres_iterations, &
                                        gmres_restart, apply_preconditioned_capacitance_shift )
        call apply_capacitance_preconditioner_transpose( uC , uii )
        uii = uC
     else
        call compute_gmres_householder( uii, uC, dimCk, gmres_tol, gmres_iterations, &
                                        gmres_restart, apply_poisson_capacitance_shift )
     endif

     ! Normalize.
     uC = uii / pnorm2( uii )

     ! Build the left kernel vector of L if asked to do so.
     if ( check_null_error ) then

        ! Apply B' to the kernel vector of C.
        BTuC = 0.0_dp
        call apply_pB_poisson_transpose( uC, BTuC )
        tmp = BTuC

        ! Divide A' into this product.
        call solve_A( tmp, uL, 1, 'T', cmplx( 0.0_dp, 0.0_dp, kind=dp ) )

        ! Normalize.
        uL = uL / pnorm2( uL )

     endif

     ! Check the error.
     if ( check_null_error ) then
        deallocate( Ax )
        allocate( Ax(1:rpk) )
        Ax = 0.0_dp
        call apply_smpm_poisson_transpose_parallel( Ax, uL )
        err = pnorm2( Ax )
        deallocate( Ax )
        allocate( Ax(1:dimCk) )
     endif

     ! Check the error.
     Ax = 0.0_dp
     call apply_poisson_capacitance_shift( Ax, uC )
     Ax   = Ax + sigma * uC
     errC = pnorm2( Ax )

     ! Store the error for dynamically adjusting the shift.
     errC_array(ii) = errC

     ! Do some notifications.
     !call notify( '      Inverse iteration number ' // caststr )
     if ( check_null_error ) then
        write( caststr, '(D17.10)' ) err
        call notify( '         ||u^TL||   : ' // caststr )
     end if
     write( caststr, '(D17.10)' ) errC
     call notify( '         ||u^TC||   : ' // caststr )
     write( caststr, '(I10)' ) gmres_iterations
     call notify( '         GMRES iter : '// caststr )
     write( caststr, '(D17.10)' ) gmres_tol
     call notify( '         GMRES error: ' // caststr )
     call notify( '  ' )

     ! If we're done, exit the loop.
     if ( errC < tol ) then
        exit
     endif

  enddo

  ! Report the error.
  write( caststr, '(D17.10)') errC
  call notify( '   Inverse iteration converged with residual: ' // caststr )

  ! Build the left kernel vector of L.
  call notify( '   Building the null vector of the Poisson operator.' )

  ! Apply B' to the kernel vector of C.
  call apply_pB_poisson_transpose( uC, uL )

  ! Divide A' into that product.
  BTuC = uL
  call solve_A( BTuC, uL, 1, 'T', cmplx( 0.0_dp, kind=dp ) )

  ! Normalize.
  uL = uL / pnorm2( uL )
  call notify( ' ' )

  deallocate( Ax, tmp )

end subroutine compute_poisson_kernel

subroutine compute_schur_right_nullity( )
! Compute the right null vector of the Schur matrix for the zero wavenumber.

  use constants, only:               nsg
  use parallel_linear_algebra, only: pnorm2
  use precision, only:               dp
  use woodbury_matrices, only:       rpk, uC_right

  implicit none

  real(kind=dp), allocatable, dimension(:) :: vC, vC_x, vC_z

  allocate( vC(1:rpk), vC_x(1:rpk), vC_z(1:rpk) )

  vC   = 1.0_dp / sqrt( real( nsg, kind=dp ) )
  vC_x = 0.0_dp
  vC_z = 0.0_dp

  call compute_gradient( vC, vC_x, vC_z )

  uC_right = 0.0_dp
  call apply_pB_poisson( vC, uC_right, vC_x, vC_z, 1 )

  uC_right = uC_right / pnorm2( uC_right )

  deallocate( vC, vC_x, vC_z )

end subroutine compute_schur_right_nullity
