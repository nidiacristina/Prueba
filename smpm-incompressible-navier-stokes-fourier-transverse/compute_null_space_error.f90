subroutine compute_null_space_error( )
! Computes the null space error of the Poisson kernel and the Schur kernel.

  use constants, only:               nsuby, sigma
  use errors, only:                  l2_poisson_kernel_error, &
                                     l2_poisson_schur_kernel_error
  use options, only:                 check_null_error
  use parallel_linear_algebra, only: pdot_product, pmaxval, pnorm2
  use precision, only:               dp
  use woodbury_matrices, only:       dimCk, rpk, uC, uL

  implicit none

  real(kind=dp), allocatable, dimension(:, :) :: Ax, y
  real(kind=dp), allocatable, dimension(:)    :: Sx

  character(len=64)                           :: caststr

  if ( check_null_error ) then

     allocate( Ax(1:rpk, 1), y(1:rpk, 1:nsuby) )
     allocate( Sx(1:dimCk) )

     call notify( 'Computing error in null space computation (slow). ' )

     ! Compute the Poisson kernel error.
     call apply_smpm_poisson_transpose_parallel( Ax, uL )

     ! Compute the Schur kernel error.
     !
     ! NOTE: sigma is used within apply_poisson_capacitance_shift() and needs to be
     !       set rather than supplied as a value.
     sigma = 0.0_dp
     call apply_poisson_capacitance_shift( Sx, uC )

     l2_poisson_kernel_error = pnorm2( reshape( Ax, (/ rpk * 1 /) ) )
     l2_poisson_schur_kernel_error = pnorm2( Sx )

     write( caststr, '(D17.10)' ) l2_poisson_kernel_error
     call notify( '   Error in null space computation: ' // caststr )

     write( caststr, '(D17.10)' ) l2_poisson_schur_kernel_error
     call notify( '   Error in capacitance null space: ' // caststr )
     call notify( ' ' )

  else
     l2_poisson_kernel_error       = 0.0_dp
     l2_poisson_schur_kernel_error = 0.0_dp
  end if

end subroutine compute_null_space_error
