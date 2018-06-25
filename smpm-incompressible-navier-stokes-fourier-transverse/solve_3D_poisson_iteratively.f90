subroutine solve_3D_poisson_iteratively( x , gmres_poisson_iterations, l2_poisson_schur_error, cmplx_part  )
!
! Uses a distributed matvec for the Poisson
! capacitance matrix to solve the capacitance equation,
! then uses the Woodbury matrix identity to extract the
! solution of the full Poisson problem.

  use constants, only:               nky, nsg, nsubx, rank, root
  use options, only:                 check_numerical_error, gmres_maxit_poisson, &
                                     gmres_tol_poisson, use_capacitance_preconditioner, &
                                     use_deflation
  use woodbury_matrices, only:       k, uC, uL, S_coarse,  &
                                     rpk, dimCk, uC_right, uC_coarse, &
                                     xC_previous_real, xC_previous_imag
  use parallel_linear_algebra,only:  pdot_product, pnorm2
  use precision, only:               dp
  use timestepping, only:            logflag

  implicit none

  real(kind=dp), dimension(1:rpk, 1:nky), intent(inout) :: x
  integer, intent(out)                                  :: gmres_poisson_iterations
  real(kind=dp), intent(out)                            :: l2_poisson_schur_error

  real(kind=dp), dimension(1:dimCk, 1:nky)              :: Ax, xC
  real(kind=dp), dimension(1:rpk, 1:nky)                :: Aib, Aib_x, Aib_z
  real(kind=dp), dimension(1:dimCk, 1:nky)              :: bC
  real(kind=dp)                                         :: gmres_tol
  character(len=32)                                     :: caststr

  integer                                               :: gmres_maxit, gmres_restart, ii
  real(kind=dp), dimension(1:rpk)                       :: v0
  real(kind=dp), dimension(1:rpk, 1:nky)                :: ExC
  real(kind=dp), dimension(1:rpk, 1:nky)                :: AiExC
  real(kind=dp), dimension(1:rpk, 1:nky)                :: b
  real(kind=dp)                                         :: denominator

  real(kind=dp), dimension(1:dimCk, 1:nky)              :: Mx, term1, term2
  real(kind=dp), dimension(1:nsubx - 1, 1:nky)          :: ZTx
  real(kind=dp), dimension(1:dimCk, 1:nky)              :: PbC

  character, intent(in)                                 :: cmplx_part ! Part of RHS to be solved
                                                                      ! R (real) ,I (imag) or N (none)

  ! Declare the function used with GMRES to solve the capacitance problem.
  external                                              :: apply_3D_poisson_schur
  external                                              :: apply_3D_poisson_schur_sparse
  external                                              :: apply_3D_preconditioned_schur
  external                                              :: apply_3D_deflated_schur

  ! Set some constants.
  gmres_tol     = gmres_tol_poisson
  gmres_maxit   = min( gmres_maxit_poisson, nky * (k-1) )
  gmres_restart = gmres_maxit

  ! Regularize the 0-wavenumber right hand side.
  x(:, 1) = x(:, 1) - pdot_product( uL, x( :, 1)) * uL
  b = x

  ! Divide A into the right hand side.
  Aib = 0.0_dp
  call solve_3D_A( x, Aib, 'N' )

  ! Multiply B into A\b to obtain the right-hand-side for the Schur problem.
  bC = 0.0_dp
  do ii = 1, nky
     call compute_gradient( Aib(:, ii), Aib_x(:, ii), Aib_z(:, ii) )
  enddo
  call apply_pB_poisson( Aib, bC, Aib_x, Aib_z, nky )

  ! Regularize the 0-wavenumber Schur problem.
  bC( :, 1 ) = bC( :, 1 ) - uC * pdot_product( uC, bC( :,  1 ) ) 

  ! Set the initial guess.
  if ( cmplx_part == 'R' ) then
     xC = xC_previous_real

  elseif ( cmplx_part == 'I' ) then
     xC = xC_previous_imag

  elseif ( cmplx_part == 'N' ) then
     xC = 0.0_dp

  endif

  ! Solve the deflated and preconditioned GMRES problem.
  call notify_cond( logflag, '            Solving the GMRES problem.' )
  if ( use_deflation ) then
     PbC = 0.0_dp
     call apply_3D_P( PbC, bC )
     gmres_restart = gmres_maxit
     call compute_gmres_householder( xC, PbC, dimCk * nky, gmres_tol, gmres_maxit, gmres_restart, &
                                     apply_3D_deflated_schur )

     ! Store the not-preconditioned solution for the next solve
     if ( cmplx_part == 'R' ) then
        xC_previous_real = xC

     elseif ( cmplx_part == 'I' ) then
        xC_previous_imag = xC

     endif

     Mx = 0.0_dp
     call apply_3D_schur_preconditioner( Mx, xC )
     xC = Mx

  elseif ( use_capacitance_preconditioner ) then     !Solve the undeflated preconditioned GMRES Problem
     call compute_gmres_householder( xC, bC, dimCk * nky, gmres_tol, gmres_maxit, &
                                     gmres_restart, apply_3D_preconditioned_schur )

     ! Store the not-preconditioned solution for the next solve
     if ( cmplx_part == 'R' ) then
        xC_previous_real = xC

     elseif ( cmplx_part == 'I' ) then
        xC_previous_imag = xC

     endif

     Mx = 0.0_dp
     call apply_3D_schur_preconditioner( Mx, xC )
     xC = Mx

  else                                               !Solve the undeflated and notpreconditioned GMRES problem
     call compute_gmres_householder( xC, bC, dimCk * nky, gmres_tol, gmres_maxit, gmres_restart, &
                                     apply_3D_poisson_schur_sparse )

     ! Store the solution for the next solve
     if ( cmplx_part == 'R' ) then
        xC_previous_real = xC

     elseif ( cmplx_part == 'I' ) then
        xC_previous_imag = xC

     endif

  endif !End of GMRES Problem

  ! Notify the user of some solver statistics.
  if ( rank == root ) then
     write( caststr, '(I10) ' ) gmres_maxit
     call notify_cond( logflag, '             Schur GMRES iterations           : ' // caststr )
     write( caststr, '(D17.10)' ) gmres_tol
     call notify_cond( logflag, '             Schur GMRES residual             : ' // caststr )
  endif

  gmres_poisson_iterations = gmres_maxit
  l2_poisson_schur_error   = gmres_tol

  ! Check to see if the GMRES iteration actually converged.
  if ( check_numerical_error ) then

     Ax = 0.0_dp
     denominator = 0.0_dp

     ! Note that we use the full operator (i.e. apply_3D_poisson_schur_sparse to check the error)
     if ( use_deflation) then
        call apply_3D_poisson_schur_sparse( Ax, xC )
        call apply_3D_P( Mx, Ax )

        denominator = pnorm2( reshape( PbC, (/ nky * dimCk /) ) )
        if ( denominator .ne. 0.0_dp ) then
           write( caststr, '(D17.10)' ) pnorm2( reshape( Mx - PbC, (/ nky * dimCk /) ) ) & 
                                      / denominator
        else
           write( caststr, '(D17.10)' ) 0.0_dp
        endif

     else ! If no deflation was used

        call apply_3D_poisson_schur_sparse( Ax, xC )

        denominator = pnorm2( reshape( bC, (/ nky * dimCk /) ) )
        if ( denominator .ne. 0.0_dp ) then
           write( caststr, '(D17.10)' ) pnorm2( reshape( Ax - bC, (/ nky * dimCk /) ) ) & 
                                      / denominator
        else
           write( caststr, '(D17.10)' ) 0.0_dp
        endif

     endif ! End deflation check
 
     if ( rank == root ) then
         call notify_cond( logflag, '            Relative error in GMRES           : ' // caststr )
     endif

  endif

  ! Solve the coarse solve needed to assemble the solution.
  if ( use_deflation ) then
     Ztx = 0.0_dp
     call apply_Z_transpose( bC, ZTx, nky )
     ZTx( :, 1 ) = ZTx( :, 1 ) - uC_coarse * dot_product( ZTx( :, 1 ), uC_coarse )

     do ii = 1, nky
        call solve_tridiagonal_system( nsubx - 1, S_coarse( :, :, ii ), ZTx(:, ii ) )
     enddo

     ! Inflate the coarse solution onto the full grid
     term1 = 0.0_dp
     call apply_Z( ZTx, term1, nky )

     ! Post-multiply the GMRES solution by the Q-projection to get the second term.
     term2 = 0.0_dp
     call apply_3D_Q( term2, xC )

     ! Assemble the Schur solution.
     xC = term1 + term2

     ! Regularize the schur solution.
     ! GAR: Do we need it? Its not present in the master branch.
     xC(:, 1) = xC(:, 1) - uC_right * pdot_product( uC_right, xC(:,1) )

  endif

  ! From the capacitance solution back out the solution to the Poisson problem.
  ExC   = 0.0_dp
  AiExC = 0.0_dp
  call apply_extension( xC, ExC, nky )
  call solve_3D_A( ExC, AiExC, 'N' )

  ! Think about potentially doing this as A\(b - ExC) instead of A\b - A\ExC.
  x = Aib - AiExC

  ! Regularize.
  v0 = 1.0_dp / sqrt( real( nsg , kind=dp ) )
  x(:, 1) = x(:, 1) - v0 * pdot_product( v0, x(:,1) ) 


end subroutine solve_3D_poisson_iteratively
                                                                 

