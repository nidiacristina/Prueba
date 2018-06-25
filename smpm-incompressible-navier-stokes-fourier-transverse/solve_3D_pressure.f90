subroutine solve_3D_pressure( p )

! Wrapper that solves the 3D poisson equation in SMPM Fourier Transverse
!
! 8 August 2017
! Gustavo Rivera

   ! Declare modules

   use constants, only:               nky, nsuby
   use errors, only:                  gmres_poisson_iterations_real, gmres_poisson_iterations_imag, &
                                      l2_poisson_error, l2_poisson_schur_error
   use options, only:                 check_numerical_error
   use precision, only:               dp
   use parallel_linear_algebra, only: pdot_product, pmaxval, pnorm2
   use timestepping, only:            logflag
   use woodbury_matrices, only:       rpk, uL

   implicit none

   ! Set the input/output variables.
   real(kind=dp), dimension(1:rpk, 1:nsuby), intent(inout) :: p
   
   ! Set the error variables
   real(kind=dp)                                           :: l2_poisson_schur_error_real
   real(kind=dp)                                           :: l2_poisson_schur_error_imag

   ! Buffers associated with pressure
   complex(kind=dp), dimension(1:rpk,1:nky)                :: Fp, Fp_rhs
   real(kind=dp), dimension(1:rpk, 1:nsuby)                :: p_rhs
   real(kind=dp), dimension(1:rpk, 1:nsuby)                :: p_rhs_real, p_rhs_imag
   real(kind=dp), dimension(1:rpk, 1:nky)                  :: Fp_real, Fp_imag

   ! Buffers associated with error
   real(kind=dp), allocatable, dimension(:,:)              :: Ax
   real(kind=dp), allocatable, dimension(:,:)              :: Ax_imag, Ax_real
   real(kind=dp), allocatable, dimension(:)                :: v0
   real(kind=dp)                                           :: denominator, thiserror

   character(len=32)                                       :: caststr
   integer                                                 :: ii

   ! Step 2.2.1: Apply the boundary conditions to the RHS.
   call apply_pressure_bc( p )
   if ( check_numerical_error ) p_rhs = p

   ! Step 2.3.2: Fourier transform the RHS.
   Fp_real = 0.0_dp
   Fp_imag = 0.0_dp
   Fp      = cmplx( 0.0_dp, 0.0_dp, kind=dp )
   call apply_fft( p, Fp )

   ! Step 2.3.3: Separate the complex RHS into real and imaginary components, on real buffers.
   do ii = 1, nky
      Fp_real( :,ii ) = real(  Fp( :,ii ), kind=dp )
      Fp_imag( :,ii ) = aimag( Fp( :,ii ) )
   enddo

   ! Step 2.3.4: Solve the real part.
   call notify_cond( logflag, ' ' )
   call notify_cond( logflag, '         Solving the Real Part ')
   call notify_cond( logflag, ' ' )

   call solve_3D_poisson_iteratively( Fp_real, gmres_poisson_iterations_real, &
                                      l2_poisson_schur_error_real, 'R' )

   ! Step 2.3.5: Solve the imaginary part.
   !             Check if there is an transverse domain first.
   if ( nsuby > 1) then
      call notify_cond( logflag, ' ' )
      call notify_cond( logflag, '         Solving the Imaginary Part ')
      call notify_cond( logflag, ' ' )

      call solve_3D_poisson_iteratively( Fp_imag, gmres_poisson_iterations_imag, &
                                         l2_poisson_schur_error_imag, 'I' )

   else ! Skip the imag solve associated with the transverse direction
      call notify_cond( logflag, ' ' )
      call notify_cond( logflag, '         Skipping Imaginary Part ')
      call notify_cond( logflag, ' ' )

      gmres_poisson_iterations_imag = 0
      l2_poisson_schur_error_imag   = 0.0_dp

   endif

   ! Store error
   l2_poisson_schur_error = cmplx( l2_poisson_schur_error_real, &
                                   l2_poisson_schur_error_imag, &
                                   kind=dp )

   ! Step 2.3.6: Assemble the solution in a complex buffer.
   do ii = 1,nky
      Fp(:,ii) = cmplx( Fp_real(:,ii), Fp_imag(:,ii), kind=dp)
   enddo

   ! Step 2.3.7: Transform the solution back to physical space.
   call apply_ifft( Fp, p )

   ! Step 2.3.8: Compute the numerical errors in the Poisson step.
   if ( check_numerical_error ) then

      ! Step 2.3.8.1: Allocate error buffers.
      allocate(v0(1:rpk))
      allocate(Ax(1:rpk,1:nsuby))
      allocate(Ax_real(1:rpk,1:nky))
      allocate(Ax_imag(1:rpk,1:nky))

      ! Step 2.3.8.2: Transform the stored RHS.
      call apply_fft(p_rhs, Fp_rhs)
      do ii = 1, nky
         p_rhs_real( :,ii ) = real(  Fp_rhs( :,ii ), kind=dp )
         p_rhs_imag( :,ii ) = aimag( Fp_rhs( :,ii ) )
      enddo

      ! Step 2.3.8.3: Regularize the stored RHS.
      p_rhs_real(:,1) = p_rhs_real(:, 1) - uL * pdot_product( uL, p_rhs_real(:, 1) )
      p_rhs_imag(:,1) = p_rhs_real(:, 1) - uL * pdot_product( uL, p_rhs_imag(:, 1) )

      ! Step 2.3.8.4: Rebuild the stored RHS.
      do ii = 1,nky
         Fp_rhs(:,ii) = cmplx( p_rhs_real(:,ii), p_rhs_imag(:,ii), kind=dp )
      enddo

      ! Step 2.3.8.5: Inverse transform the stored RHS.
      call apply_ifft(Fp_rhs, p_rhs)

      ! Step 2.3.8.6: Transform the computed pressure.
      call apply_fft(p, Fp)
      do ii = 1,nky
         Fp_real(:,ii) = real(  Fp(:,ii), kind=dp )
         Fp_imag(:,ii) = aimag( Fp (:,ii) )
      enddo

      ! Step 2.3.8.7: Apply undeflated operator on computed pressure.
      call apply_3D_poisson( Ax_real , Fp_real )
      call apply_3D_poisson( Ax_imag , Fp_imag )

      ! Step 2.3.8.8: Reassemble RHS from computed pressure.
      do ii = 1,nky
         Fp(:,ii) = cmplx( Ax_real(:, ii), Ax_imag(:,ii), kind=dp )
      enddo
      call apply_ifft(Fp, Ax)

      ! Step 2.3.8.9: Calculate L2 and Linf error || Ax - b || with L2 and Linf norms.
      denominator      = pnorm2( reshape( p_rhs     , (/ nsuby * rpk /) ) )
      if ( denominator .ne. 0.0_dp) then
         l2_poisson_error = pnorm2( reshape( Ax - p_rhs, (/ nsuby * rpk /) ) ) / &
                            denominator
      else
         l2_poisson_error = 0.0_dp
      endif
      write( caststr, '(D17.10)' ) l2_poisson_error
      call notify_cond( logflag, ' ' )
      call notify_cond( logflag, '         Rel. L2 error in pressure solve      : ' // caststr )

      thiserror = pmaxval( reshape( abs( Ax - p_rhs ), (/ nsuby * rpk /) ) )
      write( caststr, '(D17.10)' ) thiserror
      call notify_cond( logflag, '         Linf error in pressure solve         : ' // caststr )

      ! Step 2.3.8.10: Deallocate arrays.
      deallocate(v0)
      deallocate(Ax)
      deallocate(Ax_real)
      deallocate(Ax_imag)

   endif

end subroutine solve_3D_pressure
