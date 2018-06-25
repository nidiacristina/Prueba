subroutine solve_3D_A( q, Aq, trans )
!
! For an rpk x nky array representing a local grid function across
! all wavenumbers, divide this grid function by A, shifted appropriately
! for all wavenumbers.


  use constants,         only: nky
  use precision,         only: dp
  use woodbury_matrices, only: T_poisson, dimA, numA_per_rank, rpk, U_poisson
  use transverse,        only: qy

  implicit none

  real(kind=dp), dimension(1:rpk, 1:nky), intent(in)  ::  q
  real(kind=dp), dimension(1:rpk, 1:nky), intent(out) :: Aq
  character, intent(in)                               :: trans ! Either 'T' or 'N' as used in BLAS.

  ! Internal variables
  integer                                             :: iistart, iiend, info, ii, jj, kk

  complex(kind=dp), dimension( 1:rpk, 1:nky )         :: zq, zAq


  ! Multiply the transpose of U into the input.
  do ii = 1, numA_per_rank

     iistart = (ii - 1) * dimA + 1
     iiend   = (ii - 1) * dimA + dimA
     zAq( iistart:iiend, : ) = cmplx( q( iistart:iiend, : ), kind=dp )
     zq( iistart:iiend, : )  = cmplx( q( iistart:iiend, : ), kind=dp )

     ! Multiply transpose(U) into the RHS.
     call ZGEMM( 'C', 'N', dimA, nky, dimA, cmplx( 1.0_dp, 0.0_dp, kind=dp ), U_poisson( :, :, ii ), dimA, &
                 zq( iistart, 1 ), rpk, cmplx( 0.0_dp, 0.0_dp, kind=dp ), zAq( iistart, 1 ), rpk )

  enddo

  ! Solve either T or its transpose with a shift representing each wavenumber.
  do ii = 1, numA_per_rank

     iistart = (ii - 1) * dimA + 1
     iiend   = (ii - 1) * dimA + dimA

     ! Loop over the wavenumbers, solving each.
     do kk = 1, nky

        ! Shift the triangular system.
        do jj = 1, dimA
            T_poisson( jj, jj, ii ) = T_poisson( jj, jj, ii ) - qy(kk)**2
        enddo

        ! Solve the triangular system.
        if ( trans == 'N' ) then
           call ZTRTRS( 'U', 'N', 'N', dimA, 1, T_poisson( :, :, ii ), dimA, &
                        zAq( iistart, kk ), rpk, info )
        else
           call ZTRTRS( 'U', 'C', 'N', dimA, 1, T_poisson( :, :, ii ), dimA, &
                        zAq( iistart, kk ), rpk, info )
        endif

        ! Un-shift the triangular system.
        do jj = 1, dimA
            T_poisson( jj, jj, ii ) = T_poisson( jj, jj, ii ) + qy(kk)**2
        enddo

     enddo

  enddo

  ! Multiply U into the result.
  do ii = 1, numA_per_rank

     iistart = (ii - 1) * dimA + 1
     iiend   = (ii - 1) * dimA + dimA

     ! Multiply U into this solution.
     call ZGEMM( 'N', 'N', dimA, nky, dimA, cmplx( 1.0_dp, 0.0_dp, kind=dp ), U_poisson( :, :, ii ), dimA, &
                 zAq( iistart, 1 ), rpk, cmplx( 0.0_dp, 0.0_dp, kind=dp ), zq( iistart, 1 ), rpk )

     ! Store the output.
     Aq( iistart:iiend, : ) = real( zq( iistart:iiend, :  ), kind=dp )

  enddo

end subroutine solve_3D_A
                                                         

