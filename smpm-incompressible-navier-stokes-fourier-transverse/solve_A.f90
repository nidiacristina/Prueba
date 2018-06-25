subroutine solve_A( q, Aq, Nrhs, trans, diagshift )
! Solves (A - diagshift)\q or (A - diagshift)'\q and returns it in Aq, where A
! is the block-diagonal part of the Poisson matrix.

  use precision, only:         dp
  use woodbury_matrices, only: dimA, numA_per_rank, rpk, T_poisson, U_poisson

  implicit none

  real(kind=dp), dimension(1:rpk, 1:Nrhs), intent(in)  ::  q
  real(kind=dp), dimension(1:rpk, 1:Nrhs), intent(out) :: Aq
  integer, intent(in)                                  :: Nrhs
  character, intent(in)                                :: trans ! Either 'T' or 'N' as used in BLAS.
  complex(kind=dp), intent(in)                         :: diagshift

  ! Internal variables.
  integer                                              :: iistart, iiend, info, ii, jj
  complex(kind=dp), dimension(1:rpk, 1:Nrhs)           :: zq, zAq

  ! Constants
  complex(kind=dp), parameter                         :: complex_one = (1.0_dp, 0.0_dp)
  complex(kind=dp), parameter                         :: complex_zero = (0.0_dp, 0.0_dp)

  ! Divide A into each vertical column of subdomains.
  do ii = 1, numA_per_rank

     iistart               = (ii - 1) * dimA + 1
     iiend                 = (ii - 1) * dimA + dimA
     zAq(iistart:iiend, :) = cmplx( q(iistart:iiend, :), kind=dp )
     zq(iistart:iiend, :)  = cmplx( q(iistart:iiend, :), kind=dp )

     ! Multiply transpose(U) into the RHS.
     call ZGEMM( 'C', 'N', dimA, Nrhs, dimA, complex_one, U_poisson(:, :, ii), dimA, &
                 zq(iistart, 1), rpk, complex_zero, zAq(iistart, 1), rpk )

     ! Shift the triangular system.
     do jj = 1, dimA
         T_poisson(jj, jj, ii) = T_poisson(jj, jj, ii) - diagshift
     enddo

     ! Solve the triangular system.
     !
     ! NOTE: This operates on subsets of rows across all columns.
     if ( trans == 'N' ) then
        call ZTRTRS( 'U', 'N', 'N', dimA, Nrhs, T_poisson(:, :, ii), dimA, &
                     zAq(iistart, 1), rpk, info )
     else
        call ZTRTRS( 'U', 'C', 'N', dimA, Nrhs, T_poisson(:, :, ii), dimA, &
                     zAq(iistart, 1), rpk, info )
     endif

     ! Un-shift the triangular system.
     do jj = 1, dimA
         T_poisson(jj, jj, ii) = T_poisson(jj, jj, ii) + diagshift
     enddo

     ! Multiply U into this solution.
     call ZGEMM( 'N', 'N', dimA, Nrhs, dimA, complex_one, U_poisson(:, :, ii), dimA, &
                 zAq(iistart, 1), dimA, complex_zero, zq(iistart, 1), rpk )

     ! Store the output.
     Aq(iistart:iiend, :) = real( zq(iistart:iiend, :), kind=dp )

  enddo

end subroutine solve_A
