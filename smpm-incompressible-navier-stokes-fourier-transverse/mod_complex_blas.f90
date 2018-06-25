module complex_blas
! Contains wrappers for BLAS and LAPACK routines that allows for real
! matrix/complex vector algebra.

  implicit none

contains

  ! Wrap DGEMV.
  !
  ! XXX: Notice that I always assume a unit stride in X/Y.
  !      THIS WILL BREAK IF YOU PASS IN INCX/INCY != 1.
  subroutine DGEMV_complex( trans, m, n, alpha, A, LDA, X, incx, beta, Y, incy )

     use precision, only: dp

     implicit none

     character, intent(in)           :: trans
     integer, intent(in)             :: m
     integer, intent(in)             :: n
     real(kind=dp), intent(in)       :: alpha
     real(kind=dp), intent(in)       :: A(LDA,*)
     integer, intent(in)             :: LDA
     complex(kind=dp), intent(in)    :: X(*)
     integer, intent(in)             :: incx
     real(kind=dp), intent(in)       :: beta
     complex(kind=dp), intent(inout) :: Y(*)
     integer, intent(in)             :: incy

     real(kind=dp), dimension(1:n)   :: real_X, imag_X, real_Y, imag_Y

     real_X =  real( X(1:n:incx), kind=dp )
     imag_X = aimag( X(1:n:incx) )
     real_Y =  real( Y(1:n:incy), kind=dp )
     imag_Y = aimag( Y(1:n:incy) )

     ! Call DGEMV() twice, once for the real part, once for the imaginary
     ! part.
     call DGEMV( trans, m, n, alpha, A, LDA, real_X, incx, beta, real_Y, incy )
     call DGEMV( trans, m, n, alpha, A, LDA, imag_X, incx, beta, imag_Y, incy )

     ! Re-assemble the solution.
     Y(1:n:incy) = cmplx( real_Y, imag_Y, kind=dp )

  end subroutine DGEMV_complex

  ! Wrap DGETRS.
  subroutine DGETRS_complex( trans, n, nrhs, A, LDA, IPIV, B, LDB, info )

     use precision, only: dp

     implicit none

     character, intent(in)                 :: trans
     integer, intent(in)                   :: n
     integer, intent(in)                   :: nrhs
     real(kind=dp), intent(in)             :: A(LDA,*)
     integer, intent(in)                   :: LDA
     integer, dimension(1:n), intent(in)   :: IPIV
     complex(kind=dp), intent(inout)       :: B(LDB, *)
     integer, intent(in)                   :: LDB
     integer, intent(out)                  :: info

     real(kind=dp), dimension(1:n, 1:nrhs) :: real_buffer, imag_buffer

     ! Call DGETRS() twice, once for the real part and once for the imaginary
     ! part.
     real_buffer =  real( B(1:n, 1:nrhs), kind=dp )
     imag_buffer = aimag( B(1:n, 1:nrhs) )
     call DGETRS( trans, n, nrhs, A, LDA, IPIV, real_buffer, LDB, info )
     call DGETRS( trans, n, nrhs, A, LDA, IPIV, imag_buffer, LDB, info )

     ! Reassemble the solution.
     B(1:n, 1:nrhs) = cmplx( real_buffer, imag_buffer, kind=dp )

  end subroutine DGETRS_complex

end module complex_blas
