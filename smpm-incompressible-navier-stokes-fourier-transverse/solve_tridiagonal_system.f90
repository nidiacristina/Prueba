subroutine solve_tridiagonal_system( n, M, b )

   use precision, only: dp

   implicit none

   integer, intent(in)                            :: n
   real(kind=dp), dimension(1:n, 1:n), intent(in) :: M
   real(kind=dp), dimension(1:n), intent(inout)   :: b

   real(kind=dp)                                  :: q
   integer                                        :: ii
   real(kind=dp), dimension(1:n)                  :: Aii

   ! Application of the Thomas algorithm for tri-diagonal matrices.
   Aii = 0.0_dp

   ! Forward sweep.
   Aii(1) = M(1, 1)
   do ii = 2, n
      q       = M(ii, ii - 1) / Aii(ii - 1)
      Aii(ii) = M(ii, ii) - q * M(ii - 1, ii)
      b(ii)   = b(ii) - q * b(ii - 1)
   enddo

   ! Abbreviated back-substitution.
   b(n) = b(n) / Aii(n)
   do ii = n - 1, 1 , -1
      b(ii) = (b(ii) - M(ii, ii + 1) * b(ii + 1)) / Aii(ii)
   enddo

end subroutine solve_tridiagonal_system
