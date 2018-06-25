subroutine setup_filter_xz( n, points, F, wg, pf )
! Generates the matrix for spectral filtering in physical space,
! meaning x-z. 

  use precision, only: dp

  implicit none

  integer, intent(in)                         :: n
  real(kind=dp), dimension(n), intent(in)     :: points
  real(kind=dp), dimension(n, n), intent(out) :: F
  real(kind=dp), dimension(n), intent(in)     :: wg
  real(kind=dp), intent(in)                   :: pf

  ! Local variables.
  real(kind=dp), dimension(n, n)              :: B, L, M, C, W
  real(kind=dp), dimension(n)                 :: P
  integer                                     :: kp, k, km, i, j
  real(kind=dp)                               :: alpha

  ! B matrix -> Legendre polynomial evaluated in GLL points.
  do j = 1, n
     P(1) = 1.0_dp
     P(2) = points(j)
     do k = 2, n-1
        kp    = k + 1
        km    = k - 1
        P(kp) = (((2.0_dp * (k - 1)) + 1.0_dp) * P(2) * P(k) / (k + 1 - 1)) - ((k - 1) * P(km) / (k + 1 - 1))
     enddo

     do i = 1, n
        B(j, i) = P(i)
     enddo
  enddo

  C = 0.0_dp
  W = 0.0_dp
  L = 0.0_dp

  do i = 1, n-1
     C(i, i) = (i - 1.0_dp) + 0.5_dp
     W(i, i) = wg(i)
  enddo

  C(n, n) = n / 2.0_dp
  W(n, n) = wg(n)

  M = matmul( C, transpose( B ) )
  M = matmul( M, W )

  ! Exponential filter.
  alpha = -log( 1.0e-16_dp )

  do k = 1, n
     L(k, k) = exp( -alpha * ((k - 1.0_dp) / (n - 1.0_dp))**pf )
  enddo

  F = matmul( B, L )
  F = matmul( F, M )

end subroutine setup_filter_xz
