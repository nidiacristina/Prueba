subroutine quad( n, x, w, ndim )

  use precision, only: dp

  implicit none

  ! Constants.
  integer, parameter                            :: nn = 1000
  real(kind=dp), parameter                      :: small = 1.0e-30_dp

  ! Subroutine parameters.
  integer, intent(in)                           :: n
  real(kind=dp), dimension(0:ndim), intent(in)  :: x
  real(kind=dp), dimension(0:ndim), intent(out) :: w
  integer, intent(in)                           :: ndim

  ! Local variables.
  real(kind=dp), dimension(0:nn)                :: al1
  integer                                       :: k   ! loop index

  ! Also store values of 0 to ndim Legendre polynomials at all collocation
  ! points.  Determine the Gauss Quadrature weighting factors.

  do k = 0, n
     call legen( al1, n, x(k), nn )

     ! Calculate weighting factor.
     w(k) = 2.0_dp / (n * (n + 1) * al1(n) * al1(n) + small)
  enddo

  return

end subroutine quad


subroutine legen( al, n, xc, ndim )
! Calculates values of all Legendre polynomials (and immediately highest one
! in hierarchy) at a given collocation point xc.

  use precision, only: dp

  implicit none

  ! Subroutine parameters.
  real(kind=dp), dimension(0:ndim), intent(out) :: al
  integer, intent(in)                           :: n
  real(kind=dp), intent(in)                     :: xc
  integer, intent(in)                           :: ndim

  ! Local variables.
  integer                                       :: k      ! loop index
  integer                                       :: kp, km ! k plus/minus one, respectively

  al(0) = 1.0_dp
  al(1) = xc

  do k = 1, n
     kp     = k + 1
     km     = k - 1
     al(kp) = (2.0_dp * k + 1.0_dp) * xc * al(k) / kp - k * al(km) / kp
  enddo

  return

end subroutine legen
