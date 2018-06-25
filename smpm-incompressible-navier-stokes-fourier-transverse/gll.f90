! Routines to calculate the values of the Gauss-Lobatto-Legendre points within
! the interval [-1, 1].  This code is almost identical to the code presented
! in Appendix C of the book "Spectral Methods in Fluid Dynamics" Canuto et.al
! 1982.
!
! NOTE: 'n' is the number of points in each direction!!!

subroutine jacobl( n, alpha, beta, xcol, ndim )
! Computes the gauss-lobatto collocation points for Jacobi polynomials
!
!   n:              Degree of approximation (order of polynomials)
!   alpha:          Parameter in Jacobi weight
!   beta:           Parameter in Jacobi weight
!   xcol:           Output array with the collocation points
!   ndim:           Dimension of the array "xcol"
!
!  for Chebyshev-Gauss-Lobatto points use alpha=-0.5 and beta=-0.5
!  for Legendre-Gauss-Lobatto points use           0             0

  use precision, only: dp

  implicit none

  ! Subtroutine parameters.
  integer, intent(in)             :: n             ! Degree of approximation for the polynomial.
  real(kind=dp), intent(in)       :: alpha         ! Parameter in the Jacobi weight.
  real(kind=dp), intent(in)       :: beta          ! Parameter in the Jacobi weight.
  real(kind=dp), intent(out)      :: xcol          ! Array of collocation points.
  integer, intent(in)             :: ndim          ! Number of collocation points to compute.

  dimension xjac(1000), xcol(0:ndim)

  ! Local variables.
  integer                         :: np, npp       ! Order of polynomial plus one and two, respectively.
  integer                         :: nh            ! Highest degree to compute the collocation points.
                                                   ! Approximately half of the degree provided by the
                                                   ! caller.

  real(kind=dp)                   :: pnp1p, pdnp1p, pnp, pdnp, pnm1p
  real(kind=dp)                   :: pnp1m, pdnp1m, pnm, pdnm, pnm1m
  real(kind=dp)                   :: pnp1, pdnp1, pn, pdn, pnm1, pdnm1

  real(kind=dp)                   :: a, b
  real(kind=dp)                   :: dth           ! Pi times twice the degree of the polynomial.
  real(kind=dp)                   :: cd, sd        ! Cosine and sine of twice dth, respectively.
  real(kind=dp)                   :: cs, ss        ! Cosine and sine of dth, respectively.
  real(kind=dp)                   :: cssave        ! Temporary storage for cs.

  real(kind=dp)                   :: delx
  real(kind=dp)                   :: det
  real(kind=dp)                   :: recsum
  real(kind=dp)                   :: pder
  real(kind=dp)                   :: poly
  real(kind=dp)                   :: rm, rp
  real(kind=dp)                   :: x             ! Point to evaluate the collocation point?

  real(kind=dp)                   :: xjac

  integer                         :: i, j, k       ! Loop indices.
  integer                         :: jm            ! j minus one.
  integer                         :: kk

  real(kind=dp)                   :: pi

  ! Constants.
  integer, parameter              :: kstop = 100       ! Maximum iterations.
  real(kind=dp), parameter        :: eps = 1.0e-15_dp  ! Tolerance.

  ! Get pi.
  pi = 2.0_dp * acos( 0.0_dp )

  np  = n + 1

  call jacobf( np, pnp1p, pdnp1p, pnp, pdnp, pnm1p, pdnm1,  1.0_dp, alpha, beta, 1.0_dp + alpha )
  call jacobf( np, pnp1m, pdnp1m, pnm, pdnm, pnm1m, pdnm1, -1.0_dp, alpha, beta, 1.0_dp + alpha )

  det     =  pnp * pnm1m - pnm * pnm1p
  rp      = -pnp1p
  rm      = -pnp1m
  a       = (rp * pnm1m - rm * pnm1p) / det
  b       = (rm * pnp   - rp * pnm) / det
  xjac(1) = 1.0_dp
  nh      = (n + 1) / 2

  dth     = pi / (2 * n  + 1)
  cd      = cos( 2.0_dp * dth )
  sd      = sin( 2.0_dp * dth )
  cs      = cos( dth )
  ss      = sin( dth )

  do j = 2,nh
     x = cs

     do k = 1,kstop
        call jacobf( np, pnp1, pdnp1, pn, pdn, pnm1, pdnm1, x, alpha, beta, 1.0_dp + alpha )

        poly   = pnp1  + (a * pn)  + (b * pnm1)
        pder   = pdnp1 + (a * pdn) + (b * pdnm1)

        recsum = 0.0_dp
        jm     = j - 1

        do i = 1,jm
           recsum = recsum + 1.0_dp / (x - xjac(i))
        enddo

        delx = -poly / (pder - recsum * poly)
        x    = x + delx

        if( abs( delx ) .lt. eps ) exit

     enddo

     xjac(j) = x
     cssave  = (cs * cd) - (ss * sd)
     ss      = (cs * sd) + (ss * cd)
     cs      = cssave
  enddo

  xjac(np) = -1.0_dp
  npp      = n + 2

  do i = 2,nh
     xjac(npp-i) = -xjac(i)
  enddo

  ! Copy the collocation points to output array.
  if( n .ne. 2 * (n / 2) ) then
     do k = 0,n
        kk      = n - k + 1
        xcol(k) = xjac(kk)
     enddo
  else
     xjac(nh+1) = 0.0_dp

     do k = 0,n
        kk      = n - k + 1
        xcol(k) = xjac(kk)
     enddo
  endif

  return

end subroutine jacobl

subroutine jacobf( n, poly, pder, polym1, pderm1, polym2, pderm2, x, alpha, beta, rv )
! Computes the Jacobi polynomial (in this case Legendre) and its derivative of
! degree n at x
!
!   n ->    Polynomial degree
!   poly->  Magnitude of Jacobi (Legendre) polynomial evaluated in "x "
!   pder->  Magnitude of the derivative of Jacobi polynomial in "x"
!   polym1->Polynomial of degree n-1
!   pderm1->Derivative of the polynomial of degree n-1
!   polym2->Polynomial of degree n-2
!   pderm2->Derivative of the polynomial of degree n-2
!   x ->    Point in the space
!   alpha-> Jacobian weight
!   beta->  Jacobian weight
!   rv->    XXX: unclear what this parameter does (GNT)
!
!   This subroutine is an unabridged version of the original one presented in
!   Canuto et. al. (1982)

  use precision, only: dp

  implicit none

  ! Subroutine parameters.
  integer, intent(in)          :: n
  real(kind=dp), intent(inout) :: poly
  real(kind=dp), intent(inout) :: pder
  real(kind=dp), intent(out)   :: polym1
  real(kind=dp), intent(out)   :: pderm1
  real(kind=dp), intent(out)   :: polym2
  real(kind=dp), intent(out)   :: pderm2
  real(kind=dp), intent(in)    :: x
  real(kind=dp), intent(in)    :: alpha
  real(kind=dp), intent(in)    :: beta
  real(kind=dp), intent(in)    :: rv

  ! Local variables.
  real(kind=dp)                :: apb
  real(kind=dp)                :: a1, a2, b3, a3, a4
  real(kind=dp)                :: polyn, pdern
  real(kind=dp)                :: psave, pdsave, polylst, pderlst

  integer                      :: k

  apb  = alpha + beta
  poly = 1.0_dp
  pder = 0.0_dp

  if( n .eq. 0 ) return

  polylst = poly
  pderlst = pder
  poly    = rv * x
  pder    = rv

  if( n .eq. 1 ) return

  ! NOTE: These aren't strictly necessary, though gfortran complains that
  !       these variables may not be initialized before use.  Since n has to
  !       be larger than 1 at this point, they will always get assigned
  !       through at least the first trip through the loop.
  psave   = 0.0_dp
  pdsave  = 0.0_dp

  do k = 2,n
     a1      = 2.0_dp * k * (k + apb) * (2.0_dp * k + apb - 2.0_dp)
     a2      = (2.0_dp * k + apb - 1.) * (alpha**2 - beta**2)
     b3      = (2.0_dp * k + apb - 2.0_dp)
     a3      = b3 * (b3 + 1.) * (b3 + 2.0_dp)
     a4      = 2.0_dp * (K + alpha - 1.) * (k + beta - 1.) * (2.0_dp * k + apb)

     polyn   = ((a2 + a3 * x) * poly - a4 * polylst) / a1
     pdern   = ((a2 + a3 * x) * pder - a4 * pderlst + a3 * poly) / a1

     psave   = polylst
     pdsave  = pderlst
     polylst = poly
     poly    = polyn
     pderlst = pder
     pder    = pdern
  enddo

  polym1 = polylst
  pderm1 = pderlst
  polym2 = psave
  pderm2 = pdsave

  return

end subroutine jacobf
