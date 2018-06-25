subroutine derv( nterm, x, d, d2, d3, ndim )
! This code contains the routine employeed to calculate the first, second, and
! third derivative matrix at Gauss-Lobatto-Legendre grid points calculated in
! the routines jacobl and jacobf.
!
! The theory behind this code can be found in Costa & Don "On the computation
! of high order derivatives".  Applied Numerical Mathematics. 33 (2000) 151 -
! 159
!
!   List of Variables:
!   nterm ->    Polynomial degree
!   x ->        Array with Gauss-Lobatto-Legendre Points
!   d ->        First derivative matrix
!   d2 ->       Second derivetive matrix
!   d3 ->       Third derivative matrix
!   ndim ->     Number of grid points in each direction

  use precision, only: dp

  implicit none

  ! Constants.
  integer, parameter                                    :: nn = 1000

  ! Subtroutine parameters.
  integer, intent(in)                                   :: nterm
  real(kind=dp), dimension(0:ndim), intent(in)          :: x
  real(kind=dp), dimension(0:ndim, 0:ndim), intent(out) :: d
  real(kind=dp), dimension(0:ndim, 0:ndim), intent(out) :: d2
  real(kind=dp), dimension(0:ndim, 0:ndim), intent(out) :: d3
  integer, intent(in)                                   :: ndim

  ! Local variables.
  real(kind=dp), dimension(0:nn)                        :: c
  real(kind=dp)                                         :: prod, sum

  integer                                               :: j, k, l     ! Loop indices.
  real(kind=dp)                                         :: xj, xk, xl  ! Temporary variables holding
                                                                       ! the jth, kth and lth,
                                                                       ! elements of x.

  ! Introduce methodology by Costa and Don.  see Eqn (4)
  ! 1) Must first compute coefficients C_i
  do k = 0, nterm
     prod = 1.0_dp
     xk   = x(k)

     do l = 0,nterm
        xl = x(l)

        if (l .ne. k) then
           prod = prod * (xk - xl)
        else
           prod = prod
        endif
     enddo

     c(k) = prod
  enddo

  ! First Derivative:
  ! Calculate off diagonal elements of D^1.  see Eqn (6)
  do k = 0, nterm
     xk = x(k)

     do j = 0, nterm
        xj = x(j)

        if (k .ne. j) then
           d(k, j) = c(k) / (c(j) * (xk - xj))
        else
           d(k, j) = 0.0_dp
        endif
     enddo
  enddo

   ! Calculate diagonal elements now.  Diagonal element is negative row sum of
   ! off diagonal elements (See Costa & Don) Eqn(9)
   do k = 0, nterm
      sum = 0.0_dp

      do j = 0, nterm
         if (k .ne. j) then
            sum = sum + d(k, j)
         else
            sum = sum
         endif
      enddo

      d(k, k) = -sum
   enddo

   ! Second derivative:
   ! Off-diagonal elements Eqn(13)
   do k = 0, nterm
      xk = x(k)

      do j = 0, nterm
         xj = x(j)

         if (k .ne. j) then
            d2(k, j) = 2.0_dp * (d(k, k) * d(k, j) - d(k, j) / (xk - xj))
         else
            d2(k, j) = 0.0_dp
         endif
      enddo
   enddo

   ! Diagonal elements Eqn(9) again
   do k = 0, nterm
      sum = 0.0_dp

      do j = 0, nterm

         if (k .ne. j) then
            sum = sum + d2(k, j)
         else
            sum = sum
         endif
      enddo

      d2(k, k) = -sum
   enddo

   ! Third derivative:
   ! Off-diagonal elements
   do k = 0, nterm
      xk = x(k)

      do j = 0, nterm
         xj = x(j)

         if (k .ne. j) then
            d3(k, j) = 3.0_dp * (d2(k, k) * d(k, j) - d2(k, j) / (xk - xj))
         else
            d3(k, j) = 0.0_dp
         endif
      enddo
   enddo

    ! Diagonal elements
    do k = 0, nterm
       sum = 0.

       do j = 0, nterm
          if (k .ne. j) then
             sum = sum + d3(k, j)
          else
             sum = sum
          endif
       enddo

       d3(k, k) = -sum
    enddo

return

end subroutine derv
