module legendre
! This module defines the Legendre collocation points as well as the
! differentiation matrices associated to those points.

  use precision, only: dp

  implicit none
  save

  real(kind=dp), allocatable, dimension(:)    :: points    ! Array with Gauss-Lobatto-Legendre collocation points.
  real(kind=dp), allocatable, dimension(:)    :: wg        ! Weights for filtering.
  real(kind=dp), allocatable, dimension(:, :) :: filter_xz ! Filtering matrix in the x-z plane.
  real(kind=dp), allocatable, dimension(:, :) :: D, D2, D3 ! 1st, 2nd and 3rd order differentiation matrices.

end module legendre
