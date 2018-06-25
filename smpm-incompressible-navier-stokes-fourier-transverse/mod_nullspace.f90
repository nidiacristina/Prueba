module nullspace
! Contains constants related to the grid, physics, and numerics.

  use precision, only: dp

  implicit none
  save

  integer                                        :: null_dim
  real(kind=dp), allocatable, dimension(:, :, :) :: null_basis

end module nullspace
