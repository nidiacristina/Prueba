module geom
! Contains the grid.

  use precision, only: dp

  implicit none
  save

  real(kind=dp), allocatable, dimension(:, :) :: cx, cz

end module geom
