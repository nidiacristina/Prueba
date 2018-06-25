module sponge
! Contains the variables related to the sponge layer.

  use precision, only: dp

  implicit none
  save

  ! Sponge layer variable with the coefficients.
  real(kind=dp), allocatable, dimension(:, :) :: raycoeff

  ! Initialize the location of the sponge layer
  integer, dimension(1:4)                     :: sponge_layer_location = (/0, 0, 0, 0/)

  ! Set the percentage for the thickness of the sponge at the
  ! left and right wall. The default is 5% of the total domain length
  real(kind=dp)                               :: left_fraction  = 0.05_dp
  real(kind=dp)                               :: right_fraction = 0.05_dp

  ! Time scale of the absorbing layer: this parameter is unique to the simulation
  ! and needs to be computed separately. In the case of the propagating ISW, it is
  ! 2.0_dp * cph / (dxL + dxR) * (-1.0_dp) * log( eps )
  ! where dxL =  left_fraction * Lx & dxR = right_fraction * Lx
  !
  ! Initialized to 0 on purposed, to force the user to set it on the input
  ! file, if he or she is requesting a sponge.
  real(kind=dp)                               :: time_scale_in_x = 0.0_dp

end module sponge
