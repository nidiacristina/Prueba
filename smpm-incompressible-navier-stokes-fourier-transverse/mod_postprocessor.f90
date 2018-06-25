module postprocessor
! Contains the variables related to the postprocessor.

  use precision, only: dp

  implicit none
  save

  ! This is the variable that corresponds to the field time and shall be included
  ! as the time in the postprocessor field.
  real(kind=dp)                               :: field_time = 0.0_dp

  ! Request computation of streamfunction
  logical                                     :: compute_streamfunction = .false.
  
  ! Request computation of pressure
  logical                                     :: compute_pressure       = .false.

end module postprocessor
