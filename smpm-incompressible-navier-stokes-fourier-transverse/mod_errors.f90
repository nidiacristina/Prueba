module errors
! Contains errors related to the solution, including both setup and per-timestep
! solutions.

  use precision, only: dp

  implicit none
  save

  ! GMRES iteration counts for each of the physics problems.
  integer          :: gmres_diffusion_iterations_real
  integer          :: gmres_diffusion_iterations_imag
  integer          :: gmres_poisson_iterations_real
  integer          :: gmres_poisson_iterations_imag
  integer          :: gmres_viscous_iterations_real
  integer          :: gmres_viscous_iterations_imag

  ! Solver setup errors.

  ! Null space computation errors.
  real(kind=dp)    :: l2_poisson_kernel_error
  real(kind=dp)    :: l2_poisson_schur_kernel_error

  ! Per-timestep errors.

  ! Solution errors for each of the physics problems.
  complex(kind=dp) :: linf_diffusion_error
  real(kind=dp)    :: linf_divergence_error
  real(kind=dp)    :: l2_poisson_error
  complex(kind=dp) :: l2_poisson_schur_error
  complex(kind=dp) :: linf_viscous_error

  contains

    function compute_numeric_error() result( numeric_error )

      use precision, only: dp

      implicit none

      real(kind=dp) :: numeric_error

      numeric_error = (abs( linf_diffusion_error ) + &
                       linf_divergence_error + &
                       l2_poisson_error + &
                       abs( l2_poisson_schur_error ) + &
                       abs( linf_viscous_error ))

    end function compute_numeric_error

end module errors
