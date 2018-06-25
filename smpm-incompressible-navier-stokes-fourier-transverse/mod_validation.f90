module validation

  use precision, only:    dp

  implicit none
  save

  integer :: ADVECTION_TERMS                         = 1
  integer :: FFTW                                    = 2
  integer :: DIVERGENCE                              = 3
  integer :: DOUBLE_CURL                             = 4
  integer :: FILTER_Y                                = 5
  integer :: GRADIENTS                               = 6
  integer :: HELMHOLTZ_2D_SCALAR_DIRICHLET           = 7
  integer :: HELMHOLTZ_3D_SCALAR_DIRICHLET           = 8
  integer :: HELMHOLTZ_3D_SCALAR_NEUMANN             = 9
  integer :: HELMHOLTZ_3D_VECTOR_DIRICHLET           = 10
  integer :: HELMHOLTZ_3D_VECTOR_NEUMANN             = 11
  integer :: LAPLACIAN                               = 12
  integer :: POISSON_3D_ITERATIVELY                  = 13
  integer :: TRANSPORT_TERM                          = 14


  character(len=64), parameter :: validation_routines(14) =  &
           (/ 'advection_terms                        ', &
              'fftw                                   ', &
              'divergence                             ', & 
              'double_curl                            ', & 
              'filter_y                               ', &
              'gradients                              ', & 
              'helmholtz_2D_scalar_dirichlet          ', &
              'helmholtz_3D_scalar_dirichlet          ', &
              'helmholtz_3D_scalar_neumann            ', &
              'helmholtz_3D_vector_dirichlet          ', &
              'helmholtz_3D_vector_neumann            ', &
              'laplacian                              ', & 
              'poisson_3D_iteratively                 ', &
              'transport_term                         ' /)
  ! Internal state needed for validation.
  !
  ! NOTE: This should only be used for validating routines that have a fixed
  !       interface and used as values for other routines (e.g. operators that
  !       are applied within GMRES).
  !
  complex(kind=dp) :: wave        ! stores the wave number for validating GMRES routines.

end module validation
