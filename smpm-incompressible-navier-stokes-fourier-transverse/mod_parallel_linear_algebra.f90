module parallel_linear_algebra
! Contains functions for doing parallel linear algebra on distributed arrays.

  implicit none

contains

  function pnorm2( x ) result( normx )
  ! Compute a parallel distributed Euclidean norm.

     use mpi, only:       MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_SUM
     use precision, only: dp

     implicit none

     real(kind=dp), dimension(:), intent(in) :: x
     real(kind=dp)                           :: normx

     integer                                 :: ierr

     real(kind=dp), external                 :: DDOT

     normx = DDOT( size( x, 1 ), x(1), 1, x(1), 1 )
     call MPI_ALLREDUCE( MPI_IN_PLACE, normx, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
     normx = sqrt( normx )

  end function pnorm2

  function znorm2( x ) result( normx )
  ! Compute a parallel distributed complex Euclidean norm.

     use mpi, only:       MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_SUM
     use precision, only: dp

     implicit none

     complex(kind=dp), dimension(:), intent(in) :: x
     real(kind=dp)                              :: normx

     integer                                    :: ierr

     complex(kind=dp), external                 :: ZDOTC

     normx = real( ZDOTC( size( x, 1 ), x(1), 1, x(1), 1 ), kind=dp )
     call MPI_ALLREDUCE( MPI_IN_PLACE, normx, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
     normx = sqrt( normx )

  end function znorm2

  function pdot_product( x, y ) result( xdoty )
  ! Compute a parallel distributed dot product.

     use mpi, only:       MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_SUM
     use precision, only: dp

     implicit none

     real(kind=dp), dimension(:), intent(in) :: x
     real(kind=dp), dimension(:), intent(in) :: y
     real(kind=dp)                           :: xdoty

     real(kind=dp)                           :: iixdoty
     integer                                 :: ierr

     real(kind=dp), external                 :: DDOT

     iixdoty = DDOT( size( x, 1 ), x(1), 1, y(1), 1 )
     call MPI_ALLREDUCE( iixdoty, xdoty, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )

  end function pdot_product

  function zdot_product( x, y ) result( xdoty )
  ! Compute a complex parallel distributed dot product.

     use mpi, only:       MPI_COMM_WORLD, MPI_DOUBLE_COMPLEX, MPI_SUM
     use precision, only: dp

     implicit none

     complex(kind=dp), dimension(:), intent(in) :: x
     complex(kind=dp), dimension(:), intent(in) :: y
     complex(kind=dp)                           :: xdoty

     complex(kind=dp)                           :: iixdoty
     integer                                    :: ierr

     complex(kind=dp), external                 :: ZDOTC

     iixdoty = ZDOTC( size( x, 1 ), x(1), 1, y(1), 1 )
     call MPI_ALLREDUCE( iixdoty, xdoty, 1, &
                         MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr )

  end function zdot_product

  function assign_value_at( x, a, z, ndx, val )
  ! Parallel distributed assignment.

     use precision, only: dp

     implicit none

     ! NOTE: We have to order our definitions so that the compiler can define
     !       each variable when it is encountered, rather than by maintaining
     !       the order in the argument list.  If we don't do this, the build
     !       with ifort will break.
     integer, intent(in)                              :: a         ! Ownership range of this rank.
     integer, intent(in)                              :: z         ! Ownership range of this rank.
     real(kind=dp), dimension(1:z-a+1), intent(inout) :: x         ! Local block of the distributed array.
     integer, intent(in)                              :: ndx       ! Index in the full array we want to
                                                                   ! Assign a value to.
     real(kind=dp), intent(in)                        :: val       ! Value we want to assign.

     real(kind=dp), dimension(1:z-a+1)                :: assign_value_at
     integer                                          :: local_ndx

     ! If the modification is in this rank's range, do it, otherwise return
     ! the same array.
     if ( ndx >= a .and. ndx <= z ) then
        ! Get the local index and assign.
        local_ndx    = ndx - a + 1
        x(local_ndx) = val
     endif

     ! Return.
     assign_value_at = x

  end function assign_value_at

  function assign_complex_value_at( x, a, z, ndx, val )
  ! Parallel distributed complex assignment.

     use precision, only: dp

     implicit none

     ! NOTE: We have to order our definitions so that the compiler can define
     !       each variable when it is encountered, rather than by maintaining
     !       the order in the argument list.  If we don't do this, the build
     !       with ifort will break.
     integer, intent(in)                                 :: a         ! Ownership range of this rank.
     integer, intent(in)                                 :: z         ! Ownership range of this rank.
     complex(kind=dp), dimension(1:z-a+1), intent(inout) :: x         ! Local block of the distributed array.
     integer, intent(in)                                 :: ndx       ! Index in the full array we want to
                                                                      ! Assign a value to.
     complex(kind=dp), intent(in)                        :: val       ! Value we want to assign.

     complex(kind=dp), dimension(1:z-a+1)                :: assign_complex_value_at
     integer                                             :: local_ndx

     ! If the modification is in this rank's range, do it, otherwise return
     ! the same array.
     if ( ndx >= a .and. ndx <= z ) then
        ! Get the local index and assign.
        local_ndx    = ndx - a + 1
        x(local_ndx) = val
     endif

     ! Return.
     assign_complex_value_at = x

  end function assign_complex_value_at

  function pmaxval( x ) result( xmax )
  ! Compute a parallel distributed maxval().

     use mpi, only:       MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_MAX
     use precision, only: dp

     implicit none

     real(kind=dp), dimension(:), intent(in) :: x
     real(kind=dp)                           :: xmax

     real(kind=dp)                           :: iixmax
     integer                                 :: ierr

     iixmax = maxval( x )
     xmax   = 0.0_dp
     call MPI_ALLREDUCE( iixmax, xmax, 1, &
                         MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  end function pmaxval

  function pminval( x ) result( xmin )
  ! Compute a parallel distributed minval().

    use mpi, only:       MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_MIN
    use precision, only: dp

    implicit none

    real(kind=dp), dimension(:), intent(in) :: x
    real(kind=dp)                           :: xmin

    real(kind=dp)                           :: iixmin
    integer                                 :: ierr

    iixmin = minval( x )
    xmin   = 0.0_dp
    call MPI_ALLREDUCE( iixmin, xmin, 1, &
                        MPI_DOUBLE_PRECISION, MPI_Min, MPI_COMM_WORLD, ierr )

  end function pminval

  function zmaxval( x ) result( xmax )
  ! Compute a parallel distributed complex MAX, and returns the magnitude of
  ! the maximal element.

     use mpi, only:       MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_MAX
     use precision, only: dp

     implicit none

     complex(kind=dp), dimension(:), intent(in) :: x
     real(kind=dp)                              :: xmax

     real(kind=dp)                              :: iixmax
     integer                                    :: ierr

     iixmax = maxval( abs( x ) )
     xmax   = 0.0_dp
     call MPI_ALLREDUCE( iixmax, xmax, 1, &
                         MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  end function zmaxval

   function apply_nullspace_projection( n, x, u )
   ! Computes the null space projection of x out of u for the real and
   ! imaginary parts of x seperately.

     use precision, only: dp

     implicit none

     integer, intent(in)                             :: n
     complex(kind=dp), dimension(1:n), intent(inout) :: x
     real(kind=dp), dimension(1:n), intent(in)       :: u

     real(kind=dp), dimension(1:n)                   :: real_buff, imag_buff
     complex(kind=dp), dimension(1:n)                :: apply_nullspace_projection

     real_buff                  = real( x, kind=dp ) - u * pdot_product( real( x, kind=dp ), u )
     imag_buff                  = aimag( x )        - u * pdot_product( aimag( x ), u )
     apply_nullspace_projection = cmplx( real_buff, imag_buff, kind=dp )

   end function apply_nullspace_projection

end module parallel_linear_algebra
