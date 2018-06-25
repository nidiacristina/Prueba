subroutine sync_ranks( left, right, nbuffer )
! Takes the left/right arrays of dim nbuffer and sends the left to the rank on
! the left and the right to the rank on the right.  On exit, left contains
! information from the left rank, and right contains information from the
! right rank.

  use constants, only: nprocs, rank
  use mpi, only:       MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_STATUSES_IGNORE
  use precision, only: dp

  implicit none

  real(kind=dp), dimension(1:nbuffer), intent(inout) :: left
  real(kind=dp), dimension(1:nbuffer), intent(inout) :: right
  integer, intent(in)                                :: nbuffer

  real(kind=dp), dimension(1:nbuffer)                :: from_left, from_right
  integer                                            :: tag2left = 1, tag2right = 2
  integer                                            :: ierr
  integer                                            :: req2left_sent, req2left_recv, req2right_sent, req2right_recv
  integer, dimension(1:2)                            :: allbutN, allbut0
  integer, dimension(1:4)                            :: middleranks

  ! Post the pass to the right.
  if ( rank > 0 ) then
     call MPI_IRECV( from_left, nbuffer, MPI_DOUBLE_PRECISION, &
                     rank - 1, tag2right, MPI_COMM_WORLD, req2right_recv, ierr )
  endif
  if ( rank < nprocs - 1 ) then
     call MPI_ISEND( right, nbuffer, MPI_DOUBLE_PRECISION, &
                     rank + 1, tag2right, MPI_COMM_WORLD, req2right_sent, ierr )
  endif

  ! Post the pass to the left.
  if ( rank < nprocs - 1 ) then
     call MPI_IRECV( from_right, nbuffer, MPI_DOUBLE_PRECISION, &
                     rank + 1, tag2left, MPI_COMM_WORLD, req2left_recv, ierr )
  endif
  if ( rank > 0 ) then
     call MPI_ISEND( left, nbuffer, MPI_DOUBLE_PRECISION, &
                     rank - 1, tag2left, MPI_COMM_WORLD, req2left_sent, ierr )
  endif

  ! Wait on the left/right sends.
  if ( rank == 0 .and. nprocs > 1 ) then
     allbutN(1) = req2right_sent
     allbutN(2) = req2left_recv
     call MPI_WAITALL( 2, allbutN, MPI_STATUSES_IGNORE, ierr )
  elseif ( rank == nprocs - 1 .and. nprocs > 1 ) then
     allbut0(1) = req2left_sent
     allbut0(2) = req2right_recv
     call MPI_WAITALL( 2, allbut0, MPI_STATUSES_IGNORE, ierr )
  elseif ( nprocs > 1 ) then
     middleranks(1) = req2right_sent
     middleranks(2) = req2left_recv
     middleranks(3) = req2left_sent
     middleranks(4) = req2right_recv
     call MPI_WAITALL( 4, middleranks, MPI_STATUSES_IGNORE, ierr )
  endif

  ! Update the buffers.
  left  = from_left
  right = from_right

end subroutine sync_ranks

subroutine sync_ranks_complex( left, right, nbuffer )
! Takes the left/right arrays of dim nbuffer and sends the left to the rank on
! the left and the right to the rank on the right.  On exit, left contains
! information from the left rank, and right contains information from the
! right rank.

  use constants, only: nprocs, rank
  use mpi, only:       MPI_COMM_WORLD, MPI_DOUBLE_COMPLEX, MPI_STATUSES_IGNORE
  use precision, only: dp

  implicit none

  complex(kind=dp), dimension(1:nbuffer), intent(inout) :: left
  complex(kind=dp), dimension(1:nbuffer), intent(inout) :: right
  integer, intent(in)                                   :: nbuffer

  complex(kind=dp), dimension(1:nbuffer)                :: from_left, from_right
  integer                                               :: tag2left = 1, tag2right = 2
  integer                                               :: ierr
  integer                                               :: req2left_sent, req2left_recv, req2right_sent, req2right_recv
  integer, dimension(1:2)                               :: allbutN, allbut0
  integer, dimension(1:4)                               :: middleranks

  ! Post the pass to the right.
  if ( rank > 0 ) then
     call MPI_IRECV( from_left, nbuffer, MPI_DOUBLE_COMPLEX, &
                     rank - 1, tag2right, MPI_COMM_WORLD, req2right_recv, ierr )
  endif
  if ( rank < nprocs - 1 ) then
     call MPI_ISEND( right, nbuffer, MPI_DOUBLE_COMPLEX, &
                     rank + 1, tag2right, MPI_COMM_WORLD, req2right_sent, ierr )
  endif

  ! Post the pass to the left.
  if ( rank < nprocs - 1 ) then
     call MPI_IRECV( from_right, nbuffer, MPI_DOUBLE_COMPLEX, &
                     rank + 1, tag2left, MPI_COMM_WORLD, req2left_recv, ierr )
  endif
  if ( rank > 0 ) then
     call MPI_ISEND( left, nbuffer, MPI_DOUBLE_COMPLEX, &
                     rank - 1, tag2left, MPI_COMM_WORLD, req2left_sent, ierr )
  endif

  ! Wait on the left/right sends.
  if ( rank == 0 .and. nprocs > 1 ) then
     allbutN(1) = req2right_sent
     allbutN(2) = req2left_recv
     call MPI_WAITALL( 2, allbutN, MPI_STATUSES_IGNORE, ierr )
  elseif ( rank == nprocs - 1 .and. nprocs > 1 ) then
     allbut0(1) = req2left_sent
     allbut0(2) = req2right_recv
     call MPI_WAITALL( 2, allbut0, MPI_STATUSES_IGNORE, ierr )
  elseif ( nprocs > 1 ) then
     middleranks(1) = req2right_sent
     middleranks(2) = req2left_recv
     middleranks(3) = req2left_sent
     middleranks(4) = req2right_recv
     call MPI_WAITALL( 4, middleranks, MPI_STATUSES_IGNORE, ierr )
  endif

  ! Update the buffers.
  left  = from_left
  right = from_right

end subroutine sync_ranks_complex

subroutine sync_flag( flag )
! Synchronizes an integer flag on the root rank to every other ranks.  Useful
! for detecting error conditions (e.g. I/O failures) and coordinating error
! handling on ranks where the error did not occur.

  use constants, only:         root
  use mpi, only:               MPI_COMM_WORLD, MPI_INTEGER

  implicit none

  integer, intent(inout) :: flag

  integer                :: ierr

  call MPI_BCAST( flag, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr )

end subroutine sync_flag

subroutine smpm_assert( condition, message )
! Asserts that the supplied condition is true.  If it is not, the root rank
! writes the supplied message via notify() and terminates the program on every
! rank.  No synchronization is done if the condition is true.

  use constants, only: rank, root

  implicit none

  logical, intent(in) :: condition
  character(len=*)    :: message

  integer             :: ierr

  if ( condition .eqv. .false. ) then
     if ( rank == root ) call notify( message )
     call MPI_FINALIZE( ierr )
     stop
  end if

end subroutine smpm_assert

subroutine smpm_assert_code( condition, code, message )
! Asserts that the supplied condition is true.  If it is not, the root rank
! writes the supplied message via notify() and terminates the program on every
! rank.  No synchronization is done if the condition is true.  The program's
! return code is supplied by the caller.

  use constants, only: rank, root

  implicit none

  logical, intent(in)          :: condition
  integer, intent(in)          :: code
  character(len=*), intent(in) :: message

  integer                      :: ierr

  if ( condition .eqv. .false. ) then
     if ( rank == root ) call notify( message )
     call MPI_FINALIZE( ierr )
     call stop_wrapper( code )
  end if

end subroutine smpm_assert_code

subroutine stop_wrapper( stop_code )
! Wrapper around STOP that halts execution with a user-specified stop code.
! This exists to work around the silliness of Fortran 2008 which allows
! non-string stop codes, but does not allow non-constant stop codes.
!
! NOTE: As of 2017/07/23 it appears that the Fortran 2008 only allows
!       constant expressions to be accepted by STOP.  This is based on
!       the statement in "Modern Fortran explained" section 20.1.6:
!
!         "The stop statement now accepts any default integer or default
!          character scalar constant expression as the stop code, instead
!          of only simple literals."
!
!       Combined with the fact that Intel's 2017 Fortran compiler complains
!       about:
!
!         integer :: stop_code = 0
!         ...
!         STOP stop_code
!

  implicit none

  integer :: stop_code

  ! Stop with the same value provided when stop_code is in [0, 20], otherwise
  ! default to 21.
  !
  ! NOTE: The choice to support 21 distinct values is mostly arbitrary though
  !       informed by the current callers which involves the main solver and
  !       validation suite, where we could have a stop code equal to the number
  !       of validation cases (which is 14 as of 2017/07/23).  We aim to
  !       support current needs and give a little head room, which is unlikely
  !       to be exercised.  Should it be, a warning message is generated before
  !       clipping the code to something supported.
  !
  select case (stop_code)
  case (0)
     stop 0
  case (1)
     stop 1
  case (2)
     stop 2
  case (3)
     stop 3
  case (4)
     stop 4
  case (5)
     stop 5
  case (6)
     stop 6
  case (7)
     stop 7
  case (8)
     stop 8
  case (9)
     stop 9
  case (10)
     stop 10
  case (11)
     stop 11
  case (12)
     stop 12
  case (13)
     stop 13
  case (14)
     stop 14
  case (15)
     stop 15
  case (16)
     stop 16
  case (17)
     stop 17
  case (18)
     stop 18
  case (19)
     stop 19
  case (20)
     stop 20
  case default
     write( *, '(A,I0,A)' ) 'Requested a stop code of ', stop_code, &
                            ' though STOP_WRAPPER() only supports codes up to 20.'
     stop 21
  end select

end subroutine stop_wrapper
