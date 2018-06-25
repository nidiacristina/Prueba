subroutine notify( nstring )
! If this is the root rank, pre-pend a date string to the supplied message and
! write it to standard output.  The date string is of the form:
!
!    29-Mar-2013 12:41:33: nstring
!

  use constants, only: rank, root

  implicit none

  character(len=*), intent(in) :: nstring

  character(len=100)           :: date, time
  integer, dimension(8)        :: T
  character(len=10)            :: zone
  character(len=3)             :: month

  ! Nothing to do if we're not the one notifying.
  if ( (rank /= root) ) then
     return
  end if

  call DATE_AND_TIME( date, time, zone, T )

  select case (T(2))
     case (1)
        month = 'Jan'
     case (2)
        month = 'Feb'
     case (3)
        month = 'Mar'
     case (4)
        month = 'Apr'
     case (5)
        month = 'May'
     case (6)
        month = 'Jun'
     case (7)
        month = 'Jul'
     case (8)
        month = 'Aug'
     case (9)
        month = 'Sep'
     case (10)
        month = 'Oct'
     case (11)
        month = 'Nov'
     case (12)
        month = 'Dec'
   end select

   write(*,*) date(7:8), '-', month, '-', date(1:4), ' ', time(1:2), ':', time(3:4), ':', time(5:6), ':  ', trim( nstring )

end subroutine notify

subroutine notify_cond( condition, nstring )
! Conditionally calls notify() with the supplied string.  Useful for filtering
! notifications depending on how the solver was invoked.

  implicit none

  logical, intent(in)          :: condition
  character(len=*), intent(in) :: nstring

  if ( condition ) then
     call notify( nstring )
  end if

end subroutine notify_cond
