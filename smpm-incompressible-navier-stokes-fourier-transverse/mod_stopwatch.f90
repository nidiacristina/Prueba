module stopwatch
! Timing module.  Contains routines to compute time spent in a section of code
! by marking a reference point (via stopwatch_tick()) and computing the elapsed
! time (via stopwatch_get_elapsed()).  Prior to use, the module must be
! initialized (via stopwatch_initialize()).

  use precision, only: int64

  implicit none
  save

  ! Rate at which the clock advances on this system.  This is necessary to
  ! convert clock ticks into actual time.
  !
  ! NOTE: We explicitly use 64-bit integers to ensure the highest precision
  !       timer available is used by SYSTEM_CLOCK().
  !
  ! NOTE: We explicitly initialize this to 0 so we get NaN's if the module
  !       isn't initialized properly (via initialize_stopwatches()).  That
  !       should be noticeable to someone without adding overhead to check
  !       for initialization.
  integer(kind=int64) :: clock_rate = 0_int64

  contains

    subroutine stopwatch_initialize()
    ! Initializes the Stopwatches module so that its routines may be called.
    ! This *must* be called before any other subroutines or functions from this
    ! module otherwise invalid results will be computed.

      implicit none

      ! Determine how finely we can time things, though only do it once
      ! should someone accidentally call this routine multiple times.
      if (clock_rate == 0) then
         call SYSTEM_CLOCK( COUNT_RATE=clock_rate )
      end if

    end subroutine stopwatch_initialize

    subroutine stopwatch_tick( watch )
    ! Ticks the stopwatch and records the current time.  This can be used as a
    ! reference point in time so as to compute elapsed time via
    ! stopwatch_get_elapsed().  The Stopwatches module must be initialized via
    ! initialize_stopwatches() before this function may be called.

      use precision, only: int64

      implicit none

      integer(kind=int64), intent(inout) :: watch

      call SYSTEM_CLOCK( COUNT=watch )

    end subroutine stopwatch_tick

    function stopwatch_get_elapsed( before ) result( elapsed )
    ! Computes the elapsed time since the supplied stopwatch tick.  The
    ! Stopwatches module must be initialized via initialize_stopwatches()
    ! before this function may be called.

      use precision, only: dp, int64

      implicit none

      integer(kind=int64) :: before
      real(kind=dp)       :: elapsed

      integer(kind=int64) :: now

      ! get the current clock value and compute the elapsed time based on the
      ! previously recorded rate.
      call SYSTEM_CLOCK( COUNT=now )

      elapsed = real( (now - before), kind=dp ) / clock_rate

    end function stopwatch_get_elapsed

end module stopwatch
