module command_line
! Contains routines for parsing command lines.
!
! NOTE: This module is Fortran 2003.

  implicit none
  save

  contains

    subroutine get_command_line_double( argument_index, argument, success_flag )
    ! Converts a command line argument into a double precision floating point
    ! value.  Detection of a successful conversion from string to floating point
    ! is optional.

      use precision, only: dp

      implicit none

      integer, intent(in)            :: argument_index
      real(kind=dp), intent(out)     :: argument
      logical, intent(out), optional :: success_flag

      character(len=:), allocatable  :: argument_string
      integer                        :: ierr

      call get_command_line_string( argument_index, argument_string )
      read( argument_string, *, iostat=ierr ) argument

      ! Pretend we parsed zero if we didn't get a numeric argument.
      if (ierr /= 0) then
         argument = 0.0_dp
      end if

      ! Update our flag only if it was provided by the caller.
      if ( present( success_flag ) ) then
         success_flag = (ierr == 0)
      end if

    end subroutine get_command_line_double

    subroutine get_command_line_string( argument_index, argument )
    ! Acquires the command line argument specified by an index.  The memory
    ! needed to hold a non-trimmed string is allocated automatically.

      implicit none

      integer, intent(in)                        :: argument_index
      character(len=:), allocatable, intent(out) :: argument

      integer                                    :: argument_length

      ! NOTE: GET_COMMAND_ARGUMENT() doesn't take a deferred length argument so
      !       have to call it once to get the length of our argument, allocate
      !       the proper amount of space, and then call it again to get the
      !       actual argument.
      call GET_COMMAND_ARGUMENT( argument_index, length=argument_length )
      allocate( character(argument_length) :: argument )

      call GET_COMMAND_ARGUMENT( argument_index, value=argument )

    end subroutine get_command_line_string

    subroutine get_command_line_strings( start_index, end_index, strings )
    ! Acquires one or more command line arguments and returns an array of scalar
    ! character buffers, each sized to hold the largest command line argument.
    ! If the start index is greater than the end index, nothing is done and the
    ! supplied array remains unchanged.

      implicit none

      integer, intent(in)                                        :: start_index
      integer, intent(in)                                        :: end_index
      character(len=:), allocatable, dimension(:), intent(inout) :: strings

      integer                                                    :: string_index
      integer                                                    :: longest_length
      integer                                                    :: current_length

      ! There is nothing to do if we weren't given a range of indices to work
      ! with.
      if ( start_index > end_index ) then
         return
      end if

      ! Find the longest string so we can allocate our strings.
      longest_length = 0
      do string_index = start_index, end_index
         call GET_COMMAND_ARGUMENT( string_index, length=current_length )
         longest_length = MAX( longest_length, current_length )
      end do

      ! Create a buffer for all of the strings, padded out to the longest length
      ! present.
      allocate( character(longest_length) :: strings(end_index - start_index + 1) )

      do string_index = start_index, end_index
         call GET_COMMAND_ARGUMENT( string_index, &
                                    strings(string_index - start_index + 1) )
      end do

    end subroutine get_command_line_strings

end module command_line
