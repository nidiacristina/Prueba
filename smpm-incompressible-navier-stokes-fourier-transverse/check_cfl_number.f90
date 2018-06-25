subroutine check_cfl_number( dt_in, dt_out )
! Updates the maximum CFLs based on the current velocities.  If the solver has
! been configured for adaptive timestepping, the timestep is updated as
! necessary.

  use constants, only:               nsuby
  use field_variables, only:         ubc, ux0, uy0, uz0
  use mesh_deformation_maps, only:   delta_x, delta_z
  use mpi, only:                     MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_MAX
  use options, only:                 adaptive_timestep, &
                                     cflxmin_lim, cflxmax_lim, &
                                     cflymin_lim, cflymax_lim, &
                                     cflzmin_lim, cflzmax_lim
  use precision, only:               dp
  use timestepping, only:            cfl_max_x, cfl_max_y, &
                                     cfl_max_z, logflag, &
                                     time_ndx
  use transverse, only:              delta_y
  use woodbury_matrices, only:       rpk

  implicit none

  ! Set the output variable.
  real(kind=dp), intent(in)                            :: dt_in
  real(kind=dp), intent(out)                           :: dt_out

  ! Arrays that store CFL # per transverse wavemode
  real(kind=dp), allocatable, dimension(:)             :: cfl_x_transverse, cfl_y_transverse, cfl_z_transverse 

  ! Internal variables.
  real(kind=dp)                                        :: cfl_x, cfl_y, cfl_z

  ! Flags.
  logical                                              :: change_dt, exceed_limit, below_limit

  ! Timestep change factor.
  real(kind=dp)                                        :: cdt

  ! Buffers for messages written to standard output.
  character(len=64)                                    :: caststr, caststr2, caststr3
  
  ! Other variables
  integer                                              :: ierr
  integer                                              :: ii
  
  ! Default change factor to handle the case where we're not adapting our
  ! timestep.
  cdt = 1.0_dp

  ! Compute the local CFL # (i.e. per rank)
  if ( nsuby > 1 ) then

     allocate( cfl_x_transverse(1:nsuby), cfl_y_transverse(1:nsuby), cfl_z_transverse(1:nsuby) )

     ! Find Maximum CFL # per transverse wavemode
     do ii = 1,nsuby
        cfl_x_transverse( ii ) = maxval( abs( dt_in *  ( ux0(:, ii) + ubc( :, ii) )  / delta_x ) )
        cfl_y_transverse( ii ) = maxval( abs( dt_in * uy0( :, ii ) / delta_y ) )
        cfl_z_transverse( ii ) = maxval( abs( dt_in * uz0( :, ii ) / delta_z ) ) 
     enddo

     ! Store the maximum CFL Value
     cfl_x   = maxval( cfl_x_transverse ) 
     cfl_y   = maxval( cfl_y_transverse ) 
     cfl_z   = maxval( cfl_z_transverse ) 

     deallocate( cfl_x_transverse, cfl_y_transverse, cfl_z_transverse )

  else

     cfl_x   = maxval( abs( dt_in * reshape( (ux0 + ubc), (/rpk * nsuby/) ) / delta_x ) ) 
     cfl_y   = 0.0_dp
     cfl_z   = maxval( abs( dt_in * reshape(  uz0, (/rpk * nsuby/) ) / delta_z ) ) 

  endif

  ! Compute the maximum CFL #, in  x, over all ranks and distribute across all ranks 
  call MPI_ALLREDUCE( cfl_x, cfl_max_x, 1, MPI_DOUBLE_PRECISION, &
       MPI_MAX, MPI_COMM_WORLD, ierr )
  
  ! Compute the maximum CFL #, in  y, over all ranks and distribute across all ranks
  call MPI_ALLREDUCE( cfl_y, cfl_max_y, 1, MPI_DOUBLE_PRECISION, &
       MPI_MAX, MPI_COMM_WORLD, ierr )

  ! Compute the maximum CFL #, in  z, over all ranks and distribute across all ranks
  call MPI_ALLREDUCE( cfl_z, cfl_max_z, 1, MPI_DOUBLE_PRECISION, &
       MPI_MAX, MPI_COMM_WORLD, ierr )

  !Note: We use ALLREDUCE because we want each rank to have a stored value of the max cfl #

  ! Announce out the CFL information.
  write( caststr, '(D17.10)') cfl_max_x
  write( caststr2, '(D17.10)') cfl_max_y
  write( caststr3, '(D17.10)') cfl_max_z

  call notify_cond( logflag, ' ' )
  call notify_cond( logflag, '      CFL # in x: ' // caststr  )
  call notify_cond( logflag, '      CFL # in y: ' // caststr2 )
  call notify_cond( logflag, '      CFL # in z: ' // caststr3 )

  ! Determine if we need to adjust our timestep if requested.
  if ( adaptive_timestep ) then

     ! Determine if the maximum CFL is bounded in x, y and z.
     exceed_limit = ((cfl_max_z > cflzmax_lim) .or. (cfl_max_y > cflymax_lim) .or. (cfl_max_x > cflxmax_lim))
     below_limit  = ((cfl_max_z < cflzmin_lim) .or. (cfl_max_y < cflymin_lim) .or. (cfl_max_x < cflxmin_lim))
     change_dt    = (exceed_limit .or. below_limit)

     ! If CFL is out of bounds we need to compute the appropriate change
     ! factor.
     if ( change_dt ) then
        if (exceed_limit) then
           cdt = 5.0_dp/4.0_dp
           write( caststr, '(D17.10)') dt_in / cdt
           write( caststr2, '(I17.10)') time_ndx

           call notify( ' ' )
           call notify( '      TIMESTEP SIZE *REDUCED* BY 4/5 TO: ' // caststr )
           call notify( '      AT TIMESTEP: ' // caststr2 )
        endif

        if (below_limit) then
           cdt = 0.8_dp
           write( caststr, '(D17.10)') dt_in / cdt
           write( caststr2, '(I17.10)') time_ndx

           call notify( ' ' )
           call notify( '      TIMESTEP *INCREASED* BY 5/4 TO: ' // caststr )
           call notify( '      AT TIMESTEP: ' // caststr2 )
        endif

     endif

  endif

  ! Update the time step.
  dt_out = dt_in / cdt

end subroutine check_cfl_number
