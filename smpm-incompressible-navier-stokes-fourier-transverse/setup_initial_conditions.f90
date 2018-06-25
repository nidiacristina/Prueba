subroutine setup_initial_conditions()
! configure the solver's initial conditions.  if we're restarting from a
! previous checkpoint, we acquire them from the restart file, otherwise we load
! them from the init file.

  use field_variables, only: rho, rho0, rho1, rho2, rho_b, &
                             ux, ux0, ux1, ux2, ux_b, &
                             uy, uy0, uy1, uy2, uy_b, &
                             uz, uz0, uz1, uz2, uz_b
  use options, only:         apply_restart, fname_init, fname_restart, &
                             read_bcs_from_initfile
  use timestepping, only:    dt, dt1, dt2
  implicit none

  logical :: io_flag

  if ( apply_restart ) then
     call notify( 'Reading the restart file.' )

     ! Read the initial conditions from the restart file.
     call read_restartfile_data( fname_restart, io_flag )
     call smpm_assert( io_flag, 'Failed to open ' // fname_restart )

     ! Set the velocity and density field to the most recent output
     ! Note: We do this to avoid storing both fields in the restart
     !       file, since both are technically the same. (see advance
     !       field variables routine). Thus, we reduce the number of
     !       arrays in the restart file by 4.
     ux  = ux0
     uy  = uy0
     uz  = uz0
     rho = rho0

  else
     call notify( 'Reading the initialization file.' )

     ! Read the initial conditions.
     call read_initfile_data( fname_init, io_flag )
     call smpm_assert( io_flag, 'Failed to open ' // fname_init )

     ux0  = ux
     ux1  = ux
     ux2  = ux
     uy0  = uy
     uy1  = uy
     uy2  = uy
     uz0  = uz
     uz1  = uz
     uz2  = uz
     rho0 = rho
     rho1 = rho
     rho2 = rho
     dt1  = dt
     dt2  = dt

     ! Set the boundary values:
     !
     ! GAR (12/4/2017): This code implements homogeneous Neumann and 
     !                  homogeneous/nonhomogeneous Dirichlet Boundary 
     !                  conditions.
     !
     ! For Homogeneous Dirichlet Boundary Conditions:
     ! - set "read_bcs_from_initfile" flag to false
     !
     ! For Nonhomogeneous Dirichlet Boundary Conditions:
     ! - set "read_bcs_from_initfile" flag to true
     ! - this will set the boundary value to be that of the initial 
     !   fields
     !
     ! For Homogeneous Neumann Boundary conditions:
     ! - set "read_bcs_from_initfile" flag to false
     
     if ( read_bcs_from_initfile ) then
        ux_b  = ux
        uy_b  = uy
        uz_b  = uz
        rho_b = rho
     endif

  endif

end subroutine setup_initial_conditions
