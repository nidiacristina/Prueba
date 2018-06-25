subroutine enforce_velocity_bc( px, py, pz )
! Explicitly constrains the pressure gradient to conform to the boundary
! conditions of the vector viscous equation.  This is done to avoid creating
! bounary-mismatch in the vector viscous equations.

  use constants, only:             n, nprocs, nsuby, nsubz, rank
  use mesh_deformation_maps, only: nx, nz
  use precision, only:             dp
  use woodbury_matrices, only:     rpk, s
  use options, only:               bc_flag_viscous

  implicit none

  real(kind=dp), dimension(1:rpk,1:nsuby), intent( inout ) :: px
  real(kind=dp), dimension(1:rpk,1:nsuby), intent( inout ) :: py
  real(kind=dp), dimension(1:rpk,1:nsuby), intent( inout ) :: pz

  integer                                                  :: a, kk, z, gap

  ! Basically the idea is...
  !               if its a no-slip boundary,   make grad(p)     = 0.0_dp
  !               if its a free-slip boundary, make <grad(p),n> = px * nx + pz * nz
  !                                                             = 0 <==> px = nz/nx * pz
  ! Sumedh Joshi

  ! Note: While the no-slip approach can be implemented in 3D, the free-slip,
  !       cannot because there is no bounding surface in the transverse.
  !       Therefore, the idea is modified as follows:
  !               if its a no-slip boundary,   make grad(p)     = 0.0_dp
  !               if its a free-slip boundary, make <grad(p),n> = px * nx + py * 0 + pz * nz
  !                                                             = 0 <==> px = -nz/nx * pz
  ! 20 June 2017
  ! Gustavo Rivera w/ Sumedh Joshi

  ! Transverse Loop
  do kk = 1,nsuby

     ! Bottom boundary.
     if ( bc_flag_viscous(1) == 1 ) then
        a               = 1
        z               = rpk
        gap             = n * nsubz

        px(a:z:gap,kk)  = 0.0_dp
        py(a:z:gap,kk)  = 0.0_dp
        pz(a:z:gap,kk)  = 0.0_dp
     else
        a               = 1
        z               = rpk
        gap             = n * nsubz

        ! This is the bottom boundary, so nz should never be zero.  Divide by nz.
        pz(a:z:gap, kk) = -nx(a:z:gap, 1) * px(a:z:gap,kk) / nz(a:z:gap, 1)
     endif

     ! Right boundary.
     if ( rank == nprocs - 1 ) then
        if ( bc_flag_viscous(2) == 1 ) then
           a            = rpk - s + 1
           z            = rpk

           px(a:z,kk)   = 0.0_dp
           py(a:z,kk)   = 0.0_dp
           pz(a:z,kk)   = 0.0_dp
        else
           a            = rpk - s + 1
           z            = rpk

           ! This is the right boundary, so nx should never be zero. Divide by nx.
           px(a:z,kk)   = -nz(a:z, 2) * pz(a:z,kk) / nx(a:z, 2)
        endif
     endif

     ! Top boundary.
     if ( bc_flag_viscous(3) == 1 ) then
        a               = s
        z               = rpk
        gap             = s

        px(a:z:gap,kk)  = 0.0_dp
        py(a:z:gap,kk)  = 0.0_dp
        pz(a:z:gap,kk)  = 0.0_dp
     else
        a               = s
        z               = rpk
        gap             = s

        ! This is the top boundary, so nz should never be zero.  Divide by nz.
        pz(a:z:gap,kk)  = -nx(a:z:gap, 3) * px(a:z:gap,kk) / nz(a:z:gap, 3)
     endif

     ! Left boundary.
     if ( rank == 0 ) then
        if ( bc_flag_viscous(4) == 1 ) then
           a            = 1
           z            = s

           px(a:z,kk)   = 0.0_dp
           px(a:z,kk)   = 0.0_dp
           pz(a:z,kk)   = 0.0_dp
        else
           a            = 1
           z            = s

           ! This is the left boundary, so nx should never be zero. Divide by nx.
           px(a:z,kk)   = -nz(a:z, 4) * pz(a:z,kk) / nx(a:z, 4)
        endif
     endif

  enddo !End Transverse Loop
end subroutine enforce_velocity_bc
