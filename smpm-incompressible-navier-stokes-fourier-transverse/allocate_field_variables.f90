subroutine allocate_solver_variables()
! Allocates and initializes the solver's internal state.

  implicit none

  call allocate_field_variables()
  call allocate_background_current_variables()
  call allocate_deformation_maps()
  call allocate_poisson_matrices()
  call allocate_transverse_variables()

end subroutine allocate_solver_variables

subroutine allocate_field_variables()
! Allocates and initializes grid functions.

  use constants, only:         nsuby
  use field_variables, only:   Cux0, Cux1, Cux2, Cuy0, Cuy1, Cuy2, Cuz0, Cuz1, Cuz2, div_u_int, &
                               Nrho0, Nrho1, Nrho2, &
                               Nux0, Nux1, Nux2, Nuy0, Nuy1, Nuy2, Nuz0, Nuz1, Nuz2, p, px, py, pz, &
                               rho, rho0, rho1, rho2, rho_b, rho_bar, rho_bar_z, rho_bar_zz, rho_int, &
                               ux, ux0, ux1, ux2, ux_b, ux_int, &
                               uy, uy0, uy1, uy2, uy_b, uy_int, &
                               uz, uz0, uz1, uz2, uz_b, uz_int
  use geom, only:              cx, cz
  use precision, only:         dp
  use transverse, only:        cy
  use woodbury_matrices, only: rpk

  implicit none

  ! Allocate the field variables.
  allocate( ux(1:rpk, 1:nsuby), uy( 1:rpk, 1:nsuby ), uz(1:rpk, 1:nsuby), rho(1:rpk, 1:nsuby), p(1:rpk, 1:nsuby) )

  ! Allocate the grid.
  allocate( cx(1:rpk, 1:nsuby), cy( 1:rpk, 1:nsuby ), cz(1:rpk, 1:nsuby) )

  ! Allocate the stratification variables.
  allocate( rho_bar(1:rpk, 1:nsuby), rho_bar_z(1:rpk, 1:nsuby), rho_bar_zz(1:rpk, 1:nsuby) )

  ! Allocate the intermediate variables for the KIO time-splitting.
  allocate( ux0(1:rpk, 1:nsuby), ux1(1:rpk, 1:nsuby), ux2(1:rpk, 1:nsuby) )
  allocate( uy0(1:rpk, 1:nsuby), uy1(1:rpk, 1:nsuby), uy2(1:rpk, 1:nsuby) )
  allocate( uz0(1:rpk, 1:nsuby), uz1(1:rpk, 1:nsuby), uz2(1:rpk, 1:nsuby) )
  allocate( rho0(1:rpk, 1:nsuby), rho1(1:rpk, 1:nsuby), rho2(1:rpk, 1:nsuby) )
  allocate( Nux0(1:rpk, 1:nsuby), Nux1(1:rpk, 1:nsuby), Nux2(1:rpk, 1:nsuby) )
  allocate( Nuy0(1:rpk, 1:nsuby), Nuy1(1:rpk, 1:nsuby), Nuy2(1:rpk, 1:nsuby) )
  allocate( Nuz0(1:rpk, 1:nsuby), Nuz1(1:rpk, 1:nsuby), Nuz2(1:rpk, 1:nsuby) )
  allocate( Cux0(1:rpk, 1:nsuby), Cux1(1:rpk, 1:nsuby), Cux2(1:rpk, 1:nsuby) )
  allocate( Cuy0(1:rpk, 1:nsuby), Cuy1(1:rpk, 1:nsuby), Cuy2(1:rpk, 1:nsuby) )
  allocate( Cuz0(1:rpk, 1:nsuby), Cuz1(1:rpk, 1:nsuby), Cuz2(1:rpk, 1:nsuby) )
  allocate( Nrho0(1:rpk, 1:nsuby), Nrho1(1:rpk, 1:nsuby), Nrho2(1:rpk, 1:nsuby) )
  allocate( px(1:rpk, 1:nsuby), py( 1:rpk, 1:nsuby), pz(1:rpk, 1:nsuby) )
  allocate( ux_int(1:rpk, 1:nsuby), uy_int(1:rpk, 1:nsuby), uz_int(1:rpk, 1:nsuby), rho_int(1:rpk, 1:nsuby) )
  allocate( div_u_int(1:rpk, 1:nsuby) )

  ! Allocate arrays for storing the boundary condition values.
  allocate( ux_b(1:rpk, 1:nsuby), uy_b(1:rpk, 1:nsuby), uz_b(1:rpk, 1:nsuby), rho_b(1:rpk, 1:nsuby) )

  ! Initialize the intermediate time-stepping field variables to zero.
  ux0  = 0.0_dp
  ux1  = 0.0_dp
  ux2  = 0.0_dp
  Nux0 = 0.0_dp
  Nux1 = 0.0_dp
  Nux2 = 0.0_dp
  Cux0 = 0.0_dp
  Cux1 = 0.0_dp
  Cux2 = 0.0_dp

  uy0  = 0.0_dp
  uy1  = 0.0_dp
  uy2  = 0.0_dp
  Nuy0 = 0.0_dp
  Nuy1 = 0.0_dp
  Nuy2 = 0.0_dp
  Cuy0 = 0.0_dp
  Cuy1 = 0.0_dp
  Cuy2 = 0.0_dp

  uz0  = 0.0_dp
  uz1  = 0.0_dp
  uz2  = 0.0_dp
  Nuz0 = 0.0_dp
  Nuz1 = 0.0_dp
  Nuz2 = 0.0_dp
  Cuz0 = 0.0_dp
  Cuz1 = 0.0_dp
  Cuz2 = 0.0_dp

  rho0  = 0.0_dp
  rho1  = 0.0_dp
  rho2  = 0.0_dp
  Nrho0 = 0.0_dp
  Nrho1 = 0.0_dp
  Nrho2 = 0.0_dp

  ux_b  = 0.0_dp
  uy_b  = 0.0_dp
  uz_b  = 0.0_dp
  rho_b = 0.0_dp

end subroutine allocate_field_variables

subroutine allocate_background_current_variables()
! Allocates and initializes background current variables.

  use constants, only:         nsuby
  use field_variables, only:   ubc, dubcdz, dudx, dudy, dudz, &
                               dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, &
                               drhodx, drhody, drhodz
  use precision, only:         dp
  use sponge, only:            raycoeff
  use woodbury_matrices, only: rpk

  implicit none

  ! Allocate background current
  allocate( ubc(1:rpk, 1:nsuby) )
  allocate( dubcdz(1:rpk, 1:nsuby) )
  allocate( raycoeff(1:rpk, 1:nsuby))
  allocate( dudx(1:rpk,1:nsuby), dudy(1:rpk,1:nsuby), dudz(1:rpk,1:nsuby) )
  allocate( dvdx(1:rpk,1:nsuby), dvdy(1:rpk,1:nsuby), dvdz(1:rpk,1:nsuby) )
  allocate( dwdx(1:rpk,1:nsuby), dwdy(1:rpk,1:nsuby), dwdz(1:rpk,1:nsuby) )
  allocate( drhodx(1:rpk,1:nsuby), drhody(1:rpk,1:nsuby), drhodz(1:rpk,1:nsuby) )

  ! Initiliaze background current variables.
  dudx   = 0.0_dp
  dudz   = 0.0_dp
  dwdx   = 0.0_dp
  dwdz   = 0.0_dp
  drhodx = 0.0_dp
  drhodz = 0.0_dp

  dudy   = 0.0_dp
  dvdy   = 0.0_dp
  dwdy   = 0.0_dp
  drhody = 0.0_dp

end subroutine allocate_background_current_variables

subroutine allocate_deformation_maps()
! Allocates and initializes the deformation maps and spectral filtering matrix.

  use constants, only:             n, nky
  use legendre, only:              D, D2, D3, filter_xz, points, wg
  use mesh_deformation_maps, only: detJ, detJ_eta, detJ_xi, d_eta_to_laplacian, d_etaeta_to_laplacian, &
                                   d_xi_to_laplacian, d_xieta_to_laplacian, d_xixi_to_laplacian, &
                                   eta_x, eta_z, x_eta, x_etaeta, x_xi, x_xieta, x_xixi, xi_x, xi_z, &
                                   z_eta, z_etaeta, z_xi, z_xieta, z_xixi, &
                                   x_xi_n, x_eta_n, z_xi_n, z_eta_n, &
                                   x_xi_n_x, x_xi_n_z, x_eta_n_x, x_eta_n_z, &
                                   z_xi_n_x, z_xi_n_z, z_eta_n_x, z_eta_n_z, &
                                   nx, tx, nz, tz, tx_x, tx_z, tz_x, tz_z
  use precision, only:             dp
  use woodbury_matrices, only:     rpk
  use transverse, only:            filter_y

  implicit none

  ! Allocate collocation points, weights, and the filtering matrix.
  allocate( points(n), wg(n), filter_xz(1:n, 1:n) )
  allocate( filter_y(1:nky) )

  ! Allocate the spectral differentiation matrices.
  allocate( D(1:n, 1:n), D2(1:n, 1:n), D3(1:n, 1:n) )

  ! Allocate the mesh deformation maps.
  allocate( x_xi(1:rpk), x_eta(1:rpk), x_xixi(1:rpk), x_etaeta(1:rpk), x_xieta(1:rpk) )
  allocate( z_xi(1:rpk), z_eta(1:rpk), z_xixi(1:rpk), z_etaeta(1:rpk), z_xieta(1:rpk) )
  allocate( detJ(1:rpk), detJ_xi(1:rpk), detJ_eta(1:rpk) )
  allocate( d_xi_to_laplacian(1:rpk), d_eta_to_laplacian(1:rpk) )
  allocate( d_xixi_to_laplacian(1:rpk), d_etaeta_to_laplacian(1:rpk) )
  allocate( d_xieta_to_laplacian(1:rpk) )
  allocate( x_xi_n(1:rpk), x_eta_n(1:rpk), z_xi_n(1:rpk), z_eta_n(1:rpk) )
  allocate( x_xi_n_x(1:rpk), x_xi_n_z(1:rpk), x_eta_n_x(1:rpk), x_eta_n_z(1:rpk) )
  allocate( z_xi_n_x(1:rpk), z_xi_n_z(1:rpk), z_eta_n_x(1:rpk), z_eta_n_z(1:rpk) )

  ! Allocate the forward and backwards Jacobians.
  allocate( eta_x(1:rpk), xi_x(1:rpk), eta_z(1:rpk), xi_z(1:rpk) )

  ! Allocate the normal and tangential vectors used on the boundary.
  allocate( nx(1:rpk, 1:4), nz(1:rpk, 1:4), tx(1:rpk, 1:4), tz(1:rpk, 1:4) )
  allocate( tx_x(1:rpk, 1:4), tz_x(1:rpk, 1:4), tx_z(1:rpk, 1:4), tz_z(1:rpk, 1:4) )

  nx   = 0.0_dp
  nz   = 0.0_dp
  tx   = 0.0_dp
  tz   = 0.0_dp
  tx_x = 0.0_dp
  tz_z = 0.0_dp
  tx_z = 0.0_dp
  tz_z = 0.0_dp

end subroutine allocate_deformation_maps

subroutine allocate_poisson_matrices()
! Allocates and initializes Poisson matrices.

  use constants, only:         nky, nsubx
  use precision, only:         dp
  use woodbury_matrices, only: dimA, numA_per_rank, dimCk, &
                               A_poisson, A_poisson_pivot, C_poisson_block, &
                               C_coarse, C_coarse_pivot, uC, uC_right, uL, uC_coarse, &
                               xC_previous_real, xC_previous_imag, &
                               S_coarse, S_coarse_pivot, S_poisson, T_poisson, &
                               U_poisson, &
                               BJS, BJS_pivot, preC, preC_pivot, &
                               rpk, s

  implicit none

  ! Allocate the local block Poisson matrices.
  allocate( A_poisson( 1:dimA, 1:dimA, 1:numA_per_rank ) )
  allocate( A_poisson_pivot(1:dimA, numA_per_rank) )
  allocate( U_poisson(1:dimA, 1:dimA, 1:numA_per_rank) )
  allocate( T_poisson(1:dimA, 1:dimA, 1:numA_per_rank) )

  ! Allocate the XXX: matrix.
  allocate( C_poisson_block(1:4*s, 1:2*s, 1:numA_per_rank) )

  ! Allocate the Block-Jacobi Schur preconditioner matrices.
  allocate( BJS(1:4*s, 1:4*s, 1:numA_per_rank, 1:nky) )
  allocate( BJS_pivot(1:4*s, 1:numA_per_rank, 1:nky) )

  ! Allocate the Schur matrix.
  allocate( S_poisson(1:4*s, 1:2*s, 1:numA_per_rank, 1:nky) )

  ! Allocate the deflation matrix.
  allocate( S_coarse(1:nsubx - 1, 1:nsubx - 1, nky), S_coarse_pivot(1:nsubx - 1, nky) )

  ! Allocate the preconditioner matrix.
  allocate( preC(1:4*s, 1:4*s, numA_per_rank) )
  allocate( preC_pivot(1:4*s, numA_per_rank) )

  ! Allocate coarse capacitance setup matrices.
  allocate( C_coarse(1:nsubx - 1, 1:nsubx - 1), C_coarse_pivot(1:nsubx - 1)  )

  ! Allocate the right nullity of the coarse matrix.
  allocate( uC_coarse(1:nsubx - 1) )

  ! Allocate the Poisson kernel vectors.
  allocate( uC(1:dimCk), uC_right(1:dimCk), uL(1:rpk) )

  ! Allocate previous solution to the capacitance problem.
  allocate( xC_previous_real(1:dimCk,1:nky), xC_previous_imag(1:dimCk,1:nky) )

  ! Initialize the matrices.
  A_poisson  = 0.0_dp
  BJS        = 0.0_dp
  BJS_pivot  = 0
  preC       = 0.0_dp
  preC_pivot = 0.0_dp

  ! Set the initial guess for solving the capacitance problem.
  xC_previous_real = 0.0_dp
  xC_previous_imag = 0.0_dp

end subroutine allocate_poisson_matrices

subroutine allocate_transverse_variables()
! Allocates variables related to wavenumbers and FFT buffers.

  use, intrinsic :: iso_c_binding

  use constants, only:  nky, nsuby
  use transverse, only: real_buff_FFT, complex_buff_FFT, ky, &
                        p_real_buff_FFT, p_complex_buff_FFT, qy

  implicit none

  include 'fftw3.f03'

  ! Allocate wavenumber array and wavenumber grid.
  allocate( ky(1:nsuby) )
  allocate( qy(1:nky) )

  ! Allocate the FFT buffers.
  p_real_buff_FFT    = fftw_alloc_real(    int( nsuby, kind=C_SIZE_T ) )
  p_complex_buff_FFT = fftw_alloc_complex( int( nky, kind=C_SIZE_T ) )
  call C_F_POINTER( p_real_buff_FFT,       real_buff_FFT, [nsuby] )
  call C_F_POINTER( p_complex_buff_FFT, complex_buff_FFT, [nky] )

end subroutine allocate_transverse_variables
