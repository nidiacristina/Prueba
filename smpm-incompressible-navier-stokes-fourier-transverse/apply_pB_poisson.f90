subroutine apply_pB_poisson( q, Bq, q_x, q_z, Nrhs )
! Compute the matrix-vector multiply Bq, where B is the boundary matrix in the
! SMW decomposition of the Poisson operator.
!
! (q_x,q_z) - the gradient of q.
! delta     - the viscosity.
! fac       - the multiplicative factor on the penalty coefficient.

  use options, only:               facrobin_PPE
  use constants, only:             n, nsubz, pd, rank, nprocs
  use mesh_deformation_maps, only: eta_x, eta_z, x_xi_n, z_xi_n
  use precision, only:             dp
  use woodbury_matrices, only:     dimA, dimCk, numA_per_rank, rpk, s

  implicit none

  real(kind=dp), dimension(rpk, 1:Nrhs), intent(in)      :: q
  real(kind=dp), dimension(dimCk, 1:Nrhs), intent(inout) :: Bq
  real(kind=dp), dimension(rpk, 1:Nrhs), intent(in)      :: q_x
  real(kind=dp), dimension(rpk, 1:Nrhs), intent(in)      :: q_z
  integer, intent(in)                                    :: Nrhs

  real(kind=dp)                                          :: fac, delta

  real(kind=dp), dimension(1:3 * n * nsubz, 1:Nrhs)      :: fluxL, fluxR

  integer                                                :: aL, zL     ! start, stop, and gap index.
  integer                                                :: aR, zR     ! start, stop, and gap index.
  real(kind=dp)                                          :: tau, omega, kappa, odb
  integer                                                :: ii, kk
  integer                                                :: L_ndx_start, L_ndx_end, R_ndx_start, R_ndx_end
  integer                                                :: ndx_start, ndx_end
  integer                                                :: ax, zx, az, zz, a, z

  real(kind=dp), parameter                               :: alpha = 1.0_dp
  real(kind=dp), parameter                               :: beta = 1.0_dp

  ! Set some constants.
  delta = 1.0_dp
  fac   = facrobin_PPE
  omega = 2.0_dp / (pd * (pd + 1.0_dp))
  kappa = omega * alpha / beta
  odb   = omega * delta * beta
  tau   = (1.0_dp / odb) * (delta + (2.0_dp * kappa) - (2 * sqrt( kappa**2 + (delta * kappa) )))

  ! Set the output vector to zero (there is no good reason it should ever not
  ! be zero preinitialized to zero).
  Bq = 0.0_dp

  ! XXX: Now that Bq is only dimension dimCk, the L_ndx_start, etc. will have
  !      to change.  Appeal to the sparsity pattern to figure out what they
  !      will have to be.

  ! Loop over the right hand sides.
  !$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE), SHARED( Bq, eta_z, eta_x, x_xi_n, z_xi_n, q, q_x, q_z ), &
  !$OMP PRIVATE( aL, zL, aR, zR, L_ndx_start, L_ndx_end, R_ndx_start, R_ndx_end )
  do kk = 1, Nrhs

     ! Loop over the internal vertical interfaces.
     do ii = 1, numA_per_rank - 1

        ! Get the indices of vertical line of grid points on either side of
        ! this interface (within this rank).
        aL   = (n * n * nsubz) * ii - (n * nsubz - 1)
        zL   = (n * n * nsubz) * ii - 0
        aR   = zL + 1
        zR   = zL + n * nsubz

        ! NOTE: "L" means left side of the interface, "R" means right side of
        !       the interface.

        ! Get the start/stop row indices in the pB_poisson matrix of the
        ! left/right side interfaces.
        if ( rank == 0 ) then
           L_ndx_start = 2 * s * (ii - 1) + 1
           L_ndx_end   = 2 * s * (ii - 1) + s
        else
           L_ndx_start = s + 2 * s * (ii - 1) + 1
           L_ndx_end   = s + 2 * s * (ii - 1) + s
        endif
        R_ndx_start = L_ndx_end + 1
        R_ndx_end   = L_ndx_end + s

        ! Impose the C0 continuity condition on each side of the interface.
        Bq(L_ndx_start:L_ndx_end, kk) = Bq(L_ndx_start:L_ndx_end, kk) - &
                                        fac * tau * hypot(eta_x(aL:zL), eta_z(aL:zL)) * (-1.0_dp * q(aR:zR, kk))
        Bq(R_ndx_start:R_ndx_end, kk) = Bq(R_ndx_start:R_ndx_end, kk) - &
                                        fac * tau * hypot(eta_x(aR:zR) , eta_z(aR:zR)) * (-1.0_dp * q(aL:zL, kk))

        ! Impose the C1 continuity condition on each side of the interface.
        Bq(L_ndx_start:L_ndx_end, kk) = Bq(L_ndx_start:L_ndx_end, kk) - fac * delta * tau * hypot(eta_x(aL:zL), eta_z(aL:zL)) * &
                                        (-1.0_dp * ( z_xi_n(aL:zL) * q_x(aR:zR, kk) - x_xi_n(aL:zL) * q_z(aR:zR, kk) ) )
                  
        Bq(R_ndx_start:R_ndx_end, kk) = Bq(R_ndx_start:R_ndx_end, kk) - fac * delta * tau * hypot(eta_x(aR:zR) , eta_z(aR:zR)) * &
                                        (-1.0_dp * (-z_xi_n(aR:zR) * q_x(aL:zL, kk) + x_xi_n(aR:zR) * q_z(aL:zL, kk)) )

     enddo

  enddo
  !$OMP END PARALLEL DO

  ! Apply the B matrix to the inter-rank interfaces.

  ! Establish the range of the fluxes we'll need to send.
  aL = 1
  zL = n * nsubz
  aR = numA_per_rank * dimA - (n * nsubz - 1)
  zR = numA_per_rank * dimA

  ! Pack up the fluxes to send to the left.
  a                  = 1
  z                  = n * nsubz
  fluxL(a:z, 1:Nrhs) = q_x(aL:zL, 1:Nrhs)

  a                  = n * nsubz + 1
  z                  = 2 * n * nsubz
  fluxL(a:z, 1:Nrhs) = q_z(aL:zL, 1:Nrhs)

  a                  = 2 * n * nsubz + 1
  z                  = 3 * n * nsubz
  fluxL(a:z, 1:Nrhs) = q(aL:zL , 1:Nrhs)

  ! Pack up the fluxes to send to the right.
  a                  = 1
  z                  = n * nsubz
  fluxR(a:z, 1:Nrhs) = q_x(aR:zR, 1:Nrhs)

  a                  = n * nsubz + 1
  z                  = 2 * n * nsubz
  fluxR(a:z, 1:Nrhs) = q_z(aR:zR, 1:Nrhs)

  a                  = 2 * n * nsubz + 1
  z                  = 3 * n * nsubz
  fluxR(a:z, 1:Nrhs) = q(aR:zR, 1:Nrhs)

  ! Communicate with the neighboring ranks.
  call sync_ranks( fluxL, fluxR, 3 * n * nsubz * Nrhs )

  ! Apply the fluxes to the left/right that we just received.

  ! Fluxes on the left side of this rank (the right side of the interface).
  if ( rank > 0 ) then

     ! Set the ownership range of this side.
     ndx_start   = 1
     ndx_end     = s
     aR          = 1
     zR          = n * nsubz

     ! Loop over the right hand sides.
     do kk = 1, Nrhs

        ! Apply C0 continuity.
        a = 2 * n * nsubz + 1
        z = 3 * n * nsubz
        Bq(ndx_start:ndx_end, kk) = Bq(ndx_start:ndx_end, kk) - &
                                    fac * tau * hypot(eta_x(aR:zR) , eta_z(aR:zR)) * &
                                    (-1.0_dp * fluxL(a:z, kk))

        ! Apply C1 continuity.
        ax = 1
        zx = n * nsubz
        az = n * nsubz + 1
        zz = 2 * n * nsubz
        Bq(ndx_start:ndx_end, kk) = Bq(ndx_start:ndx_end, kk) - &
                                    fac * delta * tau * hypot(eta_x(aR:zR) , eta_z(aR:zR)) * &
                                    (-1.0_dp * (-z_xi_n(aR:zR) * fluxL(ax:zx, kk) + x_xi_n(aR:zR) * fluxL(az:zz, kk)) ) 

     enddo

  endif

  ! Fluxes to the right side of this rank ( the left side of this interface ).
  if ( rank < nprocs - 1 ) then

     ! Set the ownership range of this side.
     ndx_start = dimCk - s + 1
     ndx_end   = dimCk
     aL        = numA_per_rank * dimA - (n * nsubz - 1)
     zL        = numA_per_rank * dimA

     ! Loop over the right hand sides.
     do kk = 1, Nrhs

        ! Apply C0 continuity.
        a = 2 * n * nsubz + 1
        z = 3 * n * nsubz
        Bq(ndx_start:ndx_end, kk) = Bq(ndx_start:ndx_end, kk) - &
                                    fac * tau * hypot(eta_x(aL:zL) , eta_z(aL:zL)) * (-1.0_dp * fluxR(a:z, kk))

        ! Apply C1 continuity.
        ax = 1
        zx = n * nsubz
        az = n * nsubz + 1
        zz = 2 * n * nsubz
        Bq(ndx_start:ndx_end, kk) = Bq(ndx_start:ndx_end, kk) - &
                                    fac * delta * tau * hypot(eta_x(aL:zL) , eta_z(aL:zL)) * &
                                    (-1.0_dp * (z_xi_n(aL:zL) * fluxR(ax:zx, kk) - x_xi_n(aL:zL) * fluxR(az:zz, kk)) )

     enddo

  endif

end subroutine apply_pB_poisson
