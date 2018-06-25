subroutine setup_nullspace_projection
! Pre-computes some local deformation terms for helping the computation of the
! CFL number for establishing stability and setting an adaptive time-step.

  use constants, only:             n, nsubz
  use legendre, only:              D
  use nullspace, only:             null_basis, null_dim
  use mesh_deformation_maps, only: eta_x, eta_z, xi_x, xi_z
  use precision, only:             dp
  use woodbury_matrices, only:     dimA, numA_per_rank

  implicit none

  ! Master element grid spacing terms.
  real(kind=dp), allocatable, dimension(:, :) :: div_elt
  real(kind=dp), dimension(1:dimA)            :: iieta_x, iieta_z, iixi_x, iixi_z

  ! Divergence operator variables.
  real(kind=dp), allocatable, dimension(:, :) :: D_eta, D_xi

  ! SVD call variables.
  real(kind=dp), allocatable, dimension(:, :) :: svd_u, svd_v
  real(kind=dp), allocatable, dimension(:)    :: work, svd_s
  integer                                     :: Lwork, info

  ! Internal variables.
  integer                                     :: row, ii, jj, kk, a, z, ae, ze, elt_ndx

  ! Set the dimension of the nullspace.
  null_dim = n * n + 1

  ! Allocate the null space basis on the elements owned by this rank.
  allocate( null_basis(2 * n * n, null_dim, numA_per_rank * nsubz) )
  null_basis = 0.0_dp

  ! Allocate the element divergence matrix.
  allocate( div_elt(1:n*n, 1:2*n*n) )
  div_elt = 0.0_dp

  ! Build the single element eta and xi derivative matrices once.
  allocate( D_eta(1:n * n, 1: n * n) )
  allocate(  D_xi(1:n * n, 1: n * n) )
  D_eta = 0.0_dp
  D_xi  = 0.0_dp
  do ii = 1, n
     a              = (ii - 1) * n + 1
     z              = (ii - 1) * n + n
     D_xi(a:z, a:z) = D
  enddo

  ! Build the kronecker product manually.
  D_eta = 0.0_dp
  do ii = 1, n * n
     call getrow_AkronI( n, D, n, ii, D_eta(ii, :) )
  enddo

  ! Allocate variables for the call to DGESVD.
  allocate( svd_u(1, 1) )
  allocate( svd_v(1:2*n*n, 1:2*n*n) )
  allocate( svd_s(1:n*n) )
  Lwork = 6 * 2 * n * n
  allocate( work(1:Lwork) )

  ! Loop over all the elements, assembling the divergence matrix and computing
  ! its null space basis.
  elt_ndx = 0
  do ii = 1, numA_per_rank

     ! Set the ownership range of this vertical strip of elements.
     a = (ii - 1) * dimA + 1
     z = (ii - 1) * dimA + dimA

     ! Get the metric terms for this vertical strip.
     iieta_x = eta_x(a:z)
     iieta_z = eta_z(a:z)
     iixi_x  =  xi_x(a:z)
     iixi_z  =  xi_z(a:z)

     ! Permute the metric terms to get them into element first indexing.
     call permute_xi2e( n, 1, nsubz, iieta_x )
     call permute_xi2e( n, 1, nsubz, iieta_z )
     call permute_xi2e( n, 1, nsubz,  iixi_x )
     call permute_xi2e( n, 1, nsubz,  iixi_z )

     ! Loop over all elements in this vertical strip of elements.
     do jj = 1, nsubz

        ! Set the ownership range of this element.
        ae = (jj - 1) * n * n + 1
        ze = (jj - 1) * n * n + n * n

        ! Set the on-rank element index.
        elt_ndx = elt_ndx + 1

        ! Build the divergence matrix for this element.
        row     = 1
        div_elt = 0.0_dp
        do kk = ae, ze
           div_elt(row, 1:n * n)             = iieta_x(kk) * D_eta(row, :) + iixi_x(kk) * D_xi(row, :)
           div_elt(row, n * n + 1:2 * n * n) = iieta_z(kk) * D_eta(row, :) + iixi_z(kk) * D_xi(row, :)
           row = row + 1
        enddo

        ! Compute the SVD of this matrix.
        call DGESVD( 'N', 'A', n * n, 2 * n * n, div_elt, n * n, svd_s, svd_u, 1, svd_v, 2 * n * n, work, Lwork, info )
        svd_v = transpose( svd_v )

        ! Store the null vectors.
        null_basis(:, 1:null_dim, elt_ndx) = svd_v(:, 2 * n * n - null_dim + 1:2 * n * n)

     enddo

  enddo

end subroutine setup_nullspace_projection
