subroutine apply_nullspace_projection( vx, vz )
! Applies the exact nullspace projection to a vector field (vx,vz).
!
!
! NOTES:
!  - GAR (8/22/2017): This method only works for 2D runs. What I
!                     did here is to correct the arrays to properly
!                     account for the dimensions of vx and vz. We
!                     will put an if statement in the main routine
!                     to ensure that when this function gets called
!                     when the run is purely 2D. For a 3D implementation
!                     the math needs to be revised.

  use constants, only:             n, nsuby, nsubz
  use nullspace, only:             null_basis, null_dim
  use precision, only:             dp
  use woodbury_matrices, only:     dimA, numA_per_rank, rpk

  implicit none

  real(kind=dp), dimension(1:rpk,1:nsuby), intent(inout)  :: vx
  real(kind=dp), dimension(1:rpk,1:nsuby), intent(inout)  :: vz

  ! Internal variables
  integer                                         :: a, z, ae, ze, ii, ll, jj, elt_ndx
  real(kind=dp), dimension(1:dimA)                :: iivx, iivz
  real(kind=dp), dimension(1:2*n*n)               :: iiv
  real(kind=dp), dimension(1:null_dim)            :: Nv

  ! Loop over each wavenumber
  do ll = 1, nsuby

     ! Loop over all the elements in this rank, applying the null-space
     ! projection element-wise.
     do ii = 1, numA_per_rank

        ! Set the ownership range of this vertical strip of elements.
        a = (ii - 1) * dimA + 1
        z = (ii - 1) * dimA + dimA

        ! Get the velocities in this vertical strip of elements.
        iivx = vx(a:z, ll)
        iivz = vz(a:z, ll)

        ! Permute the velocity fields into element-wise indexing.
        call permute_xi2e( n, 1, nsubz, iivx )
        call permute_xi2e( n, 1, nsubz, iivz )

        ! Loop over all elements in this vertical strip of elements.
        do jj = 1, nsubz

           ! Set the ownership range of this element.
           ae = (jj - 1) * n * n + 1
           ze = (jj - 1) * n * n + n * n

           ! Set this element index.
           elt_ndx = (ii - 1) * nsubz + jj

           ! Assemble the velocity vector on this element.
           iiv                      = 0.0_dp
           iiv(1:n * n)             = iivx(ae:ze)
           iiv(n * n + 1:2 * n * n) = iivz(ae:ze)

           ! Apply the null-space projection.
           Nv  = 0.0_dp
           Nv  = matmul( transpose( null_basis(:, :, elt_ndx) ), iiv )
           iiv = matmul( null_basis(:, :, elt_ndx), Nv )

           ! Split up the velocity vector back into its components.
           iivx(ae:ze) = iiv(1:n*n)
           iivz(ae:ze) = iiv(n*n+1:2*n*n)

        enddo

        ! Permute back into xi-first indexing and store.
        call permute_e2xi( n, 1, nsubz, iivx )
        call permute_e2xi( n, 1, nsubz, iivz )
        vx(a:z, ll) = iivx
        vz(a:z, ll) = iivz

     enddo

  enddo

end subroutine apply_nullspace_projection
