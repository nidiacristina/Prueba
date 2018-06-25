subroutine apply_filter_xz( rho, F )
! Applies the filtering matrix computed by setup_filter to a target vector
! rho, in the x-z plane.  See Blackburn, JCP 186 (2003) 610 - 629 for details.

  use constants, only:         n, nsuby, nsubz
  use precision, only:         dp
  use woodbury_matrices, only: numA_per_rank, rpk

  implicit none

  real(kind=dp), dimension(1:rpk, 1:nsuby), intent(inout) :: rho
  real(kind=dp), dimension(n, n), intent(in)              :: F

  ! Local variables.
  integer                                                 :: numsub
  integer                                                 :: ns
  real(kind=dp), dimension(1:rpk)                         :: rhoN
  real(kind=dp), dimension(n, n)                          :: rhoT, rhoF1, rhoF2
  integer                                                 :: ii, i, j, k, c

  ! Set some constants.
  ns     = n * n
  numsub = nsubz * numA_per_rank

  ! Loop over each wavenumber filtering each element seperately.
  do ii = 1, nsuby

     ! Permute into element-first indexing.
     call permute_z( n, numA_per_rank, nsubz, rho(:, ii), .true. )

     rhoN = 0.0_dp
     do k = 0, numsub - 1
        rhoT  = 0.0_dp
        rhoF1 = 0.0_dp

        ! Extracting the subdomain information from global solution.
        do i = 0, n-1
           do j = 1, n
              c            = (ns * k) + (i * n) + j
              rhoT(i+1, j) = rho(c, ii)
           enddo
        enddo

        ! Now, filter the solution at the subdomain.
        rhoF1 = matmul( F, rhoT )
        rhoF2 = matmul( rhoF1, transpose( F ) )

        ! Returning the filtered solution to the global system.
        do i = 0, n - 1
           do j = 1, n
              c       = (ns * k) + (i * n) + j
              rhoN(c) = rhoN(c) + rhoF2(i + 1, j)
           enddo
        enddo
     enddo

     rho(:, ii) = rhoN

     ! Permute into xi-first indexing.
     call permute_z( n, numA_per_rank, nsubz, rho(:, ii), .false. )

  enddo

end subroutine apply_filter_xz
