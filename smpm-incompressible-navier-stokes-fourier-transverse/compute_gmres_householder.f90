subroutine compute_gmres_householder( x, b, n, tol, maxit, restart, apply_matrix )
! A matrix-free implementation of GMRES with Householder reflections via the
! algorithm given in Walker 1988.  This assumes that A is a square matrix, and
! solves via GMRES a linear system ||Ax - b|| < tol*||b||.
!
!	Inputs:
!
!		x 	          - [nx1]: Initial guess.
!		b 	          - [nx1]: Right hand side.
!		n 	          - [1x1]: The length of x and b, and the dimension of A.
!     rankdims     - [1xnprocs]: the number of grid points per rank.
!		tol	       - [1x1]: The GMRES stop tolerance (see above).
!		maxit	       - [1x1]: The maximum total number of iterations to compute.
!     restart      - [1x1]: The number of iterations after which you should restart ( restart << maxit ).
!     apply_matrix - [function]: apply_matrix( Ax, x ) returns y = A*x.
!                                The solution vector, Ax, does not have
!                                to be initialized prior to calling apply_matrix().
!
!	Outputs:
!
!		x		- [nx1]: The computed GMRES solution.
!		tol	- [1x1]: The computed GMRES residual.
!		maxit	- [1x1]: The number of iterations computed.
!
! 27 Oct 2013
! Sumedh Joshi
! Cornell University

  use constants, only:               rank, root, nprocs
  use mpi, only:                     MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, &
                                     MPI_INTEGER, MPI_IN_PLACE, MPI_SUM
  use parallel_linear_algebra, only: assign_value_at, pnorm2
  use precision, only:               dp

  implicit none

  ! Interface variable declarations.
  real(kind=dp), dimension(1:n), intent(inout)   :: x
  real(kind=dp), dimension(1:n), intent(in)      :: b
  integer, intent(in)                            :: n   ! The dimension of the block on this rank.
  real(kind=dp), intent(inout)                   :: tol
  integer, intent(inout)                         :: maxit
  integer, intent(inout)                         :: restart

  interface
     subroutine apply_matrix( Ax, x )

        use precision, only: dp

        real(kind=dp), dimension(*), intent(in)  :: x
        real(kind=dp), dimension(*), intent(out) :: Ax

     end subroutine apply_matrix
  end interface

  ! Internal variables.
  integer, dimension(1:nprocs)                   :: rankdims
  real(kind=dp), allocatable, dimension(:, :)    :: P, G, R
  integer                                        :: ii, jj, m, k, householder_ndx
  real(kind=dp), allocatable, dimension(:)       :: w, r0, x0, Ax, v, smallP, ym, ek, iiv
  real(kind=dp), allocatable, dimension(:)       :: wk, vk, ekk, tmp, v_gather, w_gather
  real(kind=dp)                                  :: eps
  real(kind=dp), dimension(1:2, 1:2)             :: givensjj
  integer                                        :: ierr
  integer                                        :: own_start, own_end
  integer, dimension(1:nprocs)                   :: recv_counts, displacements
  real(kind=dp)                                  :: are_any_nonzero
  integer                                        :: local_ndx
  integer                                        :: iteration_count
  real(kind=dp)                                  :: norm_rhs
  integer                                        :: is_converged

  ! Set machine precision.
  eps = 1e-16_dp

  ! Set the norm of the right-hand-side for relative residual computation.  If
  ! the norm is less than, set it to one.
  norm_rhs = pnorm2( b )
  if ( norm_rhs > 0.0_dp .and. norm_rhs < 1.0_dp ) then
     norm_rhs = 1.0_dp
  endif

  ! Check to see if RHS is zero.
  if ( norm_rhs < eps ) then
     x         = 0.0_dp
     maxit     = 0
     tol       = 0.0_dp
     return
  endif

  ! Get the dimension of all of the ranks.
  rankdims           = 0
  rankdims(rank + 1) = n
  call MPI_ALLGATHER( n, 1, MPI_INTEGER, &
                      rankdims, 1, MPI_INTEGER, &
                      MPI_COMM_WORLD, ierr )

  ! Make sure both the restart and the maximum number of iterations are less
  ! than the problem size for every rank.
  if ( maxit > minval( rankdims ) ) then
     maxit = minval( rankdims )
  endif

  ! Disable restart feature.
  restart = maxit

  ! Get the ownership range of this sub-block.
  own_start = 1
  do ii = 1, rank
     own_start = own_start + rankdims(ii)
  enddo
  own_end = own_start + n  - 1

  ! We'll need the displacements and receive counts for an MPI all-gather
  ! operation.
  do ii = 1, nprocs
     recv_counts(ii)   = rankdims(ii)
  enddo
  displacements(1) = 0
  do ii = 2, nprocs
     displacements(ii) = displacements(ii - 1) + recv_counts(ii - 1)
  enddo

  ! Allocate the Krylov basis vectors and the residual vectors on each process.
  allocate( P(1:n, 1:restart) )

  ! Allocate some variables.
  !
  ! NOTE: This is horrendously dumb.  I'm storing tons more data than I need
  !       to.  I need to find a better way to do dynamically size arrays in
  !       FORTRAN.  All of the arrays that are allocated from 1:maxit suffer
  !       from this problem.  This absolutely has to be fixed.
  allocate( Ax(1:n), R(1:restart, 1:restart), ek(1:n), G(1:2, 1:restart), r0(1:n) )
  allocate( smallP(1:n), v(1:n), w(1:n), x0(1:n), ym(1:restart), iiv(1:n) )
  allocate( v_gather( 1:sum(rankdims) ), w_gather( 1:sum(rankdims) ) )

  ! Allocate rank-level blocks of any arrays needed.
  allocate( wk(1:n), vk(1:n), ekk(1:n), tmp(1:n) )

  ! 1. Compute the residual for the initial guess, and find the appropriate
  !    Householder reflector.
  x0 = x;
  call apply_matrix( Ax, x0 )
  r0 = b - Ax
  v  = 0.0_dp
  if ( rank == 0 ) then
     householder_ndx = 1
  else
     householder_ndx = -1
  endif
  call get_householder_parallel( n, v, r0, householder_ndx )
  P(1:n, 1) = v
  w         = r0
  call apply_householder_parallel( n, w, P(1:n, 1) )
  w_gather( own_start:own_end ) = w
  call MPI_ALLGATHERV( w,        n, MPI_DOUBLE_PRECISION, &
                       w_gather, recv_counts, displacements, MPI_DOUBLE_PRECISION, &
                       MPI_COMM_WORLD, ierr )

  ! Commence the GMRES loop.  The comments parallel the 1988 Walker paper.

  ! Start the iteration counter.
  iteration_count = 1

  ! 2. For m=1,2,...,maxit do:
  m = 1
  do

     ! a) Evaluate the v = PmPm-1...P1AP1...Pmem
     v    = 0.0_dp
     v = assign_value_at( v, own_start, own_end, m, 1.0_dp )
     do jj = m, 1, -1
        call apply_householder_parallel( n, v, P(1:n, jj) )
     enddo
     Ax = 0.0_dp
     call apply_matrix( Ax, v )
     v  = Ax
     do jj = 1, m
        call apply_householder_parallel( n, v, P(1:n, jj) )
     enddo

     ! Check to see if all v(m+1) = ... = v(n) = 0.
     are_any_nonzero = 0.0_dp
     if ( m + 1 .le. own_end ) then
        local_ndx = max( m + 1 - own_start + 1, 1 )
        do ii = local_ndx, n
           if ( abs( v(ii) ) .ge. eps ) then
              are_any_nonzero = are_any_nonzero + 1.0_dp
           endif
        enddo
     endif
     call MPI_ALLREDUCE( MPI_IN_PLACE, are_any_nonzero, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )

     ! b) If v(m+1) = ... = v(n) = 0, then proceed to step (e).  Otherwise continue.
     !if ( are_any_nonzero > 0 .and. m < maxit ) then
     if ( are_any_nonzero > 0 .and. m < restart ) then

        ! c) Determine Pm+1 with a Householder vector having first m components
        !    0, such that Pm+1v has zero components after the (m+1)st.
        iiv = 0.0_dp
        if ( m + 1 .le. own_end ) then     ! If I find myself doing this again,
                                           ! write a function assign_value_inrange
           local_ndx        = max( m + 1 - own_start + 1, 1 )
           iiv(local_ndx:n) = v(local_ndx:n)
        endif
        if ( m + 1 .le. own_end .and. m + 1 .ge. own_start ) then
           householder_ndx = m + 1 - own_start + 1
        else
           householder_ndx = -1
        endif
        call get_householder_parallel( n, smallP, iiv, householder_ndx )
        P(1:n, m + 1) = smallP

        ! d) Overwrite v <--- Pm+1v.
        call apply_householder_parallel( n, v, P(1:n, m + 1) )

     endif

     ! Grab the iterative solution on all ranks for Givens rotation
     ! application (which is inherently serial anyways).
     !
     ! XXX: this can be optimized by only passing one entry to the left and
     !      the right, and then only the ranks which have ownership of entries
     !      1,...,m will do the rotations.  Profile with Greg to see if this
     !      optimization is worth it.
     !call MPI_ALLGATHERV( v,        n, MPI_DOUBLE_PRECISION, &
     !                     v_gather, recv_counts,  displacements, MPI_DOUBLE_PRECISION, &
     !                     MPI_COMM_WORLD, ierr )

     ! If you can assume that m < n for all ranks' n, then you can do this a
     ! lot faster.  If you can't you can do some modular arithmetic to see
     ! what each rank has to do.
     if ( rank == root ) then

       ! e) If m > 1, overwrite v <--- Jm-1 ... J1v
       if ( m > 1 ) then
          do jj = 1, m - 1

             ! Construct the Givens rotation.
             givensjj(1, 1) =  G(1, jj)
             givensjj(1, 2) =  G(2, jj)
             givensjj(2, 1) = -G(2, jj)
             givensjj(2, 2) =  G(1, jj)

             ! Apply this Givens rotation matrix.
             v(jj:jj+1) = matmul( givensjj, v(jj:jj+1) )

          enddo
       endif

     ! f) If v(m+1) = 0, proceed to step (i), otherwise continue.
     !if ( .not. (abs( v_gather( m + 1 ) ) < eps) ) then
     if ( v( m + 1 ) .ne. 0  ) then

        ! g) Determine Jm acting on components m and m+1 such that (Jmv)(m+1) is zero.
        G(1, m) = v(m)   / sqrt( v(m)**2 + v(m+1)**2 )
        G(2, m) = v(m+1) / sqrt( v(m)**2 + v(m+1)**2 )

        ! h) Overwrite v <--- Jmv and w <--- Jmw.
        givensjj(1, 1) =  G(1, m)
        givensjj(1, 2) =  G(2, m)
        givensjj(2, 1) = -G(2, m)
        givensjj(2, 2) =  G(1, m)

        v(m:m+1) = matmul( givensjj, v(m:m+1) )
        w(m:m+1) = matmul( givensjj, w(m:m+1) )

        !v = v_gather( own_start:own_end )
        !w = w_gather( own_start:own_end )

     endif
     !v = v_gather( own_start:own_end )
     !w = w_gather( own_start:own_end )

     ! i) Set the matrix Rm = ... etc.
     R(1:restart, m) = v(1:restart)

     endif

     ! All ranks check to see if root is reporting convergence.
     is_converged = 0
     if ( rank == root .and. abs( w(m+1) ) < tol * norm_rhs ) then
        is_converged = 1
     endif
     call MPI_BCAST( is_converged, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr )
     ! j) If |w(m+1)| < tol (we've reached convergence) or m = maxit
     !    (we failed to converge within the requested time), then solve for ym
     !    and overwrite x0 with xm.  otherwise increment m.
     !if ( abs( w_gather( m + 1 ) ) < tol * norm_rhs .or. m == restart .or. iteration_count .ge. maxit ) then
     if ( is_converged  == 1 .or. m == restart-1 .or. iteration_count .ge. maxit ) then

        ! 3. Solve for ym and overwrite x0 <--- xm.

        ! a) Determine ym which minimizes the least-squares problem for
        !    ||w - Rm y|| by solving the m x m upper triangular system
        !    R(1:m,1:m)ym = w(1:m).  This is just a back-substitution
        !    operation.
        ! Root does the back-substitution.
        ym = 0.0_dp
        if ( rank == root ) then
           ym(m) = w(m) / R(m, m)

           do ii = m-1, 1, -1
              ym(ii) = (w(ii) - dot_product( R(ii, ii+1:m), ym(ii+1:m) )) / R(ii, ii)
           enddo
        endif

        ! Broadcast the result to all ranks.
        call MPI_BCAST( ym, restart, &
                        MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr )

        ! b) For k = 1,...,m do:
        !    overwrite x0 <--- x0 + ym(k)(P1...Pk)ek.
        ! Every rank participates in Householder vector application.
        do k = 1, m

           ekk = 0.0_dp
           ekk = assign_value_at( ekk, own_start, own_end, k, 1.0_dp )
           do jj = k, 1, -1
              call apply_householder_parallel( n, ekk, P(1:n, jj) )
           enddo
           x0 = x0 + ym(k) * ekk

        enddo

        ! Update our solution, the current state of convergence, and break
        ! out of the loop.
        !
        ! NOTE: The original algorithm differentiated between exceeding
        !       the iteration threshold and converging, by complaining on
        !       standard output in the former, though this is not necessary in
        !       the SMPM code base.  With that, both cases are treated
        !       identically.

        ! If we've maxed out the total number of iterations or we've solved the
        ! problem, exit.  If not, restart.
        if ( iteration_count .ge. maxit .or. is_converged == 1 ) then

           x     = x0
           maxit = iteration_count
           tol   = abs( w(m+1) ) / norm_rhs
           call MPI_BCAST( w(m+1), 1, &
                           MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr )
           exit

        else

           ! Set the current solution as the initial guess.
           x = x0
           m = 1

           ! Repeat the initialization step of GMRES.
           ! 1. Compute the residual for the initial guess, and find the
           !    appropriate Householder reflector.
           x0 = x
           call apply_matrix( Ax, x0 )
           r0 = b - Ax
           v  = 0.0_dp
           if ( rank == 0 ) then
              householder_ndx = 1
           else
              householder_ndx = -1
           endif
           call get_householder_parallel( n, v, r0, householder_ndx )
           P(1:n, 1) = v
           w         = r0
           call apply_householder_parallel( n, w, P(1:n, 1) )
           w_gather( own_start:own_end ) = w
           call MPI_ALLGATHERV( w,        n, MPI_DOUBLE_PRECISION, &
                                w_gather, recv_counts, displacements, MPI_DOUBLE_PRECISION, &
                                MPI_COMM_WORLD, ierr )
        endif

     else

        ! Increment the counter and go on to the next iteration of GMRES.
        m = m + 1

        ! Increment the total iteration count.
        iteration_count = iteration_count + 1

     endif

  enddo

end subroutine compute_gmres_householder

! Ancillary subroutines used by compute_gmres_householder follow.

subroutine get_householder_parallel( n, v, x, ii )
! Get the Householder vector v for the transformation that takes x --> ||x|| e_ii

  use parallel_linear_algebra, only: pnorm2
  use precision, only:               dp

  implicit none

  integer, intent(in)                        :: n, ii ! Note that ii is relative to the current rank.  If the desired
  real(kind=dp), dimension(1:n), intent(in)  :: x     ! index falls outside of the current rank's ownership range, pass in
  real(kind=dp), dimension(1:n), intent(out) :: v     ! ii = -1.
  real(kind=dp)                              :: normx, normv

   ! Compute the Householder vector.
   normx = pnorm2( x )
   v     = x

   ! If the householder index lies in this rank, apply the reflector.
   if ( ii .ne. -1 ) then
      v(ii) = v(ii) - normx
   endif
   normv = pnorm2( v )
   v     = v / normv

end subroutine get_householder_parallel

subroutine apply_householder_parallel( n, x, v )
! Apply the Householder transformation: x <-- x - 2v(v'x)

  use parallel_linear_algebra, only: pdot_product
  use precision, only:               dp

  implicit none

  integer, intent(in)                          :: n
  real(kind=dp), dimension(1:n), intent(inout) :: x
  real(kind=dp), dimension(1:n), intent(in)    :: v

  ! Compute the dot product in parallel.
  x = x - 2 * v * pdot_product( x, v )

end subroutine apply_householder_parallel
