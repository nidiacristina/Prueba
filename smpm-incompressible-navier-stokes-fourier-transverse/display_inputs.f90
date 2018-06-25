subroutine display_inputs
! This subroutine outputs some information about the simulation being run.

  use constants, only:         n, nprocs, nsg, nsubx, nsuby, nsubz
  !$ use omp_lib, only:        omp_get_max_threads
  use options, only:           fname_init, fname_restart, fname_setup
  use woodbury_matrices, only: dimA, dimblock, k, numA_per_rank, s

  implicit none

  integer           :: number_openmp_threads
  character(len=16) :: caststr

  ! Display paths relevant to this simulation.
  call notify( '' )
  call notify( 'Paths used:' )
  call notify( '' )
  call notify( '   Initial conditions:  ' // fname_init )
  call notify( '   Restart file (read): ' // fname_restart )
  call notify( '   Setup file:          ' // fname_setup )

  ! Display GLL grid parameters.
  call notify( '' )
  call notify( 'Input information:' )
  call notify( '' )

  write( caststr,'(I8.0)' ) n
  call notify( ' n         = ' // caststr )

  write( caststr,'(I8.0)' ) nsubx
  call notify( ' nsubx     = ' // caststr )

  write( caststr,'(I8.0)' ) nsuby
  call notify( ' nsuby     = ' // caststr )

  write( caststr,'(I8.0)' ) nsubz
  call notify( ' nsubz     = ' // caststr )

  write( caststr,'(I8.0)' ) nsg * nsuby
  call notify( ' grid size = ' // caststr )

  call notify( ' ' )
  call notify( 'Decomposition parameters: ' )
  call notify( ' ' )

  write( caststr,'(I8.0)' ) k
  call notify( ' dim. of the Schur problem (S) = ' // caststr )

  write( caststr,'(I8.0)' ) s
  call notify( ' dim. of local interface   (s) = ' // caststr )

  write( caststr,'(I8.0)' ) dimblock * nsuby
  call notify( ' dim. of local grid            = ' // caststr )

  write( caststr,'(I8.0)' ) dimA * nsuby
  call notify( ' dim. of local block       (A) = ' // caststr )

  write( caststr,'(I8.0)' ) numA_per_rank
  call notify( ' num. of local blocks per rank = ' // caststr )

  call notify( '' )
  call notify( 'Parallelization parameters:' )
  call notify( '' )

  write( caststr,'(I8.0)' ) nprocs
  call notify( ' number of MPI processes  = ' // adjustl( caststr ) )

  number_openmp_threads    = 1
  !$ number_openmp_threads = omp_get_max_threads()

  write( caststr, '(I8.0)' ) number_openmp_threads
  call notify( ' number of openMP threads = ' // adjustl( caststr ) )
  call notify( ' ' )

end subroutine display_inputs
