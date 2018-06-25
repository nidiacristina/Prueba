subroutine read_setupfile( fname_setup, success_flag )
! Read a setup file, which contains all the assembled matrices and
! kernel vectors the SMPM code needs to run.

  use constants, only:         nky, nprocs, nsg, nsubx, rank
  use precision, only:         int64
  use woodbury_matrices, only: arank, arankC, BJS, BJS_pivot, dimA, dimCk, k, numA_per_rank, &
                               s, S_coarse, S_coarse_pivot, S_poisson, rpk, &
                               T_poisson, uC, uC_coarse, uL, U_poisson

  implicit none

  character(len=*), intent(in)   :: fname_setup
  logical, intent(out)           :: success_flag

  integer(kind=int64)            :: byte_offset, bytes_per_real, bytes_per_cplx, ii, jj, iioffset, counter, local_block_size
  integer(kind=int64)            :: uC_size, uL_size, uS_coarse_size, U_poisson_size, T_poisson_size, S_poisson_size, &
                                    BJS_size, BJS_pivot_size, S_coarse_size, S_coarse_pivot_size
  integer(kind=int64)            :: uC_units, uL_units, uS_coarse_units, U_poisson_units, T_poisson_units, S_poisson_units, &
                                    BJS_units, BJS_pivot_units, S_coarse_units, S_coarse_pivot_units
  integer(kind=int64)            :: uC_offset, uL_offset, uS_coarse_offset, U_poisson_offset, T_poisson_offset, S_poisson_offset, &
                                    BJS_offset, S_coarse_offset

  integer                        :: open_status

  ! The new setup file needs to store:
  !
  ! U_poisson:   dimA**2 * numA_per_rank (complex)
  ! T_poisson:   dimA**2 * numA_per_rank (complex)
  ! S_poisson:         8 * s**2 * numA_per_rank * nky
  ! BJS:       16 * s**2 * numA_per_rank * nky
  ! BJS_pivot:  4 * s * numA_per_rank * nky
  ! S_coarse:  ( nsubx - 1 )**2 * nky
  ! S_coarse_pivot: ( nsubx - 1 ) * nky

  ! Open the file for reading.
  open( 14850, access='stream', form='unformatted', status='old', action='read', &
        iostat=open_status, file=trim( fname_setup ) )

  call sync_flag( open_status )

  ! We assume that since we were able to open the file, everything is fine.
  !
  ! NOTE: We don't actually validate what we've read, nor do we validate
  !       each read statement.  Realistically there is no recovery from
  !       a failed read nor is the effort to detect it worth the effort,
  !       especially since this will be replaced with HDF5-based I/O in
  !       the future.
  success_flag = (open_status == 0)
  if ( success_flag .eqv. .false. ) then
     close( 14850 )
     return
  end if

  ! Set some constants related to word size.
  bytes_per_real =  8
  bytes_per_cplx = 16

  ! Set some constants related to block sizes per rank for each of the
  ! constitutent pieces.  Note the order below is the storage order.
  uC_units             = dimCk
  uL_units             = rpk
  U_poisson_units      = dimA**2 * numA_per_rank
  T_poisson_units      = dimA**2 * numA_per_rank
  S_poisson_units      =  8 * s**2 * numA_per_rank * nky
  BJS_units            = 16 * s**2 * numA_per_rank * nky
  BJS_pivot_units      =  4 * s    * numA_per_rank * nky
  S_coarse_units       = (nsubx - 1)**2 * nky
  S_coarse_pivot_units = (nsubx - 1)    * nky
  uS_coarse_units      = (nsubx - 1)

  ! Set the global sizes in bytes of each of the logical units being written
  ! to disk.
  uC_size              = k                        * bytes_per_real
  uL_size              = nsg                      * bytes_per_real
  U_poisson_size       = U_poisson_units * nprocs * bytes_per_cplx
  T_poisson_size       = T_poisson_units * nprocs * bytes_per_cplx
  S_poisson_size       = S_poisson_units * nprocs * bytes_per_real
  BJS_size             = BJS_units       * nprocs * bytes_per_real
  BJS_pivot_size       = BJS_pivot_units * nprocs * bytes_per_real
  S_coarse_size        = S_coarse_units           * bytes_per_real
  S_coarse_pivot_size  = S_coarse_pivot_units     * bytes_per_real
  uS_coarse_size       = uS_coarse_units          * bytes_per_real

  ! Set the offsets for each rank; this is the byte location that this rank
  ! will start writing its piece of each variable.
  byte_offset      = 1
  uC_offset        = byte_offset + (arankC - 1) * bytes_per_real

  byte_offset      = uC_size + 1
  uL_offset        = byte_offset + (arank - 1) * bytes_per_real

  byte_offset      = uC_size + uL_size + 1
  U_poisson_offset = byte_offset + (U_poisson_units * bytes_per_cplx) * rank

  byte_offset      = uC_size + uL_size + U_poisson_size + 1
  T_poisson_offset = byte_offset + (T_poisson_units * bytes_per_cplx) * rank

  byte_offset      = uC_size + uL_size + U_poisson_size + T_poisson_size + 1
  S_poisson_offset = byte_offset + (S_poisson_units * bytes_per_real) * rank

  byte_offset      = byte_offset + S_poisson_size
  byte_offset      = uC_size + uL_size + U_poisson_size + T_poisson_size + S_poisson_size + 1
  BJS_offset       = byte_offset + (BJS_units + BJS_pivot_units) * bytes_per_real * rank

  byte_offset      = uC_size + uL_size + U_poisson_size + T_poisson_size + S_poisson_size + &
                     BJS_size + BJS_pivot_size + 1
  S_coarse_offset  = byte_offset

  byte_offset      = uC_size + uL_size + U_poisson_size + T_poisson_size + S_poisson_size + &
                     BJS_size + BJS_pivot_size + S_coarse_size + S_coarse_pivot_size +  1
  uS_coarse_offset = byte_offset

  ! Read the kernel vectors to disk.
  read( 14850, pos=uC_offset ) uC
  read( 14850, pos=uL_offset ) uL

  ! Read U_poisson to disk in blocks.
  local_block_size = dimA**2 * bytes_per_cplx
  do ii = 1, numA_per_rank
     iioffset = U_poisson_offset + local_block_size * (ii - 1)
     read( 14850, pos=iioffset ) U_poisson(:, :, ii)
  enddo

  ! Read T_poisson to disk in blocks.
  local_block_size = dimA**2 * bytes_per_cplx
  do ii = 1, numA_per_rank
    iioffset = T_poisson_offset + local_block_size * (ii - 1)
    read( 14850, pos=iioffset ) T_poisson(:, :, ii)
  enddo

  ! Read S_poisson to disk in blocks.
  local_block_size = 8 * s**2 * bytes_per_real
  counter = 1
  do ii = 1, numA_per_rank
     do jj = 1, nky
        iioffset = S_poisson_offset + local_block_size * (counter - 1)
        read( 14850, pos=iioffset ) S_poisson(:, :, ii, jj)
        counter = counter + 1
     enddo
  enddo

  ! Read BJS (the block Jacobi Schur preconditioner) to disk in blocks.
  local_block_size = (16 * s**2 + 4 * s) * bytes_per_real
  counter = 1
  do ii = 1, numA_per_rank
     do jj = 1, nky
        iioffset = BJS_offset + local_block_size * (counter - 1)
        read( 14850, pos=iioffset ) BJS(:, :, ii, jj)
        iioffset = iioffset + 16 * s**2 * bytes_per_real
        read( 14850, pos=iioffset ) BJS_pivot(:, ii, jj)
        counter = counter + 1
     enddo
  enddo

  ! Read S_coarse and S_coarse_pivot.
  local_block_size = ((nsubx - 1)**2 + (nsubx - 1)) * bytes_per_real
  do jj = 1, nky
     iioffset = S_coarse_offset + local_block_size * (jj - 1)
     read( 14850, pos=iioffset ) S_coarse(:, :, jj)
     iioffset = iioffset + (nsubx - 1)**2 * bytes_per_real
     read( 14850, pos=iioffset ) S_coarse_pivot(:, jj)
  enddo

  ! Read the left kernel vector.
  read( 14850, pos=uS_coarse_offset ) uC_coarse

  ! Close the file.
  close( 14850 )

end subroutine read_setupfile
