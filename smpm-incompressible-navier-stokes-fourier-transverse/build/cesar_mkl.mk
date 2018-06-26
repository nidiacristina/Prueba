# default parameters for Cesar, an Ubuntu 14.04 workstation, at Cornell.

# gfortran is the preferred compiler.
FC90            := gfortran

# use ACML since that was the fastest thing evaluated the last time
# someone benchmarked the solver (circa mid-2015).
BLAS            := MKL
ACML_DIR        := 
OPENBLAS_DIR    := 
MKL_DIR         := /opt/intel/mkl/lib/intel64
NETLIB_DIR      := 

# FFTW3 and HDF5 have split installations beneath /usr.
FFTW_INC_DIR    := /opt/intel/mkl/include/fftw
FFTW_LIB_DIR    := /opt/intel/mkl/lib/intel64

HDF5_INC_DIR    := /opt/hdf5-1.8.17/hdf5/include
HDF5_LIB_DIR    := /opt/hdf5-1.8.17/hdf5/lib
HDF5_EXTRA_LIBS :=
