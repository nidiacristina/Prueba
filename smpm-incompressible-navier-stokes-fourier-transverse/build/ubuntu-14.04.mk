# default parameters for a Ubuntu 14.04 system.

# gfortran is the only compiler available.
FC90            := gfortran

# use OpenBLAS as it's the fastest thing available when restricting oneself to
# the stock repositories.
BLAS            := OPENBLAS
ACML_DIR        :=
ATLAS_DIR       :=
MKL_DIR         :=
NETLIB_DIR      := /usr/lib
OPENBLAS_DIR    := /usr/lib

# FFTW3 and HDF5 have split installations beneath /usr.
FFTW_INC_DIR    := /usr/include
FFTW_LIB_DIR    := /usr/lib/x86_64-gnu-linux

HDF5_INC_DIR    := /usr/include
HDF5_LIB_DIR    := /usr/lib/x86_64-gnu-linux
HDF5_EXTRA_LIBS :=
