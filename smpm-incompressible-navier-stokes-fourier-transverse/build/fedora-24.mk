# default parameters for a Fedora 24 system.

# gfortran is the only compiler available.
FC90            := gfortran

# use OpenBLAS as it's the fastest thing available when restricting oneself to
# the stock repositories.
BLAS            := OPENBLAS
ACML_DIR        :=
ATLAS_DIR       := /usr/lib64/atlas
MKL_DIR         := $(MKLROOT)
NETLIB_DIR      := /usr/lib
OPENBLAS_DIR    := /usr/lib/openblas-base

# FFTW3 and HDF5 are installed beneath /usr though Fedora's libraries can be
# either 32- or 64-bit.  we assume the latter since any non-trivial problem will
# require more than 4 GiB of RAM...
FFTW_INC_DIR    := /usr/include
FFTW_LIB_DIR    := /usr/lib64

HDF5_INC_DIR    := /usr/lib64/gfortran/modules
HDF5_LIB_DIR    := /usr/lib64
HDF5_EXTRA_LIBS :=
