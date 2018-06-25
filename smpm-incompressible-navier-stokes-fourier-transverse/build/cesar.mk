# default parameters for Cesar, an Ubuntu 14.04 workstation, at Cornell.

# gfortran is the preferred compiler.
FC90            := gfortran

# use ACML since that was the fastest thing evaluated the last time
# someone benchmarked the solver (circa mid-2015).
BLAS            := ACML
ACML_DIR        := /opt/acml5.3.1/gfortran64/lib
ATLAS_DIR       :=
MKL_DIR         :=
NETLIB_DIR      := /usr/lib
OPENBLAS_DIR    :=

# FFTW3 and HDF5 have split installations beneath /usr.
FFTW_INC_DIR    := /usr/include
FFTW_LIB_DIR    := /usr/lib/x86_64-gnu-linux

HDF5_INC_DIR    := /usr/include
HDF5_LIB_DIR    := /usr/lib/x86_64-gnu-linux
HDF5_EXTRA_LIBS :=
