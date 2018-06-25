# default parameters for Stampede 1 (Sandybridge + MIC).

# ifort is the preferred compiler.
FC90            := ifort

# do some basic sanity checks and warn the user.
ifeq ($(FC90), gfortran)
# Intel compilers don't need to load MKL, though GNU compilers do.
ifeq ($(MKLROOT),)
$(warning "MKLROOT needs to be set!  Please load the 'mkl' module.")
endif
else # ~gfortran

# TACC only maintains FFTW3 for ifort.
ifeq ($(TACC_FFTW3_INC),)
$(warning "TACC_FFTW3_INC needs to be set! Please load the 'fftw3' module.")
endif
ifeq ($(TACC_FFTW3_LIB),)
$(warning "TACC_FFTW3_LIB needs to be set! Please load the 'fftw3' module.")
endif

endif # ~gfortran

# all compilers need HDF5 loaded.
ifeq ($(TACC_HDF5_INC),)
$(warning "TACC_HDF5_INC needs to be set! Please load the 'hdf5' module.")
endif
ifeq ($(TACC_HDF5_LIB),)
$(warning "TACC_HDF5_LIB needs to be set! Please load the 'hdf5' module.")
endif

# MKL is automatically used with ifort, though if gfortran is used we have
# limited options for BLAS/LAPACK as MKL is the only thing officially supported.
BLAS            := MKL
ACML_DIR        :=
ATLAS_DIR       :=
MKL_DIR         := $(MKLROOT)/lib/intel64
NETLIB_DIR      := $(HOME)/depot/stampede-sb/netlib
OPENBLAS_DIR    :=

# FFTW3 is only supported by TACC when the Intel compilers are used (unclear
# why), so we need to point to a local installation when we use something else.
ifeq ($(FC90), ifort)
FFTW_INC_DIR    := $(TACC_FFTW3_INC)
FFTW_LIB_DIR    := $(TACC_FFTW3_LIB)
else
FFTW_INC_DIR    := $(HOME)/depot/stampede-sb/include
FFTW_LIB_DIR    := $(HOME)/depot/stampede-sb/lib
endif

# HDF5 is maintained by TACC regardless of the compiler suite used.
#
# NOTE: we need to specify additional libraries as it is compiled with
#       non-standard support.
HDF5_INC_DIR    := $(TACC_HDF5_INC)
HDF5_LIB_DIR    := $(TACC_HDF5_LIB)
HDF5_EXTRA_LIBS := -lsz -lxml2
