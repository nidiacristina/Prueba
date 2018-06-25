#pull in site-specific configuration before we set anything up, if it exists.
# those will override any defaultable variables below and saves the end user
# typing on the command line.
#
# NOTE: a number of system-specific configurations are found in build/ and
#       can either be copied to site.mk or symbolically linked.
-include site.mk

# parameters governing how we build the model.  by default, assume that we're
# not explicitly debugging and want things to run as fast as possible.
DEBUG          ?= yes
OPTIMIZE       ?= yes

# specify the BLAS/LAPACK we use.  by default, we use the reference Netlib
# library as it should be installed everywhere.
#
# NOTE: be warned, this is *SLOW* and should not be used for production.
#
BLAS           ?= NETLIB

# set the default Fortran 90 compiler, if nothing else has been defined.
#
# NOTE: this should be the compiler that is wrapped, not the wrapper itself.
#       that is, use "gfortran", "ifort", "pgf90", etc and not "mpif90",
#       "scorep", etc.
FC90           ?= gfortran

# set the default Fortran 90 compiler wrapper, if nothing else has been
# defined.
#
# NOTE: this should be the wrapper around the Fortran 90 compiler, and not
#       the Fortran 90 compiler itself.  that is, use "mpif90", "scorep", etc,
#       and not "gfortran", "ifort", "pgf90", etc.
FC90_WRAPPER   ?= mpif90

# configuration file for using Doxygen within the project.
DOXYGEN_CONFIG    := Doxyfile

# directory where Doxygen outputs are generated.
DOXYGEN_DIRECTORY := doxygen

# directory where ACML is installed.  a single-threaded implementation is
# preferred over the multi-threaded so we have absolute control over how
# the solver is parallelized.
ACML_DIR ?= /opt/acml5.3.1/gfortran64/lib

# directory where OpenBLAS' libraries are located.
#
# NOTE: this is the location of the shared and static libraries, not the
# 	    parent directory to the libraries.
OPENBLAS_DIR ?= /usr/lib/openblas-base

# directory where MKL libraries are stored.
#
# NOTE: this is not $MKLROOT, but rather one of the sub-directories beneath it
#       (e.g. $MKLROOT/lib/intel64).
#
MKL_DIR ?= /opt/intel/mkl

# directory where reference BLAS and LAPACK libraries are stored.
NETLIB_DIR ?= /usr/lib

# directory where FFTW's libraries and header files are located.  if only the
# root is set, then both the includes and library directories are set, otherwise
# they may be specified individually.
#
# NOTE: these are split as FFTW3 is commonly installed asymmetrically (Fedora,
#       Ubuntu, etc) so one cannot simply provide a root and build the necessary
#       search paths.
FFTW_ROOT    ?= /usr
FFTW_INC_DIR ?= $(FFTW_ROOT)/include
FFTW_LIB_DIR ?= $(FFTW_ROOT)/lib

# directory where HDF5's libraries and header files are located.
#
# NOTE: this is the location that they exist on cesar.
HDF5_ROOT    ?= /usr
HDF5_INC_DIR ?= $(HDF5_ROOT)/include
HDF5_LIB_DIR ?= $(HDF5_ROOT)/lib

# source code for the SMPM code.
WOODBURY_SOURCE = advance_time.f90 \
                  allocate_field_variables.f90 \
                  apply_3D_bc.f90 \
                  apply_3D_deflated_schur.f90 \
                  apply_3D_diffusive_fourier.f90 \
                  apply_3D_patching.f90 \
                  apply_3D_poisson.f90 \
                  apply_3D_preconditioned_schur.f90 \
                  apply_3D_schur_preconditioner.f90 \
                  apply_3D_poisson_schur_sparse.f90 \
                  apply_3D_vector_viscous_bc.f90 \
                  apply_3D_viscous_fourier.f90 \
                  apply_fft.f90 \
                  apply_filter_xz.f90 \
                  apply_filter_y.f90 \
                  apply_bc.f90 \
                  apply_capacitance_preconditioner.f90 \
                  apply_capacitance_preconditioner_transpose.f90 \
                  apply_extension.f90 \
                  apply_ifft.f90 \
                  apply_interfacial_averaging.f90 \
                  apply_nullspace_projection.f90 \
                  apply_patching.f90\
                  apply_pB_poisson.f90\
                  apply_pB_poisson_transpose.f90\
                  apply_pressure_bc.f90 \
                  apply_poisson_capacitance_sparse.f90 \
                  apply_poisson_capacitance_shift.f90 \
                  apply_preconditioned_capacitance.f90 \
                  apply_preconditioned_capacitance_shift.f90 \
                  apply_smpm_3D_viscous.f90 \
                  apply_smpm_3D_diffusion.f90 \
                  apply_smpm_advection.f90 \
                  apply_smpm_diffusion.f90 \
                  apply_smpm_poisson.f90 \
                  apply_smpm_poisson_transpose_parallel.f90 \
                  apply_smpm_viscous.f90 \
                  apply_smpm_transport.f90 \
                  apply_vector_viscous_bc.f90 \
                  apply_viscous_bc.f90 \
                  assemble_capacitance_block.f90 \
                  assemble_internal_patching.f90 \
                  assemble_laplacian.f90 \
                  check_cfl_number.f90 \
                  compute_3D_divergence.f90 \
                  compute_3D_gradient.f90 \
                  compute_3D_gradient_and_laplacian.f90 \
                  compute_3D_laplacian.f90 \
                  compute_3D_curl.f90 \
                  compute_curl.f90 \
                  compute_divergence.f90 \
                  compute_double_curl.f90 \
                  compute_gmres_householder.f90 \
                  compute_gmres_householder_complex.f90 \
                  compute_gradient.f90 \
                  compute_gradient_and_laplacian.f90 \
                  compute_laplacian.f90 \
                  compute_null_space_error.f90 \
                  compute_poisson_kernel.f90 \
                  compute_y_derivative.f90 \
                  deflation_matvecs.f90 \
                  derv.f90 \
                  display_inputs.f90 \
                  enforce_velocity_bc.f90 \
                  gll.f90 \
                  mod_constants.f90 \
                  mod_complex_blas.f90 \
                  mod_errors.f90 \
                  mod_field_variables.f90 \
                  mod_geom.f90 \
                  mod_io_diagnostics.f90 \
                  mod_io_field.f90 \
                  mod_io_restart.f90 \
                  mod_legendre.f90 \
                  mod_mesh_deformation_maps.f90 \
                  mod_nullspace.f90 \
                  mod_options.f90 \
                  mod_parallel_linear_algebra.f90 \
                  mod_postprocessor.f90 \
                  mod_precision.f90 \
                  mod_profiling.f90 \
                  mod_sponge.f90 \
                  mod_stopwatch.f90 \
                  mod_timestepping.f90 \
                  mod_transverse.f90 \
                  mod_woodbury_matrices.f90 \
                  main.f90 \
                  notify.f90 \
                  permutations.f90 \
                  quad.f90 \
                  read_input_file.f90 \
                  read_initfile_data.f90 \
                  read_restartfile_data.f90 \
                  read_setupfile.f90 \
                  seconds_since_epoch.f90 \
                  set_ownership_ranges.f90 \
                  set_timestepping_coefficients.f90 \
                  setup_2D_schur.f90 \
                  setup_3D_deflation.f90 \
                  setup_3D_poisson_solver.f90 \
                  setup_3D_schur_preconditioner.f90 \
                  setup_A_matrix.f90 \
                  setup_capacitance_preconditioner.f90 \
                  setup_cfl_calculation.f90 \
                  setup_deflation.f90 \
                  setup_deformation_maps.f90 \
                  setup_deformed_derivatives.f90 \
                  setup_filter_xz.f90 \
                  setup_filter_y.f90 \
                  setup_initial_conditions.f90 \
                  setup_nullspace_projection.f90 \
                  setup_poisson_operators.f90 \
                  setup_poisson_solver.f90 \
                  setup_sponge_layer.f90 \
                  setup_transverse.f90 \
                  solve_3D_A.f90 \
                  solve_3D_diffusion.f90 \
                  solve_3D_pressure.f90 \
                  solve_3D_poisson_iteratively.f90 \
                  solve_3D_viscous.f90 \
                  solve_A.f90 \
                  solve_tridiagonal_system.f90 \
                  synchronization.f90 \
                  write_fieldfile.f90 \
                  write_diagnosticsfile.f90 \
                  write_restartfile.f90 \
                  write_setupfile.f90 \
                  write_xdmf.f90

# source code for the SMPM post-processing code.
POSTPROC_SOURCE = allocate_field_variables.f90 \
                  apply_3D_deflated_schur.f90 \
                  apply_3D_poisson.f90 \
                  apply_3D_poisson_schur_sparse.f90 \
                  apply_3D_preconditioned_schur.f90 \
                  apply_3D_schur_preconditioner.f90 \
                  apply_bc.f90 \
                  apply_capacitance_preconditioner.f90 \
                  apply_capacitance_preconditioner_transpose.f90 \
                  apply_extension.f90 \
                  apply_filter_xz.f90 \
                  apply_filter_y.f90 \
                  apply_ifft.f90 \
                  apply_fft.f90 \
                  apply_patching.f90\
                  apply_pB_poisson.f90 \
                  apply_pB_poisson_transpose.f90\
                  apply_poisson_capacitance_sparse.f90 \
                  apply_poisson_capacitance_shift.f90 \
                  apply_preconditioned_capacitance.f90 \
                  apply_preconditioned_capacitance_shift.f90 \
                  apply_pressure_bc.f90 \
                  apply_smpm_advection.f90 \
                  apply_smpm_poisson.f90 \
                  apply_smpm_poisson_transpose_parallel.f90 \
                  apply_streamfunction_matrix.f90 \
                  assemble_capacitance_block.f90 \
                  assemble_internal_patching.f90 \
                  assemble_laplacian.f90 \
                  compute_3D_divergence.f90 \
                  compute_3D_gradient.f90 \
                  compute_3D_laplacian.f90 \
                  compute_3D_curl.f90 \
                  compute_curl.f90 \
                  compute_divergence.f90 \
                  compute_double_curl.f90 \
                  compute_gmres_householder.f90 \
                  compute_gradient.f90 \
                  compute_gradient_and_laplacian.f90 \
                  compute_laplacian.f90 \
                  compute_poisson_kernel.f90 \
                  compute_y_derivative.f90 \
                  deflation_matvecs.f90 \
                  derv.f90 \
                  display_inputs.f90 \
                  gll.f90 \
                  main_post_processor.f90 \
                  make_streamfunction.f90 \
                  mod_constants.f90 \
                  mod_errors.f90 \
                  mod_field_variables.f90 \
                  mod_geom.f90 \
                  mod_io_post.f90 \
                  mod_legendre.f90 \
                  mod_mesh_deformation_maps.f90 \
                  mod_options.f90 \
                  mod_parallel_linear_algebra.f90 \
                  mod_postprocessor.f90 \
                  mod_precision.f90 \
                  mod_sponge.f90 \
                  mod_woodbury_matrices.f90 \
                  mod_timestepping.f90 \
                  mod_transverse.f90 \
                  notify.f90 \
                  permutations.f90 \
                  quad.f90 \
                  read_input_file.f90 \
                  read_initfile_data.f90 \
                  read_restartfile_data.f90 \
                  read_setupfile.f90 \
                  set_ownership_ranges.f90 \
                  set_timestepping_coefficients.f90 \
                  setup_2D_schur.f90 \
                  setup_3D_deflation.f90 \
                  setup_3D_poisson_solver.f90 \
                  setup_3D_schur_preconditioner.f90 \
                  setup_A_matrix.f90 \
                  setup_filter_xz.f90 \
                  setup_filter_y.f90 \
                  setup_capacitance_preconditioner.f90 \
                  setup_cfl_calculation.f90 \
                  setup_deflation.f90 \
                  setup_deformation_maps.f90 \
                  setup_deformed_derivatives.f90 \
                  setup_poisson_operators.f90 \
                  setup_transverse.f90 \
                  solve_3D_poisson_iteratively.f90 \
                  solve_3D_pressure.f90 \
                  solve_3D_A.f90 \
                  solve_A.f90 \
                  solve_tridiagonal_system.f90 \
                  synchronization.f90 \
                  write_postfile.f90 \
                  write_postprocessor_xdmf.f90 \
                  write_setupfile.f90

VALIDATOR_SOURCE = allocate_field_variables.f90 \
                   apply_3D_deflated_schur.f90 \
                   apply_3D_diffusive_fourier.f90 \
                   apply_3D_poisson.f90 \
                   apply_3D_poisson_schur_sparse.f90 \
                   apply_3D_preconditioned_schur.f90 \
                   apply_3D_schur_preconditioner.f90 \
                   apply_3D_viscous_fourier.f90 \
                   apply_bc.f90 \
                   apply_capacitance_preconditioner_transpose.f90 \
                   apply_extension.f90 \
                   apply_fft.f90 \
                   apply_filter_y.f90 \
                   apply_ifft.f90 \
                   apply_patching.f90 \
                   apply_pB_poisson.f90 \
                   apply_pB_poisson_transpose.f90 \
                   apply_poisson_capacitance_shift.f90 \
                   apply_poisson_capacitance_sparse.f90 \
                   apply_preconditioned_capacitance_shift.f90 \
                   apply_smpm_advection.f90 \
                   apply_smpm_diffusion.f90 \
                   apply_smpm_poisson.f90 \
                   apply_smpm_poisson_transpose_parallel.f90 \
                   apply_smpm_transport.f90 \
                   apply_vector_viscous_bc.f90 \
                   apply_viscous_bc.f90 \
                   assemble_capacitance_block.f90 \
                   assemble_internal_patching.f90 \
                   assemble_laplacian.f90 \
                   compute_3D_divergence.f90 \
                   compute_3D_curl.f90 \
                   compute_3D_gradient.f90 \
                   compute_divergence.f90 \
                   compute_double_curl.f90 \
                   compute_gmres_householder.f90 \
                   compute_gmres_householder_complex.f90 \
                   compute_gradient.f90 \
                   compute_gradient_and_laplacian.f90 \
                   compute_3D_laplacian.f90 \
                   compute_laplacian.f90 \
                   compute_poisson_kernel.f90 \
                   compute_y_derivative.f90 \
                   deflation_matvecs.f90 \
                   derv.f90 \
                   display_inputs.f90 \
                   gll.f90 \
                   main_validator.f90 \
                   mod_command_line.f90 \
                   mod_complex_blas.f90 \
                   mod_constants.f90 \
                   mod_errors.f90 \
                   mod_field_variables.f90 \
                   mod_geom.f90 \
                   mod_io_diagnostics.f90 \
                   mod_io_field.f90 \
                   mod_io_restart.f90 \
                   mod_legendre.f90 \
                   mod_mesh_deformation_maps.f90 \
                   mod_nullspace.f90 \
                   mod_options.f90 \
                   mod_parallel_linear_algebra.f90 \
                   mod_postprocessor.f90 \
                   mod_sponge.f90 \
                   mod_timestepping.f90 \
                   mod_transverse.f90 \
                   mod_validation.f90 \
                   mod_woodbury_matrices.f90 \
                   notify.f90 \
                   permutations.f90 \
                   quad.f90 \
                   read_initfile_data.f90 \
                   read_input_file.f90 \
                   read_restartfile_data.f90 \
                   read_setupfile.f90 \
                   setup_2D_schur.f90 \
                   setup_3D_deflation.f90 \
                   setup_3D_poisson_solver.f90 \
                   setup_3D_schur_preconditioner.f90 \
                   setup_A_matrix.f90 \
                   setup_capacitance_preconditioner.f90 \
                   setup_cfl_calculation.f90 \
                   setup_deflation.f90 \
                   setup_deformation_maps.f90 \
                   setup_deformed_derivatives.f90 \
                   setup_filter_xz.f90 \
                   setup_filter_y.f90 \
                   setup_initial_conditions.f90 \
                   setup_nullspace_projection.f90 \
                   set_ownership_ranges.f90 \
                   setup_poisson_operators.f90 \
                   setup_sponge_layer.f90 \
                   setup_transverse.f90 \
                   solve_A.f90 \
                   solve_3D_A.f90 \
                   solve_3D_diffusion.f90 \
                   solve_3D_poisson_iteratively.f90 \
                   solve_3D_viscous.f90 \
                   solve_tridiagonal_system.f90 \
                   synchronization.f90 \
                   validate_advection_terms.f90 \
                   validate_fftw.f90 \
                   validate_divergence.f90 \
                   validate_double_curl.f90 \
                   validate_filter_y.f90 \
                   validate_gradients.f90 \
                   validate_helmholtz_2d_scalar_dirichlet.f90 \
                   validate_helmholtz_3d_scalar_dirichlet.f90 \
                   validate_helmholtz_3d_scalar_neumann.f90 \
                   validate_helmholtz_3d_vector_dirichlet.f90 \
                   validate_helmholtz_3d_vector_neumann.f90 \
                   validate_laplacian.f90 \
                   validate_poisson_3d_iteratively.f90 \
                   validate_real_valued_gmres_3D.f90 \
                   validate_transport_term.f90 \
                   write_fieldfile.f90 \
                   write_diagnosticsfile.f90 \
                   write_setupfile.f90

# Intel Fortran interface files.  ifort builds an explicit interface file for
# each sub-routine and function it encounters so that they can be checked.
# attempt to keep track of them.
#
# NOTE: we cannot (easily) build a list of the "*__genmod.*" files since there
#       is one per sub-routine/function.
#
# NOTE: create this list here, rather than in the ifort-specific configuration
#       so we don't need to specify FC90=ifort to remove these when building.
#       the 'clean' target.
IFORT_GENMOD_FILES = *__genmod.f90 *__genmod.mod

# object file lists.
WOODBURY_OBJECTS  = $(patsubst %.f90, %.o, $(WOODBURY_SOURCE))
POSTPROC_OBJECTS  = $(patsubst %.f90, %.o, $(POSTPROC_SOURCE))
VALIDATOR_OBJECTS = $(patsubst %.f90, %.o, $(VALIDATOR_SOURCE))

# build a list of Fortran module files from the source file list.  take care
# to translate all of the file names into lower case since the module files
# generated with only lower case.
WOODBURY_MODULES  := $(patsubst mod_%.f90, %.mod, $(filter mod_%.f90, \
                                                  $(shell echo $(WOODBURY_SOURCE) | tr 'A-Z' 'a-z')))
POSTPROC_MODULES  := $(patsubst mod_%.f90, %.mod, $(filter mod_%.f90, \
                                                  $(shell echo $(POSTPROC_SOURCE) | tr 'A-Z' 'a-z')))
VALIDATOR_MODULES := $(patsubst mod_%.f90, %.mod, $(filter mod_%.f90, \
                                                  $(shell echo $(VALIDATOR_SOURCE) | tr 'A-Z' 'a-z')))

# the woodbury model is built from everything in the project.
WOODBURY_CMD     = smpm_incompressible_navier_stokes
POSTPROC_CMD     = smpm_post_processor
VALIDATOR_CMD    = smpm_validator

# the solver's field file comparator.
SMPM_FIELD_COMPARE_CMD = smpm_compare_fieldfiles.py

# tags files contain indices of variables and functions found within the SMPM
# source code.  both vim and Emacs can use these to open the files where a
# symbol is defined and/or used.
VIM_TAGS                 = tags
EMACS_TAGS               = TAGS

# global, generic variables that contain each of the major targets in this
# Makefile.
SOURCE                   = $(WOODBURY_SOURCE) $(POSTPROC_SOURCE) $(VALIDATOR_SOURCE)
OBJECTS                  = $(WOODBURY_OBJECTS) $(POSTPROC_OBJECTS) $(VALIDATOR_OBJECTS)
EXECUTABLES              = $(WOODBURY_CMD) $(POSTPROC_CMD) $(VALIDATOR_CMD)
MODULES                  = $(WOODBURY_MODULES) $(POSTPROC_MODULES) $(VALIDATOR_MODULES)
TAGS_FILES               = $(VIM_TAGS) $(EMACS_TAGS)

# list of things that are created during development that can be removed, and
# easily regenerated, if necessary.
CLEAN_LIST               = $(OBJECTS) $(EXECUTABLES) $(MODULES) $(IFORT_GENMOD_FILES) \
                           $(TAGS_FILES)

# list of things that are created during development that can be removed,
# but aren't trivial to regenerate.
DISTCLEAN_LIST           = $(DOXYGEN_DIRECTORY)

# driver for the test cases, validation cases, and the benchmarks.
WOODBURY_TESTER      = run-woodbury-tests.sh
WOODBURY_VALIDATOR   = run-woodbury-validators.sh
WOODBURY_BENCHMARKER = run-woodbury-benchmarks.sh

# directory where test and validation cases are stored.
TESTING_DIRECTORY    = testing/
VALIDATION_DIRECTORY = validation/

# make sure that we're able to compile and link against FFTW.
F90FLAGS := -I$(FFTW_INC_DIR)
LIBS     := -L$(FFTW_LIB_DIR) -Wl,-rpath,$(FFTW_LIB_DIR) -lfftw3

# determine the F90 compilation and linking flags.
#
# NOTE: we use the $(filter) function throughout the compiler, and wrapper,
#       selections so that custom versions of individual tools do not require
#       special support.  in particular, things like 'gfortran-4.8.2' will
#       match 'gfortran'.
ifneq (,$(filter gfortran, $(FC90)))

# gfortran does not have BLAS/LAPACK implementations, so select which versions
# to link against.
#
# NOTE: ensure that we can find the libraries we build with at run-time by
#       specifying -rpath.
ifeq ($(BLAS), ACML)
ifeq ($(ACML_DIR),)
$(error "ACML was requested but ACML_DIR is empty.  Please set it!")
endif
LIBS += -L$(ACML_DIR) -lacml -Wl,-rpath,$(ACML_DIR)
else ifeq ($(BLAS), ACML_OMP)
ifeq ($(ACML_DIR),)
$(error "ACML was requested but ACML_DIR is empty.  Please set it!")
endif
LIBS += -L$(ACML_DIR) -lacml_mp -Wl,-rpath,$(ACML_DIR)
else ifeq ($(BLAS), ATLAS)
ifeq ($(ATLAS_DIR),)
$(error "ATLAS was requested by ATLAS_DIR is empty.  Please set it!")
endif
LIBS += -L$(ATLAS_DIR) -Wl,-rpath,$(ATLAS_DIR) -lsatlas
else ifeq ($(BLAS), OPENBLAS)
ifeq ($(OPENBLAS_DIR),)
$(error "OpenBLAS was requested but OPENBLAS_DIR is empty.  Please set it!")
endif
ifeq ($(NETLIB_DIR),)
$(error "OpenBLAS was requested but NETLIB_DIR (for LAPACK) is empty.  Please set it!")
endif
LIBS += -L$(OPENBLAS_DIR) -L$(NETLIB_DIR) -Wl,-rpath,$(OPENBLAS_DIR) -Wl,-rpath,$(NETLIB_DIR) -lopenblas -llapack
else ifeq ($(BLAS), MKL)
ifeq ($(MKL_DIR),)
$(error "MKL was requested but MKL_DIR is empty.  Please set it!")
endif
LIBS += -L$(MKL_DIR) -Wl,-rpath,$(MKL_DIR) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
else ifeq ($(BLAS), MKL_ILP64)
ifeq ($(MKL_DIR),)
$(error "MKL was requested but MKL_DIR is empty.  Please set it!")
endif
LIBS += -L$(MKL_DIR) -Wl,-rpath,$(MKL_DIR) -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
else ifeq ($(BLAS), NETLIB)
ifeq ($(NETLIB_DIR),)
$(error "Netlib was requested but NETLIB_DIR is empty.  Please set it!")
endif
# default to the system's implementation which is likely to be slow.
LIBS += -L$(NETLIB_DIR) -Wl,-rpath,$(NETLIB_DIR) -llapack -lblas
else
$(error "Unknown BLAS requested ($(BLAS))!")
endif

# compile without a maximum line length.  turn on all warnings so that
# questionable code is noticed early on.
F90FLAGS  += -c -ffixed-line-length-none -Wall -Wno-maybe-uninitialized

ifeq ($(PEDANTIC),yes)
# if we're grooming the code for nits, turn on a slew of warnings that aren't
# otherwise enabled because they generate a lot of noise with bad code.
F90FLAGS += -Wconversion -pedantic -Warray-temporaries -Wcharacter-truncation -Wsurprising -Wintrinsic-shadow -Wextra -pedantic
endif

# generate the most verbose debugging information we can in compiled objects,
# regardless of the request for debugging.  this allows for improved profiling
# at the cost of some extra diskspace.
F90FLAGS += -ggdb -g3

ifneq ($(DEBUG),no)
# enable array bounds checking and ensure that invalid floating point
# operations signal an exception, and output backtraces on halt.
F90FLAGS += -fbounds-check -fcheck=all -ffpe-trap=invalid -fbacktrace
endif

ifeq ($(OPTIMIZE),yes)
# enable fast math optimizations (which may not be strictly IEEE-754 compliant),
# unroll loops where it makes sense.
F90FLAGS += -O3 -funroll-loops -march=native -mtune=native
endif

ifneq ($(OPENMP),no)
# if we haven't been requested to disable OpenMP, enable it during compilation
# and linking.
F90FLAGS += -fopenmp
LIBS     += -fopenmp
endif

ifneq ($(PLUGIN),)
# the GCC compiler suite supports plugins used for compilation.  pass it onto
# the compilation flags if the user has supplied one.
F90FLAGS += $(PLUGIN)
endif

else # FC90 !~ gfortran

ifneq (,$(filter ifort, $(FC90)))
# BLAS and LAPACK are provided by Intel's Math Kernel Library (MKL) when using
# Intel's Fortran compiler.
LIBS += -mkl=sequential

# compile assuming ANSI aliasing rules, and warn about everything the compiler
# is concerned about.
F90FLAGS += -c -ansi-alias -warn all -mkl=sequential

ifeq ($(PEDANTIC),yes)
# if we're grooming the code for nits, turn on additional warnings that aren't
# otherwise enabled because they generate a lot of noise with bad code.
F90FLAGS += -debug all
endif

# generate the most verbose debugging information we can in compiled objects,
# regardless of the request for debugging.  this allows for improved profiling
# at the cost of some extra diskspace.
F90FLAGS += -g

ifneq ($(DEBUG),no)
# enable array bounds checking and ensure that invalid floating point
# operations signal an exception.
# XXX: fpe0 causes underflow to yield 0.0, abort on other exception
F90FLAGS += -fpe0 -debug full -check bounds -check all -check uninit -traceback
endif

ifeq ($(OPTIMIZE),yes)
# optimize for the host architecture and align arrays to 32-byte boundaries for
# maximum performance wrt vectorization.  enable inter-procedural optimizations
# across all translation units, unroll loops with a compiler-detected trip
# count, and use OpenMP for threading.
F90FLAGS += -xHost -align array32byte -ipo -unroll
else
# NOTE: we explicitly disable optimizations as Intel's compilers default to -O2
#       when nothing else has been specified.  this guards against unintentional
#       optimization and makes OPTIMIZE=no more intuitive.
F90FLAGS += -O0
endif

ifneq ($(OPENMP),no)
# if we haven't been requested to disable OpenMP, enable it during compilation
# and linking.
F90FLAGS += -qopenmp
LIBS     += -qopenmp
endif

else # FC90 !~ ifort
$(warning Unknown F90 compiler - $(FC90)!)
endif
endif

# determine how to compile and link against the HDF5 libraries.
F90FLAGS += -I$(HDF5_INC_DIR)
LIBS     += -L$(HDF5_LIB_DIR) -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz $(HDF5_EXTRA_LIBS) -ldl

# determine which Fortran 90 wrapper we're using.
ifneq (,$(filter mpif90, $(FC90_WRAPPER)))
# nothing is needed for MPI.  all compilation flags are passed directly to
# $(FC90).
else # FC90_WRAPPER !~ mpif90
ifneq (,$(filter scorep, $(FC90_WRAPPER)))
override FC90_WRAPPER += --mpi --openmp mpif90

# Score-P assumes that OpenMP directives always generate OpenMP-based code.
# since we only enable OpenMP when requested, we need to turn it on in the
# non-requested case so that we do not get cryptic link errors about
# mismatches in thread-level storage (TLS)-related symbols.
ifneq ($(OPENMP),yes)

# NOTE: unfortunately, we have to handle each of the different types of
#       compilers so that we ensure OpenMP is properly enabled.
ifneq (,$(filter gfortran, $(FC90)))
F90FLAGS += -fopenmp
LIBS     += -fopenmp
else
ifneq (,$(filter ifort, $(FC90)))
F90FLAGS += -openmp
LIBS     += -openmp
endif
endif
endif

else
$(warning Unknown F90 wrapper - $(FC90_WRAPPER))
endif
endif

# if we're explicitly debugging, run the test harness in debugging mode.
ifeq ($(DEBUG),yes)
TEST_HARNESS_FLAGS += -d
endif

# if the user has specified a custom MPI launcher, tell the harnesses about
# it.
ifneq ($(MPI_LAUNCHER),)
TEST_HARNESS_FLAGS      += -L $(MPI_LAUNCHER)
VALIDATOR_HARNESS_FLAGS += -L $(MPI_LAUNCHER)
BENCHMARK_HARNESS_FLAGS += -L $(MPI_LAUNCHER)
endif

# if the user has requested a specific set of MPI ranks to run benchmarks and
# tests with, pass it to the harnesses.
ifneq ($(MPI_RANKS),)
TEST_HARNESS_FLAGS      += -p $(MPI_RANKS)
VALIDATOR_HARNESS_FLAGS += -p $(MPI_RANKS)
BENCHMARK_HARNESS_FLAGS += -M $(MPI_RANKS)
endif

# if the user has requested a specific set of OpenMP threads to run benchmarks
# with, pass it to the benchmark harness.
#
# NOTE: the test framework does not support directly specifying a number of
#       of OpenMP threads to use.  the OMP_NUM_THREADS environment variable
#       directly influences its execution state.
ifneq ($(OMP_THREADS),)
BENCHMARK_HARNESS_FLAGS += -O $(OMP_THREADS)
endif

# request the harnesses to execute a subset of each test.  this allows offline
# preparation of test cases to reduce run-time dependencies.
ifneq ($(TEST_ACTIONS),)
TEST_HARNESS_FLAGS      += -A $(TEST_ACTIONS)
VALIDATOR_HARNESS_FLAGS += -A $(TEST_ACTIONS)
endif

####

# ensure our unconditional targets cannot be tripped up by the existence of
# files on disk.
.PHONY: print_environment

# by default, build the model and its validator.
#
# NOTE: ideally the print_environment dependency should be an explicit
#       dependency of each of the executables we built so it is printed
#       regardless of whether the default target is built or not.
#       unfortunately, that results in a forced rebuild of each of the
#       executables the dependency is added to.
#
all: print_environment $(WOODBURY_CMD) $(POSTPROC_CMD) $(VALIDATOR_CMD) tags-files

# link the executables.
$(WOODBURY_CMD): $(WOODBURY_OBJECTS)
	$(FC90_WRAPPER) -o $@ $(WOODBURY_OBJECTS) $(LIBS) 

$(POSTPROC_CMD): $(POSTPROC_OBJECTS)
	$(FC90_WRAPPER) -o $@ $(POSTPROC_OBJECTS) $(LIBS) 

$(VALIDATOR_CMD): $(VALIDATOR_OBJECTS)
	$(FC90_WRAPPER) -o $@ $(VALIDATOR_OBJECTS) $(LIBS) 

%.o: %.f90
	$(FC90_WRAPPER) -o $@ $(F90FLAGS) $<

%.o: %.F90
	$(FC90_WRAPPER) -o $@ $(F90FLAGS) $<

clean clena:
	-rm -rf $(CLEAN_LIST)

distclean distclena:
	-rm -rf $(DISTCLEAN_LIST)

tags-files: $(TAGS_FILES)

$(VIM_TAGS): $(SOURCE)
	ctags -o $@ $^

$(EMACS_TAGS): $(SOURCE)
	etags -o $@ $^

test: $(WOODBURY_CMD)
	@$(VERBOSE) ./$(WOODBURY_TESTER) $(TEST_HARNESS_FLAGS) $(shell pwd)/$(WOODBURY_CMD) $(shell pwd)/$(SMPM_FIELD_COMPARE_CMD) $(TESTING_DIRECTORY) $(TEST_CASES)

validate: $(WOODBURY_CMD) $(VALIDATOR_CMD)
	@$(VERBOSE) ./$(WOODBURY_VALIDATOR) $(VALIDATOR_HARNESS_FLAGS) $(shell pwd)/$(VALIDATOR_CMD) $(VALIDATION_DIRECTORY) $(VALIDATION_CASES)

benchmark: $(WOODBURY_CMD)
	@echo "Not implemented!"
#	@$(VERBOSE) ./$(WOODBURY_BENCHMARKER) $(BENCHMARK_HARNESS_FLAGS) $(shell pwd)/$(WOODBURY_CMD) $(shell pwd)/$(SMPM_LOG_VIEWER_CMD) $(TESTING_DIRECTORY) $(TEST_CASES)

$(DOXYGEN_DIRECTORY):
	@doxygen $(DOXYGEN_CONFIG)

print_environment:
	@echo "Current build parameters:"
	@echo
	@echo "    Fortran compiler:  $(FC90)"
	@echo "         MPI wrapper:  $(FC90_WRAPPER)"
	@echo
	@echo "           Optimized:  $(OPTIMIZE)"
	@echo "       Debug symbols:  $(DEBUG)"
ifneq ($(OPENMP), no)
	@echo "              OpenMP:  Enabled"
else
	@echo "              OpenMP:  Disabled"
endif
	@echo
ifeq ($(BLAS), ACML)
	@echo "                ACML:  $(ACML_DIR)"
else ifeq ($(BLAS), ACML_MP)
	@echo " Multi-threaded ACML:  $(ACML_DIR)"
else ifeq ($(BLAS), MKL)
	@echo "      Sequential MKL:  $(MKL_DIR)"
else ifeq ($(BLAS), NETLIB)
	@echo "    Reference Netlib:  $(NETLIB_DIR)"
	$(warning "You're linking against the *SLOW* reference BLAS/LAPACK implementation!")
	$(warning "Consider changing the BLAS variable!")
else ifeq ($(BLAS), OPENBLAS)
	@echo "            OpenBLAS:  $(OPENBLAS_DIR)"
endif
	@echo
	@echo "        HDF5 headers:  $(HDF5_INC_DIR)"
	@echo "      HDF5 libraries:  $(HDF5_LIB_DIR)"
ifneq ($(HDF5_EXTRA_LIBS),)
	@echo "Extra HDF5 libraries:  $(HDF5_EXTRA_LIBS)"
else
	@echo "Extra HDF5 libraries:  None"
endif
	@echo
	@echo "       FFTW3 headers:  $(FFTW_INC_DIR)"
	@echo "     FFTW3 libraries:  $(FFTW_LIB_DIR)"
	@echo
	@echo
	@echo "    NOTE: This configuration only applies to files compiled and executables"
	@echo "          linked with this invocation of make and may not reflect the"
	@echo "          configuration of existing files."
	@echo

# explicit dependencies between objects and the modules they require for
# compilation.
advance_time.o: mod_field_variables.o mod_timestepping.o
allocate_field_variables.o: mod_constants.o mod_field_variables.o mod_geom.o mod_legendre.o mod_mesh_deformation_maps.o mod_precision.o mod_sponge.o mod_transverse.o mod_woodbury_matrices.o
apply_3D_bc.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
apply_3D_deflated_schur.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
apply_3D_diffusive_fourier.o: mod_constants.o mod_options.o mod_parallel_linear_algebra.o mod_precision.o mod_timestepping.o mod_transverse.o mod_woodbury_matrices.o
apply_3D_patching.o: mod_constants.o mod_mesh_deformation_maps.o mod_precision.o mod_woodbury_matrices.o
apply_3D_poisson.o: mod_constants.o mod_options.o mod_precision.o mod_transverse.o mod_woodbury_matrices.o
apply_3D_poisson_schur_sparse.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
apply_3D_preconditioned_schur.f90: mod_constants.o mod_precision.o mod_woodbury_matrices.o
apply_3D_schur_preconditioner.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
apply_3D_vector_viscous_bc.o: mod_constants.o mod_mesh_deformation_maps.o mod_precision.o mod_woodbury_matrices.o
apply_3D_viscous_fourier.o: mod_constants.o mod_options.o mod_precision.o mod_timestepping.o mod_transverse.o mod_woodbury_matrices.o
apply_bc.o: mod_constants.o mod_mesh_deformation_maps.o mod_precision.o mod_woodbury_matrices.o
apply_capacitance_preconditioner.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
apply_capacitance_preconditioner_transpose.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
apply_extension.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
apply_fft.o: mod_constants.o mod_precision.o mod_transverse.o mod_woodbury_matrices.o
apply_filter_xz.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
apply_filter_y.o: mod_constants.o mod_precision.o mod_transverse.o mod_woodbury_matrices.o
apply_ifft.o: mod_constants.o mod_precision.o mod_transverse.o mod_woodbury_matrices.o
apply_interfacial_averaging.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
apply_nullspace_projection.o: mod_constants.o mod_nullspace.o mod_precision.o mod_woodbury_matrices.o
apply_pB_poisson.o: mod_constants.o mod_mesh_deformation_maps.o mod_options.o mod_precision.o mod_woodbury_matrices.o
apply_pB_poisson_transpose.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
apply_patching.o: mod_constants.o mod_mesh_deformation_maps.o mod_precision.o mod_woodbury_matrices.o
apply_poisson_capacitance_shift.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
apply_poisson_capacitance_sparse.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
apply_preconditioned_capacitance.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
apply_preconditioned_capacitance_shift.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
apply_pressure_bc.o: mod_constants.o mod_field_variables.o mod_mesh_deformation_maps.o mod_precision.o mod_timestepping.o mod_woodbury_matrices.o
apply_smpm_3D_diffusion.o: mod_constants.o mod_options.o mod_precision.o mod_timestepping.o mod_woodbury_matrices.o
apply_smpm_3D_viscous.o: mod_constants.o mod_options.o mod_precision.o mod_timestepping.o mod_woodbury_matrices.o
apply_smpm_advection.o: mod_constants.o mod_field_variables.o mod_mesh_deformation_maps.o mod_precision.o mod_woodbury_matrices.o
apply_smpm_diffusion.o: mod_constants.o mod_options.o mod_precision.o mod_timestepping.o mod_woodbury_matrices.o
apply_smpm_poisson.o: mod_options.o mod_precision.o mod_woodbury_matrices.o
apply_smpm_poisson_transpose_parallel.o: mod_constants.o mod_parallel_linear_algebra.o mod_precision.o mod_woodbury_matrices.o
apply_smpm_transport.o: mod_constants.o mod_field_variables.o mod_mesh_deformation_maps.o mod_precision.o mod_woodbury_matrices.o
apply_smpm_viscous.o: mod_constants.o mod_options.o mod_precision.o mod_timestepping.o mod_woodbury_matrices.o
apply_streamfunction_matrix.o: mod_precision.o mod_woodbury_matrices.o
apply_vector_viscous_bc.o: mod_constants.o mod_mesh_deformation_maps.o mod_precision.o mod_woodbury_matrices.o
apply_viscous_bc.o: mod_constants.o mod_mesh_deformation_maps.o mod_precision.o mod_woodbury_matrices.o
assemble_capacitance_block.o: mod_constants.o mod_options.o mod_precision.o mod_woodbury_matrices.o
assemble_internal_patching.o: mod_constants.o mod_mesh_deformation_maps.o mod_precision.o mod_woodbury_matrices.o
assemble_laplacian.o: mod_constants.o mod_legendre.o mod_mesh_deformation_maps.o mod_options.o mod_precision.o mod_woodbury_matrices.o
check_cfl_number.o: mod_constants.o mod_field_variables.o mod_mesh_deformation_maps.o mod_options.o mod_precision.o mod_timestepping.o mod_transverse.o mod_woodbury_matrices.o
compute_3D_curl.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
compute_3D_divergence.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
compute_3D_gradient.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
compute_3D_gradient_and_laplacian.o: mod_constants.o mod_precision.o mod_transverse.o mod_woodbury_matrices.o
compute_3D_laplacian.o: mod_constants.o mod_precision.o mod_transverse.o mod_woodbury_matrices.o
compute_curl.o: mod_precision.o mod_woodbury_matrices.o
compute_divergence.o: mod_constants.o mod_legendre.o mod_mesh_deformation_maps.o mod_precision.o mod_woodbury_matrices.o
compute_double_curl.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
compute_gmres_householder.o: mod_constants.o mod_parallel_linear_algebra.o mod_precision.o
compute_gmres_householder_complex.o: mod_constants.o mod_parallel_linear_algebra.o mod_precision.o
compute_gradient.o: mod_constants.o mod_legendre.o mod_mesh_deformation_maps.o mod_precision.o mod_woodbury_matrices.o
compute_gradient_and_laplacian.o: mod_constants.o mod_legendre.o mod_mesh_deformation_maps.o mod_precision.o mod_woodbury_matrices.o
compute_laplacian.o: mod_constants.o mod_legendre.o mod_mesh_deformation_maps.o mod_precision.o mod_woodbury_matrices.o
compute_null_space_error.o: mod_constants.o mod_errors.o mod_options.o mod_parallel_linear_algebra.o mod_precision.o mod_woodbury_matrices.o
compute_poisson_kernel.o: mod_constants.o mod_options.o mod_parallel_linear_algebra.o mod_precision.o mod_woodbury_matrices.o
compute_y_derivative.o: mod_constants.o mod_precision.o mod_transverse.o mod_woodbury_matrices.o
deflation_matvecs.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
derv.o: mod_precision.o
display_inputs.o: mod_constants.o mod_options.o mod_woodbury_matrices.o
enforce_velocity_bc.o: mod_constants.o mod_mesh_deformation_maps.o mod_options.o mod_precision.o mod_woodbury_matrices.o
gll.o: mod_precision.o
main.o: mod_constants.o mod_errors.o mod_field_variables.o mod_legendre.o mod_options.o mod_parallel_linear_algebra.o mod_precision.o mod_profiling.o mod_sponge.o mod_stopwatch.o mod_timestepping.o mod_transverse.o mod_woodbury_matrices.o
main_post_processor.o: mod_constants.o mod_field_variables.o mod_legendre.o mod_options.o mod_parallel_linear_algebra.o mod_precision.o mod_woodbury_matrices.o
main_validator.o: mod_command_line.o mod_constants.o mod_errors.o mod_field_variables.o mod_options.o mod_precision.o mod_validation.o mod_timestepping.o mod_woodbury_matrices.o
make_streamfunction.o: mod_constants.o mod_legendre.o mod_options.o mod_parallel_linear_algebra.o mod_precision.o mod_woodbury_matrices.o
mod_command_line.o: mod_precision.o
mod_complex_blas.o: mod_precision.o
mod_constants.o: mod_precision.o
mod_errors.o: mod_precision.o
mod_field_variables.o: mod_precision.o
mod_geom.o: mod_precision.o
mod_io_diagnostics.o:
mod_io_field.o:
mod_io_post.o:
mod_io_restart.o:
mod_legendre.o: mod_precision.o
mod_mesh_deformation_maps.o: mod_precision.o
mod_nullspace.o: mod_precision.o
mod_options.o: mod_precision.o
mod_parallel_linear_algebra.o: mod_precision.o
mod_postprocessor.o: mod_precision.o
mod_precision.o:
mod_profiling.o: mod_precision.o
mod_sponge.o: mod_precision.o
mod_stopwatch.o: mod_precision.o
mod_timestepping.o: mod_precision.o
mod_transverse.o: mod_precision.o
mod_validation.o: mod_precision.o
mod_woodbury_matrices.o: mod_precision.o
notify.o: mod_constants.o
permutations.o: mod_precision.o
quad.o: mod_precision.o
read_initfile_data.o: mod_constants.o mod_field_variables.o mod_geom.o mod_io_field.o mod_precision.o mod_transverse.o mod_woodbury_matrices.o
read_input_file.o: mod_constants.o mod_options.o mod_sponge.o mod_timestepping.o mod_woodbury_matrices.o
read_restartfile_data.o: mod_constants.o mod_field_variables.o mod_geom.o mod_io_restart.o mod_precision.o mod_timestepping.o mod_transverse.o
read_setupfile.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
seconds_since_epoch.o: mod_precision.o
set_ownership_ranges.o: mod_constants.o mod_woodbury_matrices.o
set_timestepping_coefficients.o: mod_options.o mod_precision.o mod_timestepping.o
setup_2D_schur.o: mod_constants.o mod_mesh_deformation_maps.o mod_precision.o mod_woodbury_matrices.o
setup_3D_deflation.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
setup_3D_poisson_solver.o: mod_constants.o mod_precision.o mod_transverse.o mod_woodbury_matrices.o
setup_3D_schur_preconditioner.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
setup_A_matrix.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
setup_capacitance_preconditioner.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
setup_cfl_calculation.o: mod_constants.o mod_legendre.o mod_mesh_deformation_maps.o mod_precision.o mod_woodbury_matrices.o
setup_deflation.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
setup_deformation_maps.o: mod_constants.o mod_legendre.o mod_options.o mod_precision.o
setup_deformed_derivatives.o: mod_constants.o mod_geom.o mod_legendre.o mod_mesh_deformation_maps.o mod_precision.o mod_woodbury_matrices.o
setup_filter_xz.o: mod_precision.o
setup_filter_y.o: mod_constants.o mod_options.o mod_precision.o mod_transverse.o
setup_initial_conditions.o: mod_field_variables.o mod_options.o mod_timestepping.o
setup_nullspace_projection.o: mod_constants.o mod_legendre.o mod_mesh_deformation_maps.o mod_nullspace.o mod_precision.o mod_woodbury_matrices.o
setup_poisson_operators.o: mod_precision.o mod_woodbury_matrices.o
setup_poisson_solver.o: mod_constants.o mod_mesh_deformation_maps.o mod_options.o mod_precision.o mod_woodbury_matrices.o
setup_sponge_layer.o: mod_constants.o mod_geom.o mod_parallel_linear_algebra.o mod_precision.o mod_sponge.o mod_woodbury_matrices.o
setup_transverse.o: mod_constants.o mod_precision.o mod_transverse.o
solve_3D_A.o: mod_constants.o mod_precision.o mod_transverse.o mod_woodbury_matrices.o
solve_3D_diffusion.o: mod_constants.o mod_errors.o mod_options.o mod_parallel_linear_algebra.o mod_precision.o mod_timestepping.o mod_woodbury_matrices.o
solve_3D_pressure.o: mod_constants.o mod_options.o mod_precision.o mod_parallel_linear_algebra.o mod_timestepping.o mod_woodbury_matrices.o
solve_3D_poisson_iteratively.o: mod_constants.o mod_options.o mod_parallel_linear_algebra.o mod_precision.o mod_timestepping.o mod_woodbury_matrices.o
solve_3D_viscous.o: mod_constants.o mod_errors.o mod_field_variables.o mod_options.o mod_parallel_linear_algebra.o mod_precision.o mod_timestepping.o mod_woodbury_matrices.o
solve_A.o: mod_precision.o mod_woodbury_matrices.o
solve_tridiagonal_system.o: mod_precision.o
synchronization.o: mod_constants.o mod_precision.o
validate_advection_terms.o: mod_constants.o mod_field_variables.o mod_geom.o mod_parallel_linear_algebra.o mod_precision.o mod_transverse.o mod_woodbury_matrices.o
validate_complex_gmres_3D.o: mod_constants.o mod_geom.o mod_options.o mod_parallel_linear_algebra.o mod_precision.o mod_timestepping.o mod_transverse.o mod_woodbury_matrices.o
validate_divergence.o: mod_constants.o mod_field_variables.o mod_geom.o mod_parallel_linear_algebra.o mod_precision.o mod_transverse.o mod_woodbury_matrices.o
validate_double_curl.o: mod_constants.o mod_field_variables.o mod_geom.o mod_parallel_linear_algebra.o mod_precision.o mod_transverse.o mod_woodbury_matrices.o
validate_fftw.o: mod_constants.o mod_parallel_linear_algebra.o mod_precision.o mod_woodbury_matrices.o
validate_filter_y.o: mod_constants.o mod_parallel_linear_algebra.o mod_precision.o mod_transverse.o mod_woodbury_matrices.o
validate_gradients.o: mod_constants.o mod_geom.o mod_parallel_linear_algebra.o mod_precision.o mod_transverse.o mod_woodbury_matrices.o
validate_helmholtz_2d_scalar_dirichlet.o: mod_constants.o mod_geom.o mod_options.o mod_parallel_linear_algebra.o mod_precision.o mod_timestepping.o mod_woodbury_matrices.o
validate_helmholtz_3d_scalar_dirichlet.o: mod_constants.o mod_geom.o mod_parallel_linear_algebra.o mod_precision.o mod_timestepping.o mod_transverse.o mod_woodbury_matrices.o
validate_helmholtz_3d_scalar_neumann.o: mod_constants.o mod_geom.o mod_parallel_linear_algebra.o mod_precision.o mod_timestepping.o mod_transverse.o mod_woodbury_matrices.o
validate_helmholtz_3d_vector_dirichlet.o: mod_constants.o mod_geom.o mod_parallel_linear_algebra.o mod_precision.o mod_timestepping.o mod_transverse.o mod_woodbury_matrices.o
validate_helmholtz_3d_vector_neumann.o: mod_constants.o mod_geom.o mod_parallel_linear_algebra.o mod_precision.o mod_timestepping.o mod_transverse.o mod_woodbury_matrices.o
validate_laplacian.o: mod_constants.o mod_geom.o mod_parallel_linear_algebra.o mod_precision.o mod_transverse.o mod_woodbury_matrices.o
validate_poisson_3d_iteratively.o: mod_constants.o mod_geom.o mod_parallel_linear_algebra.o mod_precision.o mod_transverse.o mod_woodbury_matrices.o
validate_real_valued_gmres_3D.o: mod_constants.o mod_geom.o mod_options.o mod_parallel_linear_algebra.o mod_precision.o mod_timestepping.o mod_transverse.o mod_validation.o mod_woodbury_matrices.o
validate_transport_term.o: mod_constants.o mod_field_variables.o mod_geom.o mod_parallel_linear_algebra.o mod_precision.o mod_transverse.o mod_woodbury_matrices.o
validate_viscous_solve_global.o: mod_constants.o mod_geom.o mod_parallel_linear_algebra.o mod_precision.o mod_timestepping.o mod_transverse.o mod_woodbury_matrices.o
validate_viscous_solve_global_real_valued.o: mod_constants.o mod_field_variables.o mod_geom.o mod_options.o mod_parallel_linear_algebra.o mod_precision.o mod_timestepping.o mod_transverse.o mod_validation.o mod_woodbury_matrices.o
validate_viscous_solve_single.o: mod_constants.o mod_field_variables.o mod_geom.o mod_options.o mod_parallel_linear_algebra.o mod_precision.o mod_timestepping.o mod_transverse.o mod_woodbury_matrices.o
write_fieldfile.o: mod_constants.o mod_field_variables.o mod_geom.o mod_io_field.o mod_precision.o mod_timestepping.o mod_transverse.o mod_woodbury_matrices.o
write_diagnosticsfile.o: mod_constants.o mod_errors.o mod_io_diagnostics.o mod_options.o mod_precision.o mod_timestepping.o mod_woodbury_matrices.o
write_postfile.o: mod_constants.o mod_field_variables.o mod_geom.o mod_io_post.o mod_precision.o mod_transverse.o mod_woodbury_matrices.o
write_postprocessor_xdmf.o: mod_constants.o mod_io_post.o
write_restartfile.o: mod_constants.o mod_field_variables.o mod_geom.o mod_io_restart.o mod_precision.o mod_timestepping.o mod_transverse.o
write_setupfile.o: mod_constants.o mod_precision.o mod_woodbury_matrices.o
write_xdmf.o: mod_constants.o mod_io_field.o mod_io_restart.o
