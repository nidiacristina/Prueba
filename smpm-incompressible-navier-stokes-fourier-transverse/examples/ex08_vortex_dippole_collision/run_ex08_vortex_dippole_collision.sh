#!/bin/bash
#
# This is a script to run the lid driven cavity test case.

# July 2018.

# Set some file names.
fname_in="vortex_dippole_collision_in" # The name of the input file.
tag_out="vortex_dippole_collision"     # All files generated will have this as a prefix.

# Set the number of MPI ranks and threads.
nranks=16   # Number of MPI ranks.
nthreads=1  # Number of OMP threads (currently the code is not threaded so this parameter is unused).

### Nothing below this line should need to be changed by the user. ###

# Increase the stack size (for bigger runs we may have to overflow the 8kB stack).
ulimit -s unlimited

# Add to path.
export LD_LIBRARY_PATH=/usr/local/lib:/usr/lib

# Generate some file names.
fname_out=$tag_out"_out"         # The output file.
fname_err=$tag_out"_err"         # The error log file.
fname_log=$tag_out"_log"         # The log file.
fname_diag=$tag_out"_diag"       # The diagnostics file.
fname_restart=$tag_out"_restart" # The restart file.

# Set the number of threads to one (currently the SMPM code is not threaded).
export OMP_NUM_THREADS=1

# Run the SMPM code.
/usr/bin/mpiexec -np $nranks ../../smpm_incompressible_navier_stokes $fname_in $fname_out $fname_restart $fname_diag 2>$fname_err | tee $fname_log
