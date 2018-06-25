% Example 02: Builds the inputs for the lock-exchange problem.  All files are
% written to the current directory.  A plot of the initial density is generated
% after input files are created.

% Add the SMPM Matlab directory to your path (just in case you haven't already
% done this).
addpath( genpath( '../../' ) );

[input_file_name, inputs] = create_lock_exchange_inputs();
init                      = smpm_read_initfile( inputs.fname_init );

% rename some of our longer variables.  we explicitly pull the first
% transverse slice to point out that we technically have 3D data.
n     = inputs.n;
nsubx = inputs.nsubx;
nsubz = inputs.nsubz;
x     = init.grid.x(:, :, 1);
z     = init.grid.z(:, :, 1);
rho   = init.ic.rho(:, :, 1);

% Visualize the initial condition with the mesh overlaid.
contourf( reshape( x, nsubz * n, nsubx * n ), ...
          reshape( z, nsubz * n, nsubx * n ), ...
          reshape( rho, nsubz * n , nsubx * n ) );
smpm_visualize_mesh( n, nsubx, nsubz, x, z, gcf );
   %
   % smpm_visualize_mesh will draw just the element boundaries as lines on any
   % MATLAB figure window you ask.

