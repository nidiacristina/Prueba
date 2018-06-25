% Example 07: Builds the input files for the tank-scale internal wave
% propagation problem with a background current. This is the same wave as
% in Michael Dunphy's DJLES case_ubg.m

% Add the SMPM Matlab directory to your path (just in case you haven't already
% done this).
addpath( genpath( '../../' ) );

[input_file_name, inputs] = create_tankscale_wave_with_bg_inputs();
init                      = smpm_read_initfile( inputs.fname_init );

% Visualize Initial Conditions
n     = inputs.n;
nsubx = inputs.nsubx;
nsubz = inputs.nsubz;
x     = init.grid.x(:, :, 1);
z     = init.grid.z(:, :, 1);
rho   = init.ic.rho(:, :, 1);
ux    = init.ic.ux(:, :, 1);
uz    = init.ic.uz(:, :, 1);

% Visualize the initial condition with the mesh overlaid.
contourf( reshape( x, nsubz * n, nsubx * n ), ...
          reshape( z, nsubz * n, nsubx * n ), ...
          reshape( rho, nsubz * n , nsubx * n ) );
smpm_visualize_mesh( n, nsubx, nsubz, x, z, gcf );
