% Example 03: Builds the inputs for a shoaling nonlinear internal wave problem.
% All files are written to the current directory.  A plot of the internal wave's
% characteristics (velocities and density) is generated after input files are
% created.

% Add the SMPM Matlab directory to your path (just in case you haven't already
% done this).
addpath( genpath( '../../' ) );

[input_file_name, inputs] = create_shoaling_internal_wave_inputs();
init                      = smpm_read_initfile( inputs.fname_init );

% rename some of our longer variables.  we explicitly pull the first
% transverse slice to point out that we technically have 3D data.
n        = inputs.n;
nsubx    = inputs.nsubx;
nsubz    = inputs.nsubz;
x        = init.grid.x(:, :, 1);
z        = init.grid.z(:, :, 1);
ux0      = init.ic.ux(:, :, 1);
uz0      = init.ic.uz(:, :, 1);
ubc      = init.ic.ubc(:, :, 1);
rho0     = init.ic.rho(:, :, 1);
rho0_bar = init.ic.rho_bar(:, :, 1);

% Let's visualize this initial condition along with its mesh in three sub-figures.

% Open a MATLAB figure window.
figure;

% Plot the total density.
subplot( 3, 1, 1 );
contourf( x, z, rho0 + rho0_bar );
title( 'density' );
smpm_visualize_mesh( n, nsubx, nsubz, x, z, gcf );

% Plot the total density.
subplot( 3, 1, 2 );
contourf( x, z, ux0);
title( 'u-velocity' );
smpm_visualize_mesh( n, nsubx, nsubz, x, z, gcf );

% Plot the total density.
subplot( 3, 1, 3 );
contourf( x, z, uz0 );
title( 'w-velocity' );
smpm_visualize_mesh( n, nsubx, nsubz, x, z, gcf );
xlabel( 'meters' );

return
