function [input_file_name, inputs] = create_tankscale_breaking_wave_inputs()
% [input_file_name, inputs] = create_tankscale_breaking_wave_inputs()
%
% Generates the input and initial conditions files to run a breaking,
% tank-scale shoaling internal wave based on South China Sea bathymetry near
% Dongsha Island using initial conditions generated by the DJL solver code.
%
% Takes no arguments.
%
% Returns 2 values:
%
%   input_file_name - Full path to the input file name created.
%   inputs          - Structure containing all of the parameters written to
%                     input_file_name.
%

run_name = 'tankscale_breaking_wave';

% construct the input and initial conditions file names.
input_file_name = sprintf( '%s/%s_in', ...
                           pwd, run_name );
init_file_name  = sprintf( '%s_init.h5', run_name );

% Set some constants related to the grid.
n     = 15;   % The number of GLL points per direction per element.
nsubx = 100;  % The number of x elements.
nsuby = 1;    % The number of y elements.
nsubz = 10;   % The number of z elements.
Lx    = 8;    % The length of the domain in the x direction.
Ly    = 1;    % The length of the domain in the y direction.
Lz    = 0.35; % The length of the domain in the z direction.

% Build the mesh using the SMPM cartesian mesh builder.
[xm, zm] = smpm_build_cartesian_mesh( n, nsubx, nsubz, [2, 2 + Lx], [-Lz, 0] );

% Set the reference density.  We assume fresh water.
rho_0 = 1000.;

% Load the breaking wave data for building an initial conditions file.
n_mine = n;
load( sprintf( '%s/data/tankscale_wave', pwd ) );
n = n_mine;

% Extend the fields for interpolation.
nx = size( x_full, 2 );
nz = size( x_full, 1 );

x2   = [-1e6 * ones( nz, 1 ) x_full 1e6*ones( nz, 1 )];
z2   = [z_full(:, 1), z_full, z_full(:, end) ];
rho2 = [rho_full(:, 1), rho_full, rho_full(:, end)];
ux2  = [u_full(:, 1), u_full, u_full(:, end)];
uz2  = [w_full(:, 1), w_full, w_full(:, end)];

x2   = [x2(1, :); x2; x2(end, :)];
z2   = [10 * ones( 1, nx+2 ); z2; -1e3*ones( 1, nx+2 )];
rho2 = [rho2(1, :); rho2; rho2(end, :)];
ux2  = [ux2(1, :);  ux2;  ux2(end, :)];
uz2  = [uz2(1, :);  uz2;  uz2(end, :)];

% Interpolate the fields onto the meshed grid.
rho0 = 1000 * interp2( x2, z2, rho2, xm, zm );
ux0  = interp2( x2, z2, ux2,  xm, zm );
uy0  = 0;
uz0  = interp2( x2, z2, uz2,  xm, zm );
p0   = 0 * rho0;

% Build the background density field for Boussinesq.
rho_bar    = 1000 * rho_full(:, 1) - rho_0;
rho_bar_z  = gradient(rho_bar,z_full(:,1)); 
z_bar      = z_full(:, 1);
z_bar      = [10; z_bar; -1e6];
rho_bar    = [rho_bar(1); rho_bar; rho_bar(end)];
rho_bar_z  = [rho_bar_z(1); rho_bar_z; rho_bar_z(end)];
rho0_bar   = zeros(size(rho0));
rho0_bar_z = zeros(size(rho0));
for ii = 1:n*nsubx
    iistart                   = (ii - 1) * nsubz * n + 1;
    iiend                     = (ii - 1) * nsubz * n + nsubz * n;
    rho0_bar(iistart:iiend)   = interp1( z_bar, rho_bar, zm(iistart:iiend) );
    rho0_bar_z(iistart:iiend) = interp1( z_bar, rho_bar_z, zm(iistart:iiend) );
end

% replicate our fields across our Y dimension.
ux0        = repmat( ux0,        [nsuby, 1] );
uz0        = repmat( uz0,        [nsuby, 1] );
rho0       = repmat( rho0,       [nsuby, 1] );
rho0_bar   = repmat( rho0_bar,   [nsuby, 1] );
rho0_bar_z = repmat( rho0_bar_z, [nsuby, 1] );

% Set the background currents to zero.
ubc0    = zeros(size(ux0));
dubcdz0 = zeros(size(ux0));

% Artificially increase the velocities to initialize wave breaking.
ux0 = ux0 * 1.25;
uz0 = uz0 * 1.25;

% Extrude the mesh into three dimensions.
dy = (Ly - 0 )/(nsuby);
ym = 0 + dy*(0:nsuby-1);
%ym           = linspace( 0, Ly, nsuby );
[xm, ym, zm] = smpm_extrude_mesh( n, nsubx, nsubz, xm, ym, zm );

% Write the mesh and init files to disk.
smpm_write_initfile( n, nsubx, nsuby, nsubz, xm, ym, zm, ...
                     rho0 - rho0_bar - rho_0, ux0, uy0, uz0, ubc0, dubcdz0, rho0_bar, ...
                     rho0_bar_z, init_file_name );

% Estimate the maximum time-step possible.
dx     = abs( diff( xm ) );
dz     = abs( diff( zm ) );
dx_min = min( dx(dx > 1e-6) );
dz_min = min( dz(dz > 1e-6) );
dt_x   = dx_min / max( abs( ux0 ) );
dt_z   = dz_min / max( abs( uz0 ) );
dt     = min( dt_x, dt_z ) / 3;

% Set the maximum time.
c_wave = c;
L      = max( xm(:) );
T      = L / c_wave;

%% Create a structure that will contain all the inputs for the SMPM code.
inputs                                = struct();

% Give this run a name and specify its initial conditions file.
inputs.fname_runname                  = run_name;
inputs.fname_init                     = init_file_name;

% Set the number of GLL points and the number of x, y, and z elements.
inputs.n                              = n;
inputs.nsubx                          = nsubx;
inputs.nsuby                          = nsuby;
inputs.nsubz                          = nsubz;

% Set the final time and the time step.
inputs.dt                             = dt;
inputs.tend                           = T;

% Set the multiplicative factors on the penalty terms.
inputs.facrobin                       = 10000;
inputs.facrobin_ppe                   = 10000;

% Set the viscosity (nu) and the diffusivity (nu_d).
%
% XXX: Should these be 1e-6?  They were 1e-3 in the last NLIW shoaling work SMJ
%      did.
inputs.nu                             = 1e-6;
inputs.nu_d                           = 1e-6;

% Set the reference density.
inputs.rho_0                          = rho_0;

% Set the order of the exponential filter.
inputs.filter_order_xz                = 11;
inputs.filter_order_y                 = 11;

% Set the GMRES tolerance and maximum iterations for the poisson and viscous
% solvers.
inputs.gmres_tol_poisson              = 1e-7;
inputs.gmres_tol_viscous              = 1e-7;
inputs.gmres_maxit_viscous            = 500;
inputs.gmres_maxit_poisson            = 3000;

% NOTE: Make sure that the restart parameter is set to the maximum number of
%       GMRES iterations so that we don't actually allow GMRES restarting.
inputs.gmres_restart_viscous          = 500;
inputs.gmres_restart_poisson          = 3000;

% Set the boundary conditions.  1 means Dirichlet and 2 means Neumann.  The
% 4x1 array of boundary conditions represents the [bottom, right, top, left]
% boundaries.
%
% Set the boundary conditions to be free-slip for the viscous terms,
% and Neumann for the diffusion operator.
inputs.bc_viscous                     = [2, 2, 2, 2];
inputs.bc_diffusion                   = [2, 2, 2, 2];

% Set some execution options; all of these are optional.
inputs.check_null_error               = logical( 0 );
inputs.check_numerical_error          = logical( 1 );
inputs.read_bcs_from_initfile         = logical( 0 );

inputs.do_interfacial_averaging       = logical( 1 );
inputs.use_capacitance_preconditioner = logical( 1 );
inputs.use_deflation                  = logical( 1 );

% Specify that you want to write out data and see a status update every 100th
% time-step.
inputs.timesteps_between_logs         = 100;
inputs.timesteps_between_writes       = 100;

% Write the input file to disk using the SMPM input file writer.
smpm_write_inputfile( input_file_name, inputs );

end
