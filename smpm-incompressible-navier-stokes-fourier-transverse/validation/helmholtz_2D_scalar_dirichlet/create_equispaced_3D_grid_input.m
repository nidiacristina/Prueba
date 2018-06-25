function [input_file_name, inputs] = create_equispaced_3D_grid_input()
% [input_file_name, inputs] = create_create_equispaced_3D_grid_input()
%
% Generates the input and initial condition files to run a small 3D box with an
% equi-spaced grid.  Both files are created within the current directory.
%
% Takes no arguments.
%
% Returns 2 values:
%
%   input_file_name - Full path to the input file name created.
%   inputs          - Structure containing all of the parameters written to
%                     input_file_name.
%

run_name = 'equispaced_3D_grid';

% construct the input and initial conditions file names.
input_file_name = sprintf( '%s/%s_in', ...
                           pwd, run_name );
init_file_name  = sprintf( '%s_init.h5', run_name );

% Set some constants related to the grid.
n     = 8;     % The number of GLL points per direction per element.
nsubx = 8;    % The number of x elements.
nsuby = 1;    % The number of y elements.
nsubz = 8;    % The number of z elements.
Lx    = 1.0;   % The length of the domain in the x direction.
Ly    = 1.0;   % The length of the domain in the y direction.
Lz    = 1.0;   % The length of the domain in the z direction.

% Build the mesh using the SMPM cartesian mesh builder.
[x, z] = smpm_build_cartesian_mesh( n, nsubx, nsubz, [0, Lx], [0, Lz] );

% Extrude the mesh into three dimensions.
y         = linspace( 0, Ly, nsuby );
[x, y, z] = smpm_extrude_mesh( n, nsubx, nsubz, x, y, z );

% Set the reference density.  We assume fresh water.
rho_0 = 1000.;

% Set the initial conditions for the stationary box problem.
ux0          = 0;         % The x-velocity is zero throughout the domain.
uy0          = 0;         % The y-velocity is zero throughout the domain.
uz0          = 0;         % The z-velocity is zero throughout the domain.
ubc0         = 0;         % The background current is zero throughout the domain.
dubcdz0      = 0;         % The background current's derivative is zero throughout the domain.
rho0         = 0;         % The initial Boussinesq density is zero everywhere.
rho0_bar     = 0;         % There is no background stratified density.
rho0_bar_dz  = 0;         % There background stratified density's derivative is zero throughout the domain.
s            = n * nsubz; % The number of grid points in the vertical direction.

% Write the initial conditions file to disk.
smpm_write_initfile( n, nsubx, nsuby, nsubz, x, y, z, rho0, ux0, uy0, uz0, ubc0, dubcdz0, rho0_bar, rho0_bar_dz, ...
                     init_file_name );

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
inputs.tend                           = 0.03;
inputs.dt                             = 3e-4;

% Set the multiplicative factors on the penalty terms.
inputs.facrobin                       = 1000;
inputs.facrobin_ppe                   = 1;

% Set the viscosity (nu) and the diffusivity (nu_d).
inputs.nu                             = 0.0040;
inputs.nu_d                           = 0.0040;

% Set the reference density.
inputs.rho_0                          = rho_0;

% Set the order of the exponential filter.
inputs.filter_order_xz                = 11;
inputs.filter_order_y                 = 11;

% Set the GMRES tolerance and maximum iterations for the poisson and viscous
% solvers.
inputs.gmres_tol_poisson              = 1e-6;
inputs.gmres_tol_viscous              = 1e-6;
inputs.gmres_maxit_poisson            = 1000;
inputs.gmres_maxit_viscous            = 200;

% Set the boundary conditions.  1 means Dirichlet and 2 means Neumann.  The
% 4x1 array of boundary conditions represents the [bottom, right, top, left]
% boundaries.
inputs.bc_viscous                     = [1, 1, 1, 1];
inputs.bc_diffusion                   = [1, 1, 1, 1];

% Set some execution options; all of these are optional.
inputs.check_null_error               = logical( 1 );
inputs.check_numerical_error          = logical( 1 );

% This option tells the code that the values of the initial conditions on the
% boundaries are meant to be read as persistent boundary conditions.  For the
% lid-driven cavity problem this means that the lid keeps translating to the
% right for t > 0.
inputs.read_bcs_from_initfile         = logical( 1 );

inputs.use_capacitance_preconditioner = logical( 1 );
inputs.use_deflation                  = logical( 1 );

% Specify that you want to write out every time-step of data.
inputs.timesteps_between_writes       = 1;
inputs.timesteps_between_logs         = 1;

% Write the input file to disk using the SMPM input file writer.
smpm_write_inputfile( input_file_name, inputs );

end
