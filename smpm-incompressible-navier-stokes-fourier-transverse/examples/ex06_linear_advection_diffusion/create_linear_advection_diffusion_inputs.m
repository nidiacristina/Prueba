function [input_file_name, inputs] = create_linear_advection_diffusion_inputs()
% [input_file_name, inputs] = create_linear_advection_diffusion_inputs()
%
% Generates the input and initial condition files to run a small 3D Gaussian
% Bump Problem.  Both files are created within the current directory.
%
% Takes no arguments.
%
% Returns 2 values:
%
%   input_file_name - Full path to the input file name created.
%   inputs          - Structure containing all of the parameters written to
%                     input_file_name.
%

run_name = 'linear_advection_diffusion';

% construct the input and initial conditions file names.
input_file_name = sprintf( '%s/%s_in', ...
                           pwd, run_name );
init_file_name  = sprintf( '%s_init.h5', run_name );

% Set some constants related to the grid.
n     = 16;    % The number of GLL points per direction per element.
nsubx = 4;    % The number of x elements.
nsuby = 128;     % The number of y elements.
nsubz = 4;    % The number of z elements.
Lx    = 1.0;   % The length of the domain in the x direction.
Ly    = 1.0;   % The length of the domain in the y direction.
Lz    = 1.0;   % The length of the domain in the z direction.

% Build the mesh using the SMPM cartesian mesh builder.
[x, z] = smpm_build_cartesian_mesh( n, nsubx, nsubz, [0, Lx], [0, Lz] );

% Extrude the mesh into three dimensions.
dy = (Ly - 0 )/(nsuby);
y = 0 + dy*(0:nsuby-1);
%y           = linspace( 0, Ly, nsuby );
[x, y, z] = smpm_extrude_mesh( n, nsubx, nsubz, x, y, z );

% Set the reference density.  We assume a cavity filled with fresh water.
rho_0 = 1000;

% Construct Gaussian Function in the scalar field
% Set mean
mu = 0.5;

% Set standard deviation
sigma = 0.1;

% Construct Gaussian Bump with amplitude of 1.
rho0 = exp( - ( (x - mu).^2 ./ ( 2* sigma.^2) +...
    (y - mu).^2 ./ ( 2* sigma.^2) + (z - mu).^2 ./ ( 2* sigma.^2) )) ; 

% Set the constant velocity in the transverse
uy0 = 1.0;

% Set the remaining initial field values to zero.
ux0        = 0;
uz0        = 0;
ubc0       = 0;
dubcdz0    = 0;
rho0_bar   = 0;
rho0_bar_z = 0;

% Write the initial conditions file.
smpm_write_initfile( n, nsubx, nsuby, nsubz, x, y, z, rho0, ux0, uy0, uz0, ubc0, dubcdz0, rho0_bar, ...
                     rho0_bar_z, init_file_name );


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
inputs.tend                           = 20.0;
inputs.dt                             = 0.001; % A sample value

% Set the multiplicative factors on the penalty terms.
inputs.facrobin                       = 10000;
inputs.facrobin_ppe                   = 10000;

% Set the viscosity (nu) and the diffusivity (nu_d).
inputs.nu                             = 1e-3;
inputs.nu_d                           = 1e-3;

% Set the reference density.
inputs.rho_0                          = rho_0;

% Set the order of the exponential filter.
inputs.filter_order_xz                = 11;
inputs.filter_order_y                 = 11;

% Set the GMRES tolerance and maximum iterations for the poisson and viscous
% solvers.
inputs.gmres_tol_poisson              = 1e-8;
inputs.gmres_tol_viscous              = 1e-8;
inputs.gmres_maxit_poisson            = 1000;
inputs.gmres_maxit_viscous            = 1000;

% Set the boundary conditions.  1 means Dirichlet and 2 means Neumann.  The
% 4x1 array of boundary conditions represents the [bottom, right, top, left]
% boundaries.
inputs.bc_viscous                     = [2, 2, 2, 2];
inputs.bc_diffusion                   = [2, 2, 2, 2];

% Set some execution options; all of these are optional.
inputs.check_null_error               = logical( 0 );
inputs.check_numerical_error          = logical( 1 );
inputs.read_bcs_from_initfile         = logical( 0 );

inputs.do_interfacial_averaging       = logical( 1 );
inputs.use_capacitance_preconditioner = logical( 1 );
inputs.use_deflation                  = logical( 1 );

% Set the portions of INSE we want to solve
inputs.solve_momentum_equation        = logical( 0 );
inputs.use_gravity_force              = logical( 0 );
inputs.use_scalar_diffusion           = logical( 1 );
inputs.use_scalar_transport           = logical( 1 );

% Specify that you want to write out data every 1/4th of a second.
inputs.timesteps_between_writes       = round(0.25/inputs.dt);
inputs.timesteps_between_logs         = round(0.25/inputs.dt);

% Write the input file to disk using the SMPM input file writer.
smpm_write_inputfile( input_file_name, inputs );

end
