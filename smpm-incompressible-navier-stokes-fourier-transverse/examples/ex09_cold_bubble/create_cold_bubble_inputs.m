function [input_file_name, inputs] = create_cold_bubble_inputs(varargin)
% [input_file_name, inputs] = create_vortex_dippole_collision_inputs()
%
% Generates the input and initial conditions files to run a stable,
% vortex dippole collision based on the sample provided by Clercx & Bruneau.
%
%
% Takes none or 4 arguments.
%   n              - Number of GLL points per direction, per subdomain.
%   nsubx             - Number of subdomains in the x-direction.
%   nsuby             - Number of transverse grid points.
%   nsubz             - Number of subdomains in the z-direction.
%
% Returns 2 values:
%
%   input_file_name - Full path to the input file name created.
%   inputs          - Structure containing all of the parameters written to
%                     input_file_name.
%

run_name = 'cold_bubble';

% construct the input and initial conditions file names.
input_file_name = sprintf( '%s/%s_in', ...
                           pwd, run_name );
init_file_name  = sprintf( '%s_init.h5', run_name );

if nargin == 4
n     = varargin{1};
nsubx = varargin{2};
nsuby = varargin{3};
nsubz = varargin{4};
else    
% Set some constants related to the grid.
n     = 10;     % The number of GLL points per direction per element.
nsubx = 24;    % The number of x elements.
nsuby = 1;     % The number of y elements.
nsubz = 24;    % The number of z elements.

end

Lx    = 19.2;   % The length of the domain in the x direction.
Ly    = 1.0;   % The length of the domain in the y direction.
Lz    = 4.8;   % The length of the domain in the z direction.

%%%%%%

% Build the mesh using the SMPM cartesian mesh builder.
[x, z] = smpm_build_cartesian_mesh( n, nsubx, nsubz, [0, Lx], [0, Lz] );

% Extrude the mesh into three dimensions.
dy = (Ly - 0 )/(nsuby);
y = 0 + dy*(0:nsuby-1);
%y           = linspace( 0, Ly, nsuby );
[x, y, z] = smpm_extrude_mesh( n, nsubx, nsubz, x, y, z );

% Set the reference density.  We assume a cavity filled with fresh water.
rho_0 = 1000;

% Build the initial density profile.
%
% This is the initial function that characterizes the density
%Parameter used to calculate initial velocity 
xc=0.0;
xr=4.0;
zc=3.0;
zr=2.0;
z2=flipud(z);
z2=z;
L=(((x-xc)*xr^-1).^2+((z2-zc)*zr^-1).^2).^0.5;
Dt=zeros(size(L));
Dt(L<=1)=-15*(cos(pi*L(L<=1))+1)./2;
Dt(L>1)=0;

%alpha=26.85;
alpha=1;
rho0      = rho_0*(1-alpha*Dt);

% Set the remaining initial field values to zero.
%
% NOTE: Notice that you can set them to scalar values and the code will
%       correctly interpret them as a constant vector.
ux0        = 0;
uy0        = 0;
uz0        = 0;
ubc0       = 0;
dubcdz0    = 0;
rho0_bar   = 0;
rho0_bar_z = 0;

%%%%%%

% Write the initial conditions file to disk.
smpm_write_initfile( n, nsubx, nsuby, nsubz, x, y, z, rho0, ux0, uy0, uz0, ubc0, dubcdz0, ...,
                     rho0_bar, rho0_bar_z, init_file_name );

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
inputs.tend                           = 2;
inputs.dt                             = 1e-4;

% Set the multiplicative factors on the penalty terms.
inputs.facrobin                       = 10000;
inputs.facrobin_ppe                   = 100;

% Set the viscosity (nu) and the diffusivity (nu_d).
%inputs.nu                             = 0.0004;
inputs.nu                             = 0.0002;
inputs.nu_d                           = 0.0002;

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
inputs.timesteps_between_writes       = round((inputs.tend/inputs.dt)/60,0);
%inputs.timesteps_between_writes       = 1;
inputs.timesteps_between_logs         = round((inputs.tend/inputs.dt)/25,0);
inputs.timesteps_between_restarts     = round((inputs.tend/inputs.dt)/10);

% Write the input file to disk using the SMPM input file writer.
smpm_write_inputfile( input_file_name, inputs );

end
