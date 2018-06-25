function [input_file_name, inputs] = create_tankscale_wave_with_bg_inputs()
% [input_file_name, inputs] = create_tankscale_wave_with_bg_inputs()
%
% Generates the input and initial conditions files to run a stable,
% tank-scale shoaling internal wave with a background current based 
% on the sample provided in DJLES by Michael Dunphy.
%
% Takes no arguments.
%
% Returns 2 values:
%
%   input_file_name - Full path to the input file name created.
%   inputs          - Structure containing all of the parameters written to
%                     input_file_name.
%

run_name = 'tankscale_wave_with_bg';

% construct the input and initial conditions file names.
input_file_name = sprintf( '%s/%s_in', ...
                           pwd, run_name );
init_file_name  = sprintf( '%s_init.h5', run_name );

% Set some constants related to the grid. All lengths in meters.
n     = 10;     % The number of GLL points per direction per element.
nsubx = 80;     % The number of x elements.
nsuby = 1;      % The number of y elements.
nsubz = 10;     % The number of z elements.
Lx    = 24;     % The length of the domain in the x direction.
Ly    = 1;      % The length of the domain in the y direction.
Lz    = 0.1996; % The length of the domain in the z direction.

% Build the mesh using the SMPM cartesian mesh builder.
[xm, zm] = smpm_build_cartesian_mesh( n, nsubx, nsubz, [0, Lx], [-Lz, 0] );

% Extrude the mesh into three dimensions.
dy = (Ly - 0 )/nsuby;
ym = 0 + dy*(0:nsuby-1);
[xm, ym, zm] = smpm_extrude_mesh( n, nsubx, nsubz, xm, ym, zm );

% Get some secondary constants for rough computations with domain size.
Lx = max(xm(:)) - min(xm(:));
Lz = max(zm(:)) - min(zm(:));

% Load the breaking wave data for building an initial conditions file.
n_mine = n;
load( 'data/data.mat' );
n = n_mine;

% Set the reference density.  We assume fresh water.
rho_0 = 1000.;

% Grab the Background Velocity Values
ubc_DJL = ubc;
dubcdz_DJL = dubcdz;

% Determine domain length of DJL
LDJL_MAX = max(X(:));
LDJL_MIN = min(X(:));

% Extend the fields for interpolation.
% Err, scratch that.  The fields should be zero outside of
% the region they inhabit.
nx = size(X,2);
nz = size(X,1);
        
% Read In Full Density Profile
rho_full = density.*rho_0;

% Grab the Background Density
rho_bar = repmat(rho_full(:,1) - rho_0, 1, size(rho_full,2) );

% Compute derivative of background density profile
[~,rho_bar_z] = gradient(rho_bar,X,Z);

% Compute Perturbation Density
rho_p = (rho_full - rho_0) - rho_bar;
 
% Next, window these functions so that they are zero at the boundaries to
% avoid waves.
W = tukeywin( nx, 0.5 );
W = repmat( W', nz, 1 );
rho_p    = rho_p .* W;
u_full   = uwave .* W;
w_full   = w .* W;

% Interpolate with splines, extrapolating only zero values.
% CAREFUL WITH EXTRAPOLATING: SUMEDH HAD IT AND IT SCREWS UP THE 
% VELOCITY FIELD
rho0 = interp2( X, Z,rho_p, xm, zm, 'spline' );
ux0  = interp2( X, Z,u_full, xm, zm, 'spline' );
uy0  = 0;
uz0  = interp2( X, Z,w_full, xm, zm, 'spline' );
        
% Interpolate Background Velocity Variables
ubc0 = zeros(size(xm(:)));
dubcdz0 = zeros(size(xm(:)));
        
% Zero the field outside the DJL Domain. 
% Note: the interpolation produces error near the boundaries of the
% computational domain and thus without the extrapolation option, we
% need to do this. 
ux0(find(xm > LDJL_MAX))  = 0;
uz0(find(xm > LDJL_MAX))  = 0;
rho0(find(xm > LDJL_MAX)) = 0;
ux0(find(xm < LDJL_MIN))  = 0;
uz0(find(xm < LDJL_MIN))  = 0;
rho0(find(xm < LDJL_MIN)) = 0;

% Add a noise floor, eliminating noise in the interpolation and computation
% of the DJL solution.
minlevel = 1e-8;
rho0( abs( rho0 ) < minlevel ) = 0.0;
ux0( abs( ux0 )  < minlevel )  = 0.0;
uz0( abs( uz0 )  < minlevel )  = 0.0;

% Reshape Grid for generating background fields
xm = xm(:);
zm = zm(:);

% Build the background density field for Boussinesq.
z_bar        = Z(:,1);
z_bar        = [ 10; z_bar; -1e6];
dummyrho     = rho_bar(:,1);
dummyrho_z   = rho_bar_z(:,1);
dummyubc     = ubc_DJL(:,1);
dummydubcdz  = dubcdz_DJL(:,1);
dummyrho     = [ dummyrho(1); dummyrho; dummyrho(end) ];
dummyrho_z   = [ dummyrho_z(1); dummyrho_z; dummyrho_z(end) ];
dummyubc     = [ dummyubc(1); dummyubc; dummyubc(end) ];
dummydubcdz  = [ dummydubcdz(1); dummydubcdz; dummydubcdz(end) ];
rho0_bar     = zeros(size(rho0));
rho0_bar_z   = zeros(size(rho0));
ubc0         = zeros(size(rho0));
dubcdz0      = zeros(size(rho0));
for ii = 1:n*nsubx
    iistart = ( ii - 1 ) * nsubz * n + 1;
    iiend   = ( ii - 1 ) * nsubz * n + nsubz * n;
    rho0_bar( iistart:iiend )   = interp1( z_bar, dummyrho, zm(iistart:iiend), 'spline' );
    rho0_bar_z( iistart:iiend ) = interp1( z_bar, dummyrho_z, zm(iistart:iiend), 'spline' );
    ubc0( iistart:iiend )       = interp1( z_bar, dummyubc, zm(iistart:iiend), 'spline' );
    dubcdz0( iistart:iiend )    = interp1( z_bar, dummydubcdz, zm(iistart:iiend), 'spline' );
end

% replicate our fields across our Y dimension.
ux0      = repmat( ux0,      [nsuby, 1] );
uz0      = repmat( uz0,      [nsuby, 1] );
rho0     = repmat( rho0,     [nsuby, 1] );
rho0_bar = repmat( rho0_bar, [nsuby, 1] );
ubc0     = repmat( ubc0,     [nsuby, 1] );
dubcdz0  = repmat( dubcdz0,  [nsuby, 1] );

% Write the mesh and init files to disk.
smpm_write_initfile( n, nsubx, nsuby, nsubz, xm, ym, zm, ...
                     rho0, ux0, uy0, uz0, ubc0, dubcdz0, rho0_bar, ...
                     rho0_bar_z, init_file_name );

% Estimate the maximum time-step possible.
dx     = abs( diff( xm(:) ) );
dz     = abs( diff( zm(:) ) );
dx_min = min( dx(dx > 1e-6) );
dz_min = min( dz(dz > 1e-6) );
dt_x   = dx_min / max( abs( ux0(:) ) );
dt_z   = dz_min / max( abs( uz0(:) ) );
dt     = min( dt_x, dt_z ) / 6;

% Set the maximum time.
c_wave = c;
L      = max( xm(:) ) - min( xm(:) );
T      = L / c_wave;

% Sponge Layer time scale
eps = 1e-6;
time_scale_sponge = 2*c/(.05*L + .05*L)*(-1)*log(eps);

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
inputs.facrobin                       = 1000;
inputs.facrobin_ppe                   = 1000;

% Set the viscosity (nu) and the diffusivity (nu_d).
%
% Inviscid mode (nu = 0)
inputs.nu                             = 0;
inputs.nu_d                           = 0;

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

% Specify that you want to write out data and see a status update every 1s
inputs.timesteps_between_logs         = round(0.5/inputs.dt);
inputs.timesteps_between_writes       = round(1.0/inputs.dt);
inputs.timesteps_between_restarts     = round(10.0/inputs.dt);

% Add sponge layer parameters
inputs.apply_sponge_layer             = logical( 1 );
inputs.sponge_layer_location          = [0 1 0 1];
inputs.time_scale_sponge_in_x         = time_scale_sponge;

% Add adaptive timestep
inputs.adaptive_timestep              = logical( 1 );

% Write the input file to disk using the SMPM input file writer.
smpm_write_inputfile( input_file_name, inputs );

end
