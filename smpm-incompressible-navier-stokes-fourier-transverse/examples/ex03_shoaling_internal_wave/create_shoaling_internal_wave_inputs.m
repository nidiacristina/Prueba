function [input_file_name, inputs] = create_shoaling_internal_wave_inputs()
% [input_file_name, inputs] = create_shoaling_internal_wave_inputs()
%
% Generates the input and initial conditions files to run a shoaling internal
% wave over bathymetry near Dongsha Island in the South China Sea with initial
% conditions created using the DJL solver code.
%
% Takes no arguments.
%
% Returns 2 values:
%
%   input_file_name - Full path to the input file name created.
%   inputs          - Structure containing all of the parameters written to
%                     input_file_name.
%

run_name = 'dongsha';

% construct the input and initial conditions file names.
input_file_name = sprintf( '%s/%s_in', ...
                           pwd, run_name );
init_file_name  = sprintf( '%s_init.h5', run_name );

% Set some constants related to the grid.
n     = 15;
nsubx = 200;
nsuby = 1;
nsubz = 12;

% Set the constant fluid density of water. We assume warm(-ish) salt water in
% the ocean.
rho_0 = 1020.;

% Load the bathymetric data and create a slice for meshing.

% Load the bathymetric data.
load( sprintf( '%s/data/Dongsha_bathy.mat', pwd ) );

% Extract a slice in latitudes and longitudes.

% Lien_2012 track. Latitude is in Degrees N.  Longitude is in Degrees E.
%
% NOTE: these coordinates are just made up; the coordinates reported in the Lien2012
%       paper do not result in bathymetry that looks anything like what they show
%       in the paper.  These coordinates present a smooth-ish bathy that somewhat
%       resembles the Dongsha track presented in Lien2012.
lat0 = 21 - 4 / 60;
lat1 = 21 - 4 / 60;
lon0 = 117.5;
lon1 = 117.0;

% Fix the Depths and Lats and Lons that were just loaded from the field bathymetric data.
Depth2 = flipud( Depth' );
lat    = linspace( 20, 22, size( Depth2, 1 ) );
lon    = linspace( 116.5, 119, size( Depth2, 2 ) );

% Convert lats and longs to kilometers.
%
% NOTE: magic constant comes from Stack Overflow
%
%       http://stackoverflow.com/questions/639695/how-to-convert-latitude-or-longitude-to-meters
%
magic  = 0.000008998719243599958;
width  = sqrt( (max( lon(:) ) - min( lon(:) )).^2  ) / magic / 1000;
height = sqrt( (max( lat(:) ) - min( lat(:) )).^2  ) / magic / 1000;
x      = linspace( 0, 1, size( Depth2, 2 ) ) * width;
dy = (1 - 0 )/(size( Depth2, 1 ) );
y = 0 + dy*(0:size( Depth2, 1 )-1) .* height;
%y      = linspace( 0, 1, size( Depth2, 1 ) ) * height
z      = Depth2;
[X, Y] = meshgrid( x, y );

% Extract a bathymetric slice to test build a shoaling grid.

% Convert the interpolation track from lat/lon to x/z if available.
if ( exist( 'lon0', 'var' ) )
    x0 = interp1( lon, x, lon0, 'spline' );
    x1 = interp1( lon, x, lon1, 'spline' );
    y0 = interp1( lat, y, lat0, 'spline' );
    y1 = interp1( lat, y, lat1, 'spline' );
else
    x0 = 201;
    y0 = 210;
    x0 = 111;
    y0 = 120;
    x1 = 50.9;
    y1 = 81.9;
end

% Set up the interpolation parameters.
nint  = 200;
xi    = linspace( x0, x1, nint );
yi    = linspace( y0, y1, nint );
zi    = interp2( x, y, z, xi, yi, 'cubic' );

% Compute the along-track distance.
d = sqrt( (x1 - x0).^2 + (y1 - y0).^2 );
d = linspace( 0, d, nint );

% Add in a false plateau whose length is equal to the NLIW wavelength.
% This is done so the wave can propagate in a flat bottom for a short
% period of time at the beginning of the simulation.
n_old = n;
load( sprintf( '%s/data/Dongsha_DJL_solution.mat', pwd ) );

n              = n_old;
plateau_buffer = 0.05; % Fraction of minimum plateau size to add as a buffer.
plateau_length = (max( x_full(:) ) * (1.0 + plateau_buffer)) / 1000;
d              = [d(1), d + plateau_length];
zi             = [zi(1), zi];

% Build the bathymetric mesh with the bottom profile zi as extracted from the field data.
d        = d * 1000;
[xm, zm] = smpm_build_bathy_mesh( n, nsubx, nsubz, d, -zi );

% Get some secondary constants for rough computations with domain size.
Lx = max( xm ) - min( xm );
Ly = 100.0;
Lz = max( zm ) - min( zm );

% Extrude the mesh into three dimensions.
ym        = linspace( 0, Ly, nsuby );
[x, y, z] = smpm_extrude_mesh( n, nsubx, nsubz, xm, ym, zm );

% Load the Dongsha internal wave run for building an initial conditions file.
% This data was produced by the DJL code.

% Load the data.
n_mine = n;
load( sprintf( '%s/data/Dongsha_DJL_solution.mat', pwd ) );
n = n_mine;

% Extend the fields for interpolation.
% Err, scratch that.  The fields should be zero outside of
% the region they inhabit.
nx = size( x_full, 2 );
nz = size( x_full, 1 );

% Use an inclusion-like map to project the field from the DJL solver's
% grid to the SMPM code's grid.

% Step zero, convert density to standard units (standard defined by my
% SMPM code).
rho_full = rho_full * 1000;

% First, subtract the background density field to build the Boussinesq
% perturbed density (we'll interpolate onto this).
rho_bar  = rho_full(:, 1) - rho_0;
rho_full = rho_full - repmat( rho_bar, 1, size( rho_full, 2 ) ) - rho_0;

% Next, window these functions so that they are zero at the boundaries
% to avoid waves.
W        = tukeywin( nx, 0.5 );
W        = repmat( W', nz, 1 );
rho_full = rho_full .* W;
u_full   = u_full .* W;
w_full   = w_full .* W;

% Compute gradient of background density
rho_bar_z = gradient(rho_bar,z_full(:,1));

% Interpolate with splines, extrapolating only zero values.
rho0 = interp2( x_full, z_full, rho_full, xm, zm, 'cubic', 0.0 );
ux0  = interp2( x_full, z_full,   u_full, xm, zm, 'cubic', 0.0 );
uy0  = 0;
uz0  = interp2( x_full, z_full,   w_full, xm, zm, 'cubic', 0.0 );

% Add a noise floor, eliminating noise in the interpolation and
% computation of the DJL solution.
minlevel = 1e-8;
rho0( abs( rho0 ) < minlevel ) = 0.0;
ux0( abs( ux0 )  < minlevel ) = 0.0;
uz0( abs( uz0 )  < minlevel ) = 0.0;

% Build the background density field for Boussinesq.
z_bar       = z_full(:, 1);
z_bar       = [10; z_bar; -1e6];
rho_bar     = [rho_bar(1); rho_bar; rho_bar(end)];
rho_bar_z   = [rho_bar_z(1); rho_bar_z; rho_bar_z(end)];
rho0_bar    = zeros(size(rho0));
rho0_bar_z  = zeros(size(rho0));
for ii = 1:n*nsubx
    iistart = (ii - 1) * nsubz * n + 1;
    iiend   = (ii - 1) * nsubz * n + nsubz * n;
    rho0_bar(iistart:iiend)   = interp1( z_bar, rho_bar,   zm(iistart:iiend), 'spline' );
    rho0_bar_z(iistart:iiend) = interp1( z_bar, rho_bar_z, zm(iistart:iiend), 'spline' );
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

% Estimate the maximum time-step possible.
dx     = abs( diff( xm ) );
dz     = abs( diff( zm ) );
dx_min = min( dx( dx > 1e-6 ) );
dz_min = min( dz( dz > 1e-6 ) );
dt_x   = dx_min ./ max( abs( ux0(:) ) );
dt_z   = dz_min ./ max( abs( uz0(:) ) );
dt     = min( dt_x, dt_z ) / 2;

% Set the maximum time-step.
c_wave = 6.0;
c_wave = c;
L      = max( xm );
T      = L / c_wave;

% Write the initial condition file to disk.
smpm_write_initfile( n, nsubx, nsuby, nsubz, x, y, z, ...
                     rho0, ux0, uy0, uz0, ubc0, dubcdz0, rho0_bar, ...
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
inputs.dt                             = dt;
inputs.tend                           = T;

% Set the multiplicative factors on the penalty terms.
inputs.facrobin                       = 1;
inputs.facrobin_ppe                   = 10000;

% Set the viscosity (nu) and the diffusivity (nu_d).
%
% NOTE: These are of course not the actual values of ocean
%       viscosity/diffusivity.
inputs.nu                             = 1e-3;
inputs.nu_d                           = 1e-3;

% Set the reference density.
inputs.rho_0                          = rho_0;

% Set the order of the exponential filter.
inputs.filter_order_xz                = 11;
inputs.filter_order_y                 = 11;

% Set the GMRES tolerance and maximum iterations for the poisson and viscous
% solvers.
inputs.gmres_tol_poisson              = 1e-9;
inputs.gmres_tol_viscous              = 1e-9;
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

% Specify that you want to write out data every 100th time-step.
inputs.timesteps_between_writes       = 100;

% Write the input file to disk using the SMPM input file writer.
smpm_write_inputfile( input_file_name, inputs );

end
