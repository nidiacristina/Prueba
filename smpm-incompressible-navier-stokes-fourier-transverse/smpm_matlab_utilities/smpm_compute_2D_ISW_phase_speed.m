function [c] = smpm_compute_2D_ISW_phase_speed( data , input )
%[c] = smpm_compute_2D_ISW_phase_speed( data, input )
%
% Computes the phase speed of a Propagating Internal Solitary Wave via the
% method of Fringer & Street (2003), by minimizaing the velocity components
% normal to the lines of constant density. For more information, please
% refer to their work on Page 330.
%
% Takes 2 input argument:
%   
%   field - Field structure, as read by smpm_read_fieldfile, containing
%           the field data (i.e. grid variables, field variables, time, 
%           etc...) for the run of interest.
%
%   input - Input structure, as read by smpm_read_inputfile, containing all
%           the input information for the run of interest.
%
% Returns 1 value:
%   
%   c     - The wave phase speed in m/s.
%
% May 2018
% Gustavo Rivera

% Grab the field and grid properties
x   = data.grid.x;
y   = data.grid.y;
z   = data.grid.z;
n   = data.grid.n;
mx  = data.grid.mx;
my  = data.grid.my;
mz  = data.grid.mz;

% Grab Field Data
ux  = data.field.ux + data.ic.ubc;
uz  = data.field.uz;
rho = data.field.rho + data.ic.rho_bar + input.rho_0;

% Differentiate the density field
[rho_x, ~, rho_z] = smpm_compute_gradient( rho, n, mx, my, mz, x, y, z);

% Construct the integral kernels
numerator   = ( ux .* rho_x + uz .* rho_z ) .* rho_x;
denominator = rho_x.^2;

% Compute the phase speed
c = integrate_grid_function_2D(numerator(:), x(:), z(:), n, mx, mz) ./...
        integrate_grid_function_2D( denominator(:), x(:), z(:), n, mx, mz);

end

