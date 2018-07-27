function [Ens] = compute_enstrophy(data)
%function [ux_E,uy_E,uz_E] = compute_energy(data)
%
% Computes enstrophy of the flow field according to Clercx & Bruneau (2003)
% 
% Takes 1 input argument:
%   
%   data - Field structure, as read by smpm_read_fieldfile, containing
%           the field data (i.e. grid variables, field variables, time, 
%           etc...) for the run of interest.
%
% Returns 1 value:
%   
%   Ens     - Escalar resulting from enstrophy computation
%
%
% Jul 2018
% NCRG

% 1) Grab the field and grid properties
x   = data.grid.x;
y   = data.grid.y;
z   = data.grid.z;
n   = data.grid.n;
mx  = data.grid.mx;
my  = data.grid.my;
mz  = data.grid.mz;

ux  = data.field.ux + data.ic.ubc;
uy  = data.field.uy;
uz  = data.field.uz;

% 2) Compute vorticity
% Differentiate the density field
[~, ~, dux_z ] = smpm_compute_gradient( ux, n, mx, my, mz, x, y, z);
[duz_x, ~, ~ ] = smpm_compute_gradient( uz, n, mx, my, mz, x, y, z);
vorty=duz_x-dux_z;          % Compute vorticity

% 3) Compute the integral - Gauss quadrature
ff=vorty.^2;                 %Discrete function to be integrated
% 3.1) Integral in the x direction
[~, wi, ~] = lglnodes( double( n ) - 1 ); %GLL weight and points
a=x(1);                      %Limits of integration
b=x(end);                    %Limits of integration
w=((b-a)/2)*wi;                %Mapping weights
int_x=ff*repmat(w,mz,1);     %Compute the integral in x direction
% 3.2) Integral in the z direction
a=z(1);                      %Limits of integration
b=z(end);                    %Limits of integration
w=((b-a)/2)*wi;                %Mapping weights
int_z=int_x'*repmat(w,mz,1); %Compute the integral in z direction

Ens=0.5*int_z;


end

