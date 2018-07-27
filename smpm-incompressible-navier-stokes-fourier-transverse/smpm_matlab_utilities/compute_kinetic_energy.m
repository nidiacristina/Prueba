function [kE] = compute_kinetic_energy(data)
%function [ux_E,uy_E,uz_E] = compute_energy(data)
%
% Computes kinetic energy solving the integral over the domain throught
% Gaussian Quadrature
%
% Takes 1 input argument:
%   
%   data - Field structure, as read by smpm_read_fieldfile, containing
%           the field data (i.e. grid variables, field variables, time, 
%           etc...) for the run of interest.
%
% Returns 1 value:
%   
%   kE     - Scalar value resulting from Kinetic Energy computation
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

% 2) Compute the integral - Gauss quadrature
ff=ux.^2+uy.^2+uz.^2;        %Discrete function to be integrated
% 2.1) Integral in the x direction
[~, wi, ~] = lglnodes( double( n ) - 1 ); %GLL weight and points
a=x(1);                      %Limits of integration
b=x(end);                    %Limits of integration
w=((b-a)/2)*wi;                %Mapping weights
int_x=ff*repmat(w,mz,1);     %Compute the integral in x direction
% 2.2) Integral in the z direction
a=z(1);                      %Limits of integration
b=z(end);                    %Limits of integration
w=((b-a)/2)*wi;                %Mapping weights
int_z=int_x'*repmat(w,mz,1); %Compute the integral in z direction


KE2 = integrate_grid_function_2D( ff, x, z, n, mx, mz-1 );
kE=0.5*int_z;

end

