function [APE, KE] = smpm_compute_wave_energy( Ea, Ek, z, x, n, mx, mz )
% [APE, KE] = smpm_compute_wave_energy( Ea, Ek, z, x, n, mx, mz )
%
% Computes the available potential energy (APE) and kinetic energy (KE)
% by integrating the energy density per unit Volume Across all domain
% using the GLL Quadrature Integral from Sumedh Joshi.
%
%   Ea is taken to be rho_w*g*eta where eta is some reference length.
%   Ek is taken to be rho_0/2*(u^2 + w^2) where u and w are the wave
%   induced velocity fields.
%
% For more information please refer to Lamb 2002 and Lamb & Nguyen 2009.
%
% Takes 7 arguments:
%
%   Ek  - Column vector, of length n^2 * mx * mz x 1, containing the
%         potential energy density per unit volume.
%   Ek  - Column vector, of length n^2 * mx * mz x 1, containing the
%         kinetic energy density per unit volume.
%   x   - Column vector, of length n^2 * mx * mz x 1, containing the x-
%         coordinates of the generated mesh.
%   z   - Column vector, of length n^2 * mx * mz x 1, containing the z-
%         coordinates of the generated mesh.
%   n   - Number of GLL points per direction, per subdomain.
%   mx  - Number of subdomains in the x-direction.
%   mz  - Number of subdomains in the z-direction.
%
% Returns 2 values:
%
%   APE  - Available potential energy [J/m] in the domain. (XXX: size)
%   KE   - Kinetic energy [J/m] in the domain. (XXX: size)
%

% 4 March 2016
% Gustavo Rivera

% Compute Available Potential Energy
APE = integrate_grid_function_2D( Ea, x, z, n, mx, mz );

% Compute Kinetic Energy
KE  = integrate_grid_function_2D( Ek, x, z, n, mx, mz );

end
