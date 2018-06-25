function raycoeff = smpm_build_sponge_layer( xm, dxL, dxR, c0, dt )
% raycoeff = smpm_build_bathy_sponge_layer( xm, dxL, dxR, c0, dt )
%
% Generate a sponge layer in the streamwise direction based on the methodology
% found in Ammar's Thesis. The procedure can also be found in forcing.f of
% Peter Diamessis's 2D SMPM code as programmed by Tak Sakai.
%
% Takes 5 arguments:
%
%   x   - Column vector, of length n^2 * mx * mz x 1, containing the x
%         coordinates of the generated mesh.  Coordinates are indexed xi-first.
%   dxL - Thickness of sponge layer at the left.
%   dxR - Thickness of sponge layer at the right.
%   c0  - Scalar wave celerity.
%   dt  - Scalar timestep for the run.
%
% Returns 1 value:
%
%   raycoeff - Column vector, of length n^2 * mx * mz x 1, containing the
%              value of the Rayleigh Coefficient of the sponge layer
%              throughout the domain.
%

% 8 June 2016
% Gustavo Rivera

% Allocate Rayleight coefficient array.
raycoeff = zeros( size( xm(:) ) );

% Compute length of the domain.
Lx = max( xm ) - min( xm );

% Set tolerance.
eps = 1.0e-6;

% Set time scale of layer thickness.
kappax = (2. * c0) / (dxL + dxR) * (-log( eps ));

% Set the beginning coordinate.
x0 = xm(1);

% Set the end of layer.
xend = x0 + Lx;

% Build the Rayleight coefficients throughout the domain.
for i = 1:length( xm )

    xx = xm(i);

    % Construct layer.
    if xx < x0 + dxL
        eta0        = x0 + dxL;
        eta         = pi./2 .* (xx - eta0) / dxL;
        raycoeff(i) = kappax .* sin( eta ).^2;
    elseif xx > xend - dxR
        eta0        = xend - dxR;
        eta         = pi./2 .* (xx - eta0) / dxR;
        raycoeff(i) = kappax .* sin( eta ).^2;
    else
        raycoeff(i) = 0.0;
    end

    if raycoeff(i) > 0.5 / dt
        raycoeff(i) = 0.5 / dt;
    end

end

end
