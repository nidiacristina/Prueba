function [x, z] = smpm_build_cartesian_mesh( n, mx, mz, xlims, zlims )
% [x, z] = smpm_build_cartesian_mesh( n, mx, mz, xlims, zlims )
%
% Generate a cartesian (undeformed) mesh in xi-first indexing for reading by
% the SMPM solver code.
%
% Takes 5 arguments:
%
%   n     - Number of GLL points per direction, per subdomain.
%   mx    - Number of subdomains in the x-direction.
%   mz    - Number of subdomains in the z-direction.
%   xlims - Row vector, of length 1 x 2, specifying bounds for the
%           x-dimension.
%   zlims - Row vector, of length 1 x 2, specifying bounds for the
%           z-dimension.
%
% Returns 2 values:
%
%   x - Column vector, of length n^2 * mx * mz x 1, containing the x
%       coordinates of the generated mesh.  Coordinates are indexed xi-first.
%   z - Column vector, of length n^2 * mx * mz x 1, containing the z
%       coordinates of the generated mesh.  Coordinates are indexed xi-first.
%
% 7 Mar 2013
% Sumedh Joshi

% Get the GLL points.
xi  = lglnodes( n - 1 );
xi  = sort( xi );
eta = xi;

% Build the 1D multidomain mesh in xi.
Lz     = (zlims(2) - zlims(1)) / mz;
xi     = (xi + 1) * Lz / 2 - (zlims(2) - zlims(1)) / 2;
XI     = repmat( xi, mz, 1 );
offset = repmat( [0:mz-1]' * Lz, 1, n )';
offset = offset(:);
XI     = XI + offset;

% Build the 1D multidomain mesh in eta.
Lx     = (xlims(2) - xlims(1)) / mx;
eta    = (eta + 1) * Lx / 2 - (xlims(2) - xlims(1)) / 2;
ETA    = repmat( eta, mx, 1 );
offset = repmat( [0:mx-1]' * Lx, 1, n )';
offset = offset(:);
ETA    = ETA + offset;

% Make the first entry of both grid vectors zero.
ETA = ETA - ETA(1);
XI  = XI  - XI(1);

% Do the meshgrid.
[x z] = meshgrid( ETA, XI );
x     = x(:);
z     = z(:);

% Map to the desired starting location.
x = x + xlims(1);
z = z + zlims(1);

%    % Add back in the offsets.
%    x = x + (xlims(2) - xlims(1)) / 2;
%    z = z + (zlims(2) - zlims(1)) / 2;
end
