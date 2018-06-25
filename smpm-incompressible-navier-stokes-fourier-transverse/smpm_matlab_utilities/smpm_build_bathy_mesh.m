function [x, z] = smpm_build_bathy_mesh( n, mx, mz, xbathy, zbathy )
% [x, z] = smpm_build_bathy_mesh( n, mx, mz, xbathy, zbathy )
%
% Generate a bathymetry mesh in xi-first indexing for reading by the SMPM
% solver code.  Note that the ocean surface will be place at z = 0, and so the
% bathymetry should have all z < 0.
%
% Takes 5 arguments:
%
%   n      - Number of GLL points per direction, per subdomain.
%   mx     - Number of subdomains in the x-direction.
%   mz     - Number of subdomains in the z-direction.
%   xbathy - Row vector, of length 1 x r, containing x-coordinates of
%            bathymetry.
%   zbathy - Row vector, of length 1 x r, containing z-coordinates of
%            bathymetry.
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

% Build the 1D multidomain mesh in eta.
xlims   = [min( xbathy ), max( xbathy )];
Lx      = (xlims(2) - xlims(1)) / mx;
eta     = (eta + 1) * Lx / 2;
eta_bot = repmat( eta, mx, 1 );
offset  = repmat( [0:mx-1]' * Lx, 1, n )';
offset  = offset(:);
eta_bot = eta_bot + offset + xlims(1);

% Interpolate to get the grid on the bottom.
xi_bot = interp1( xbathy, zbathy, eta_bot );

% Get the grid on the top.
xi_top  = 0 * xi_bot;
eta_top = eta_bot;

% Loop over columns, filling in with a xi grid.
z = zeros( n * mz, n * mx );
for ii = 1:n*mx

    zlims    = [xi_bot(ii), xi_top(ii)];
    Lz       = (zlims(2) - zlims(1)) / mz;
    iixi     = xi;
    iixi     = (iixi + 1) * Lz / 2 + xi_bot(ii);
    iixi     = repmat( iixi, mz, 1 );
    offset   = repmat( [0:mz-1]' * Lz, 1, n )';
    offset   = offset(:);
    iixi     = iixi + offset;
    z(:, ii) = iixi;

end

% Build x.
x = repmat( eta_bot', n * mz, 1 );
x = reshape( x, n * n * mz * mx, 1 );

% Build z.
z = z(:);

end
