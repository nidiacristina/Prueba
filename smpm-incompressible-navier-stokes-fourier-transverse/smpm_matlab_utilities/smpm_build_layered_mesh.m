function [x, z] = smpm_build_layered_mesh( n, mx, mz, xlayer, zlayer )
% [x, z] = smpm_build_layered_mesh( n, mx, mz, xlayer, zlayer, mesh_file_name )
%
% Generate a mesh in xi-first indexing for reading by the SMPM solver code.
%
% This mesh creator asks the user to specify all of the horizontal layers
% between subdomains in the global mesh, and interpolates those based on
% the mesh resolution parameters given.
%
% NOTE: The zlayer matrix should have zlayer(1, :) be the bottom-most layer.
%
% Takes 5 arguments:
%
%   n      - Number of GLL points per direction, per subdomain.
%   mx     - Number of subdomains in the x-direction.
%   mz     - Number of subdomains in the z-direction.
%   xlayer - Row vector, of length 1 x r, specifying the x coordinates of all
%            the layers.
%   zlayer - Matrix, of size (mz + 1) x r, specifying the z coordinates of all
%            the layers, for each x coordinate.
%
% Returns 2 values:
%
%   x - Column vector, of length n^2 * mx * mz x 1, containing the x
%       coordinates of the generated mesh.  Coordinates are indexed xi-first.
%   z - Column vector, of length n^2 * mx * mz x 1, containing the z
%       coordinates of the generated mesh.  Coordinates are indexed xi-first.
%

% 25 Mar 2013
% Sumedh Joshi

% Do some error checking.
if size( zlayer, 1 )  ~= mz + 1
    error( 'Number of z layers does not equal mz + 1' );
end

% Get the GLL points.
xi  = lglnodes( n - 1 );
xi  = sort( xi );
eta = xi;

% Build the 1D multidomain mesh in eta.
xlims   = [min( xlayer ), max( xlayer )];
Lx      = (xlims(2) - xlims(1)) / mx;
eta     = (eta + 1) * Lx / 2;
eta_bot = repmat( eta, mx, 1 );
offset  = repmat( [0:mx-1]' * Lx, 1, n )';
offset  = offset(:);
x_int   = eta_bot + offset + xlayer(1);

% Interpolate to get the grid at each layer.
z_int = zeros( mz + 1, n * mx );
for ii = 1:mz + 1
    z_int(ii, :) = interp1( xlayer, zlayer(ii, :), x_int );
end

% Loop over the x-grid, filling in vertical layers.
z = zeros( n * n * mx * mz, 1 );
x = zeros( n * n * mx * mz, 1 );
for ii = 1:length( x_int )

    % Loop over this vertical grid, filling in between the layers.
    for jj = 1:mz

        % Get the start and stop indices.
        jjstart = (ii - 1) * n * mz + (jj - 1) * n + 1;
        jjend   = (ii - 1) * n * mz + (jj - 1) * n + n;

        % Grab the limits in the vertical in this column.
        zlims = [z_int(jj, ii), z_int(jj+1, ii)];
        jjLz  = zlims(2) - zlims(1);

        % Stretch the xi coordinate in this column to match the
        % boundaries.
        jjz = ((xi + 1) / 2) * jjLz + zlims(1);

        % Store.
        x(jjstart:jjend) = x_int(ii);
        z(jjstart:jjend) = jjz;

    end
end

end
