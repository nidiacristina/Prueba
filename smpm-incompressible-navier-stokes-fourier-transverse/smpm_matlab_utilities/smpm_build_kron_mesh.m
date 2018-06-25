function [x, z] = smpm_build_kron_mesh( n, xlims, zlims )
% [x, z] = smpm_build_kron_mesh( n, xlims, zlims )
%
% Generate a mesh made by taking zlims \kron xlims, in xi-first indexing, for
% reading by the SMPM solver code.
%
% Takes 3 arguments:
%
%   n     - Number of GLL points per direction, per subdomain.
%   xlims - Row vector, of length 1 x (mx + 1), specifying the bounds
%           for the x-dimension.
%   zlims - Row vector, of length 1 x (mz + 1), specifying the bounds
%           for the z-dimension.
%
% Returns 2 values:
%
%   x - Column vector, of length n^2 * mx * mz x 1, containing the x
%       coordinates of the generated mesh.  Coordinates are indexed xi-first.
%   z - Column vector, of length n^2 * mx * mz x 1, containing the z
%       coordinates of the generated mesh.  Coordinates are indexed xi-first.
%
% 7 May 2014
% Sumedh Joshi

% Get some constants.
mx = length( xlims ) - 1;
mz = length( zlims ) - 1;

% Get the GLL points.
xi  = lglnodes( n - 1 );
xi  = sort( xi );
eta = xi;

% Build the 1D multidomain mesh in xi.
XI = zeros( mz * n, 1 );
for ii = 1:mz

    % Grab the subdomain boundaries.
    iiz0 = zlims(ii);
    iiz1 = zlims(ii + 1);
    iihz = iiz1 - iiz0;

    % Build the local grid.
    iiz = (eta + 1) * iihz / 2 + iiz0;

    % Put the local grid in the right place.
    iistart           = (ii - 1) * n + 1;
    iiend             = (ii - 1) * n + n;
    XI(iistart:iiend) = iiz;

end

% Build the 1D multidomain mesh in eta.
ETA = zeros( mz * n, 1 );
for ii = 1:mx

    % Grab the subdomain boundaries.
    iix0 = xlims(ii);
    iix1 = xlims(ii + 1);
    iihx = iix1 - iix0;

    % Build the local grid.
    iix = (eta + 1) * iihx / 2 + iix0;

    % Put the local grid in the right place.
    iistart            = (ii - 1) * n + 1;
    iiend              = (ii - 1) * n + n;
    ETA(iistart:iiend) = iix;

end

% Do the meshgrid.
[x, z] = meshgrid( ETA, XI );
x      = x(:);
z      = z(:);

end
