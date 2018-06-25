function [x, y, z] = smpm_extrude_mesh( n, mx, mz, x, y, z )
% [x, y, z] = smpm_extrude_mesh( n, mx, mz, x, y, z )
%
% Generate a 3D mesh in z/x/y-first indexing for reading by the SMPM solver
% code.
%
% NOTE: The zlayer matrix should have zlayer(1, :) be the bottom-most layer.
%
% Takes 6 arguments:
%
%   n  - Number of collocation points.
%   mx - Number of domains in the X dimension.
%   mz - Number of domains in the Z dimension.
%   x  - 2D array, as built by SMPM meshing tools, representing the
%        y-coordinate.
%   y  - 1D array representing the y-coordinate.
%   z  - 2D array, as built by SMPM meshing tools, representing the
%        y-coordinate.
%
% Returns 3 values:
%
%   x - Matrix, of dim [n * mz, n * mx, my], containing the x
%       coordinates of the generated mesh.  Coordinates are indexed xi-first.
%   y - Matrix, of dim [n * mz, n * mx, my], containing the y
%       coordinates of the generated mesh.  Coordinates are indexed xi-first.
%   z - Matrix, of dim [n * mz, n * mx, my], containing the z
%       coordinates of the generated mesh.  Coordinates are indexed xi-first.
%

% 16 Jan 2016
% Sumedh Joshi

% Reshape the grid.
my = length( y );
x  = reshape( x, [mz * n, mx * n, 1] );
z  = reshape( z, [mz * n, mx * n, 1] );

% Extrude the 2D grid.
x  = repmat( x, [1, 1, my] );
z  = repmat( z, [1, 1, my] );
y  = reshape( y, [1, 1, my] );
y  = repmat( y, [mz * n, mx * n, 1] );
