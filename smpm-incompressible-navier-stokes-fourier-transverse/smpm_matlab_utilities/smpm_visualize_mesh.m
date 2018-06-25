function h = smpm_visualize_mesh( varargin )
% h = smpm_visualize_mesh( mesh_file_name[, h] )
% h = smpm_visualize_mesh( n, mx, mz, x, z[, h] )
% h = smpm_visualize_mesh( data.grid[, h] );
%
% Generate a visualization of a meshfile that plots all subdomain and external
% interfaces but does not plot the GLL grid.
%
% Takes 2 arguments:
%
%   mesh_file_name - Name of the init file that contains the desired mesh.
%   h              - Optional figure handle to plot onto.
%
% -OR-
%
% Takes 6 arguments:
%
%   n  - GLL points per direction per element.
%   mx - number of x subdomains.
%   mz - number of z subdomains.
%   x  - x-coordinates as read in by smpm_read_initfile.
%   z  - z-coordinates as read in by smpm_read_initfile.
%   h  - figure handle to plot into (optional).
%
% Returns 1 value:
%
%   h - Figure handle to the figure generated.
%
% 23 Jan 2015
% Sumedh Joshi

% Parse arguments.
switch nargin
   case 1
      if ~isstruct( varargin{1} )
         mesh_file_name = varargin{1};
         [data]         = smpm_read_initfile( mesh_file_name );
         n = data.grid.n; mx = data.grid.mx; mz = data.grid.mz; x = data.grid.x; z = data.grid.z;
      else
         grid = varargin{1};
         n = grid.n; mx = grid.mx; mz = grid.mz; x = grid.x; z = grid.z;
      end
      h = figure;
   case 2
      mesh_file_name = varargin{1};
      h              = varargin{2};
      if ~isstruct( varargin{1} )
         [n, mx, mz, x, z] = smpm_read_initfile( mesh_file_name );
         n = data.grid.n; mx = data.grid.mx; mz = data.grid.mz; x = data.grid.x; z = data.grid.z;
      else
         grid = varargin{1};
         n = grid.n; mx = grid.mx; mz = grid.mz; x = grid.x; z = grid.z;
      end
   case 5
      n  = varargin{1};
      mx = varargin{2};
      mz = varargin{3};
      x  = varargin{4};
      z  = varargin{5};
      h  = figure;
   case 6
      n  = varargin{1};
      mx = varargin{2};
      mz = varargin{3};
      x  = varargin{4};
      z  = varargin{5};
      h  = varargin{6};
end

% Setup the figure.
figure( h );
hold on;
whiteness = 0.0;

% Loop over the logically horizontal interfaces.
for ii = 1:mz+1

    % Get the grid points.
    a   = (ii - 1) * n + 1;
    iix = x(a:n*mz:end);
    iiz = z(a:n*mz:end);

    % Plot this horizontal line.
    plot( iix, iiz, 'color', whiteness * [1 1 1] );

end

% Loop over the logically vertical interfaces.
for ii = 1:mx

    % Get the grid points.
    iistart = (ii - 1) * n * n * mz + 1;
    iiend   = (ii - 1) * n * n * mz + n * mz;
    iix     = x(iistart:iiend);
    iiz     = z(iistart:iiend);

    % Plot this vertical line.
    plot( iix, iiz, 'color', whiteness * [1 1 1] );

end

% Plot the right-most and top-most outer boundaries.
r = n * n * mx * mz;
s = n * mz;

% Right-most.
iix = x(r - s + 1:r);
iiz = z(r - s + 1:r);
plot( iix, iiz, 'color', whiteness * [1 1 1] );

% Top-most.
iix = x(s:s:r);
iiz = z(s:s:r);
plot( iix, iiz, 'color', whiteness * [1 1 1] );

end
