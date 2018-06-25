function data = smpm_read_initfile( init_file_name )
% data = smpm_read_initfile( init_file_name )
%
% Read an initial conditions file and returns its contents.
%
% Takes 1 argument:
%
%   init_file_name - String indicating the initial conditions file to read
%                    from disk.
%
% Returns 1 value:
%
%   data        - Structure with fields specified as below.
%     .grid     - Struct with fields specified as below.
%       n           - Number of GLL points per direction, per subdomain.
%       mx          - Number of subdomains in the x-direction.
%       my          - Number of transverse grid points.
%       mz          - Number of subdomains in the z-direction.
%       x           - Vector, of length n^2 * mx * mz x 1, containing the
%                     x-coordinates of the mesh associated with the field
%                      variables.
%       y           - Vector, of length my, containing the y-coordinates.
%       z           - Vector, of length n^2 * mx * mz x 1, containing the
%                 z-coordinates of the mesh associated with the field
%                 variables.
%     .ic       - Struct with fields specified as below.
%        rho0        - Matrix, of dim mz * n by mx * n by my containing the
%                      initial density at each grid point.
%        ux0         - Matrix, of dim mz * n by mx * n by my containing the
%                      initial x-velocity at each grid point.
%        uy0         - Matrix, of dim mz * n by mx * n by my containing the
%                      initial y-velocity at each grid point.
%        uz0         - Matrix, of dim mz * n by mx * n by my containing the
%                      initial z-velocity at each grid point.
%        rho0_bar    - Matrix, of dim mz * n by mx * n containing the
%                      initial background stratification at each grid point.
%

% 7 Mar 2013
% Sumedh Joshi

% Read the grid information and the grid.
data.grid.n  = h5read( init_file_name, '/grid/n' );
data.grid.mx = h5read( init_file_name, '/grid/mx' );
data.grid.my = h5read( init_file_name, '/grid/my' );
data.grid.mz = h5read( init_file_name, '/grid/mz' );
data.grid.x  = h5read( init_file_name, '/grid/x' );
data.grid.y  = h5read( init_file_name, '/grid/y' );
data.grid.z  = h5read( init_file_name, '/grid/z' );

% Read the initial conditions.
data.ic.ux        = h5read( init_file_name, '/ic/ux' );
data.ic.uy        = h5read( init_file_name, '/ic/uy' );
data.ic.uz        = h5read( init_file_name, '/ic/uz' );
data.ic.ubc       = h5read( init_file_name, '/ic/ubc' );
data.ic.dubcdz    = h5read( init_file_name, '/ic/dubcdz' );
data.ic.rho       = h5read( init_file_name, '/ic/rho' );
data.ic.rho_bar   = h5read( init_file_name, '/ic/rho_bar' );
data.ic.rho_bar_z = h5read( init_file_name, '/ic/rho_bar_z' );

end
