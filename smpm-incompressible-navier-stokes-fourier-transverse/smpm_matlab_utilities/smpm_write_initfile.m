function smpm_write_initfile( n, mx, my, mz, x, y, z, rho0, ux0, uy0, uz0, ubc0, dubcdz0, rho0_bar, rho0_bar_z, init_file_name )
% smpm_write_initfile( n, mx, my, mz, x, y, z, rho0, ux0, uy0, uz0, ubc0, dubcdz0, rho0_bar, rho0_bar_z, init_file_name )
%
% Writes an initial conditions file containing the mesh and conditions
% specified 3D multidomain.  Each condition may be specified as a constant
% value across the grid, or on a per grid point basis (in xi-first indexing).
%
% Takes 15 arguments:
%
%   n              - Number of GLL points per direction, per subdomain.
%   mx             - Number of subdomains in the x-direction.
%   my             - Number of transverse grid points.
%   mz             - Number of subdomains in the z-direction.
%   x              - Vector, of length n^2 * mx * mz * my, containing the
%                    x-coordinates of the mesh associated with the field
%                    variables.
%   y              - Vector, of length n^2 * mx * mz * my, containing the
%                    y-coordinates of the mesh.
%   z              - Vector, of length n^2 * mx * mz * my, containing the
%                    z-coordinates of the mesh associated with the field
%                    variables.
%   rho0           - Scalar, or vector, specifying the initial density of
%                    the field.  Dimension n^2 * mx * mz * my.
%   ux0            - Scalar, or vector, specifying the initial velocity in
%                    the x-direction of the field. Dimension n^2 * mx * mz * my.
%   uy0            - Scalar, or vector, specifying the initial velocity in
%                    the y-direction of the field. Dimension n^2 * mx * mz * my.
%   uz0            - Scalar, or vector, specifying the initial velocity in
%                    the z-direction of the field. Dimension n^2 * mx * mz * my.
%                    the x-direction of the field.
%   ubc            - Scalar, or vector, specifying the background velocity in
%                    the x-direction of the field.  Dimension n^2 * mx * mz * my.
%   dubcdz         - Scalar, or vector, specifying the background velocity
%                    derivative in the x-direction of the field.  Dimension
%                    n^2 * mx * mz * my.
%   rho0_bar       - Scalar, or vector, specifying the initial background
%                    density of the field. Dimension n^2 * mx * mz * my.
%   rho0_bar_z     - Scalar, or vector, specifying the derivative of the initial 
%                    background density of the field. Dimension n^2 * mx * mz * my.
%   init_file_name - String indicating the file to save the generated initial
%                    conditions in.
%
% Returns nothing.
%

% 7 Mar 2013
% Sumedh Joshi

    % remove an existing initial conditions file.  we start fresh, rather than
    % update, so as to ensure that the dimensions of each parameter is
    % correct.  the HDF5 interface will happily write a subset of data in an
    % existing dataset which is very much not what we want.
    if exist( init_file_name, 'file' )
        delete( init_file_name );
    end

    % Create and write all the constants.
    h5create( init_file_name, '/grid/n',  1, 'Datatype', 'int32' );
    h5create( init_file_name, '/grid/mx', 1, 'Datatype', 'int32' );
    h5create( init_file_name, '/grid/my', 1, 'Datatype', 'int32' );
    h5create( init_file_name, '/grid/mz', 1, 'Datatype', 'int32' );

    h5write( init_file_name, '/grid/n', n );
    h5write( init_file_name, '/grid/mx', mx );
    h5write( init_file_name, '/grid/my', my );hdf5read
    h5write( init_file_name, '/grid/mz', mz );

    % Write the grid.
    h5create( init_file_name, '/grid/x', double( [mz*n, mx*n, my] ) );
    h5create( init_file_name, '/grid/y', double( [mz*n, mx*n, my] ) );
    h5create( init_file_name, '/grid/z', double( [mz*n, mx*n, my] ) );

    h5write( init_file_name, '/grid/x', reshape( x, [mz*n, mx*n, my] ) );
    h5write( init_file_name, '/grid/y', reshape( y, [mz*n, mx*n, my] ) );
    h5write( init_file_name, '/grid/z', reshape( z, [mz*n, mx*n, my] ) );

    % If any of the passed in parameters are scalars, expand them as a constant vector onto the grid.
    if numel( ux0 ) == 1
      ux0 = ux0 * ones( mz*n, mx*n, my );
    end
    if numel( uy0 ) == 1
      uy0 = uy0 * ones( mz*n, mx*n, my );
    end
    if numel( uz0 ) == 1
      uz0 = uz0 * ones( mz*n, mx*n, my );
    end
    if numel( ubc0 ) == 1
      ubc0 = ubc0 * ones( mz*n, mx*n, my );
    end
    if numel( dubcdz0 ) == 1
      dubcdz0 = dubcdz0 * ones( mz*n, mx*n, my );
    end
    if numel( rho0 ) == 1
      rho0 = rho0 * ones( mz*n, mx*n, my );
    end
    if numel( rho0_bar ) == 1
      rho0_bar = rho0_bar * ones( mz*n, mx*n, my );
    end
    if numel( rho0_bar_z ) == 1
      rho0_bar_z = rho0_bar_z * ones( mz*n, mx*n, my );
    end

    % Write the initial conditions.
    h5create( init_file_name, '/ic/rho',       double( [ mz*n, mx*n, my ] ) );
    h5create( init_file_name, '/ic/ux',        double( [ mz*n, mx*n, my ] ) );
    h5create( init_file_name, '/ic/uy',        double( [ mz*n, mx*n, my ] ) );
    h5create( init_file_name, '/ic/uz',        double( [ mz*n, mx*n, my ] ) );
    h5create( init_file_name, '/ic/ubc',       double( [ mz*n, mx*n, my ] ) );
    h5create( init_file_name, '/ic/dubcdz',    double( [ mz*n, mx*n, my ] ) );
    h5create( init_file_name, '/ic/rho_bar',   double( [ mz*n, mx*n, my ] ) );
    h5create( init_file_name, '/ic/rho_bar_z', double( [ mz*n, mx*n, my ] ) );

    h5write( init_file_name, '/ic/rho',       reshape( rho0, mz*n, mx*n, my ) );
    h5write( init_file_name, '/ic/ux',        reshape( ux0,  mz*n, mx*n, my ) );
    h5write( init_file_name, '/ic/uy',        reshape( uy0,  mz*n, mx*n, my ) );
    h5write( init_file_name, '/ic/uz',        reshape( uz0,  mz*n, mx*n, my ) );
    h5write( init_file_name, '/ic/ubc',       reshape( ubc0, mz*n, mx*n, my ) );
    h5write( init_file_name, '/ic/dubcdz',    reshape( dubcdz0, mz*n, mx*n, my ) );
    h5write( init_file_name, '/ic/rho_bar',   reshape( rho0_bar, mz*n, mx*n, my ) );
    h5write( init_file_name, '/ic/rho_bar_z', reshape( rho0_bar_z, mz*n, mx*n, my ) );

end
