function smpm_create_file_from_checkpoint( data , ic_file_name )
% smpm_create_file_from_checkpoint( data, ic_file_name )
%
% Writes a file containing the restart conditions obtained from the 
% checkpoint of interest via smpm_read_restart.m.  
%
% Takes 2 arguments:
%   data              - Struct of variable name/values to write (c.f. as
%                       returned by smpm_read_checkpoint_from_restart).
%   ic_file_name      - String indicating the file to save the generated 
%                       conditions in.
%
% Returns nothing.
%
%
% 27 May 2017
% Gustavo Rivera w/ Greg Thomsen & Sumedh Joshi

    % remove an existing ic conditions file, assuming its named the same
    % as ic_file_name. We start fresh, rather than
    % update, so as to ensure that the dimensions of each parameter is
    % correct.  the HDF5 interface will happily write a subset of data in an
    % existing dataset which is very much not what we want.
    if exist( ic_file_name, 'file' )
        delete( ic_file_name );
    end

    % Create and write all the constants.
    h5create( ic_file_name, '/grid/n',  1, 'Datatype', 'int32' );
    h5create( ic_file_name, '/grid/mx', 1, 'Datatype', 'int32' );
    h5create( ic_file_name, '/grid/my', 1, 'Datatype', 'int32' );
    h5create( ic_file_name, '/grid/mz', 1, 'Datatype', 'int32' );

    h5write( ic_file_name, '/grid/n' , data.grid.n );
    h5write( ic_file_name, '/grid/mx', data.grid.mx );
    h5write( ic_file_name, '/grid/my', data.grid.my );
    h5write( ic_file_name, '/grid/mz', data.grid.mz );

    % Write the grid.
    h5create( ic_file_name, '/grid/x', double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
    h5create( ic_file_name, '/grid/y', double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
    h5create( ic_file_name, '/grid/z', double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );

    h5write( ic_file_name, '/grid/x', reshape( data.grid.x, [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ));
    h5write( ic_file_name, '/grid/y', reshape( data.grid.y, [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ));
    h5write( ic_file_name, '/grid/z', reshape( data.grid.z, [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ));
    
    % Write the timestep scalars
    h5create( ic_file_name, '/ic/dt',  1  );
    h5create( ic_file_name, '/ic/dt1', 1  );
    h5create( ic_file_name, '/ic/dt2', 1  );
    h5create( ic_file_name, '/ic/time', 1  );
    
    h5write( ic_file_name, '/ic/dt' , data.field.dt );
    h5write( ic_file_name, '/ic/dt1', data.field.dt1 );
    h5write( ic_file_name, '/ic/dt2', data.field.dt2 );
    h5write( ic_file_name, '/ic/time' , data.field.time );
   
    % If any of the passed in parameters are scalars, expand them as a constant vector onto the grid.

    % Write the ic conditions.   
    h5create( ic_file_name, '/ic/rho0',     double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] )  );
    h5create( ic_file_name, '/ic/ux0',      double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] )  );
    h5create( ic_file_name, '/ic/uy0',      double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] )  );
    h5create( ic_file_name, '/ic/uz0',      double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] )  );
    
    h5create( ic_file_name, '/ic/rho1',     double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] )  );
    h5create( ic_file_name, '/ic/ux1',      double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] )  );
    h5create( ic_file_name, '/ic/uy1',      double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] )  );
    h5create( ic_file_name, '/ic/uz1',      double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] )  );
    
    h5create( ic_file_name, '/ic/rho2',     double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] )  );
    h5create( ic_file_name, '/ic/ux2',      double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] )  );
    h5create( ic_file_name, '/ic/uy2',      double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] )  );
    h5create( ic_file_name, '/ic/uz2',      double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] )  );
    
    h5create( ic_file_name, '/constant_fields/ubc',       double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] )  );
    h5create( ic_file_name, '/constant_fields/dubcdz',    double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] )  );
    h5create( ic_file_name, '/constant_fields/rho_bar',   double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] )  );
    h5create( ic_file_name, '/constant_fields/rho_bar_z', double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] )  );
    h5create( ic_file_name, '/constant_fields/ux_b',      double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] )  );
    h5create( ic_file_name, '/constant_fields/uy_b',      double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] )  );
    h5create( ic_file_name, '/constant_fields/uz_b',      double( [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] )  );
   
    h5write( ic_file_name, '/ic/rho0',     reshape( data.field.rho0, [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
    h5write( ic_file_name, '/ic/ux0',      reshape( data.field.ux0,  [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
    h5write( ic_file_name, '/ic/uy0',      reshape( data.field.uy0,  [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
    h5write( ic_file_name, '/ic/uz0',      reshape( data.field.uz0,  [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
    
    h5write( ic_file_name, '/ic/rho1',     reshape( data.field.rho1, [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
    h5write( ic_file_name, '/ic/ux1',      reshape( data.field.ux1,  [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
    h5write( ic_file_name, '/ic/uy1',      reshape( data.field.uy1,  [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
    h5write( ic_file_name, '/ic/uz1',      reshape( data.field.uz1,  [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
    
    h5write( ic_file_name, '/ic/rho2',     reshape( data.field.rho2, [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
    h5write( ic_file_name, '/ic/ux2',      reshape( data.field.ux2,  [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
    h5write( ic_file_name, '/ic/uy2',      reshape( data.field.uy2,  [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
    h5write( ic_file_name, '/ic/uz2',      reshape( data.field.uz2,  [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
    
    h5write( ic_file_name, '/constant_fields/ubc',       reshape( data.field.ubc,      [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
    h5write( ic_file_name, '/constant_fields/dubcdz',    reshape( data.field.dubcdz,   [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
    h5write( ic_file_name, '/constant_fields/rho_bar',   reshape( data.field.rho_bar,  [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
    h5write( ic_file_name, '/constant_fields/rho_bar_z', reshape( data.field.rho_bar_z,[data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
    h5write( ic_file_name, '/constant_fields/ux_b',      reshape( data.field.ux_b,     [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
    h5write( ic_file_name, '/constant_fields/uy_b',      reshape( data.field.uy_b,     [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
    h5write( ic_file_name, '/constant_fields/uz_b',      reshape( data.field.uz_b,     [data.grid.mz*data.grid.n, data.grid.mx*data.grid.n, data.grid.my] ) );
end
