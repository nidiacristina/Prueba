function status = smpm_compare_fieldfiles( field_file_name1, field_file_name2, tolerance, verbosity_flag, display_flag )
% status = smpm_compare_fieldfiles( field_file_name1, field_file_name2, tolerance, verbosity_flag, display_flag )
%
% Reads the contents of two field files and determines if they are the same to
% within the relative tolerance provided.  Optionally plots the differences as
% ratios to identify structure in the differing field variables.
%
% Takes 5 arguments:
%
%   field_file_name1 - String indicating the first field file to compare.
%   field_file_name2 - String indicating the second field file to compare.
%   tolerance        - Optional, absolute tolerance used in determining
%                      equality between field variables in field_file_name1
%                      and field_file_name2.  If omitted, defaults to 1e-14.
%   verbosity_flag   - Optional flag indicating whether the function should
%                      display comparison results to standard out.  If omitted,
%                      defaults to true so that the comparison operation is
%                      verbose.
%   display_flag     - Optional flag indicating whether the function should
%                      display figures of the field variables at each timestep.
%                      If omitted, defaults to false so that no figures are
%                      generated.
%
% Returns 1 value:
%
%   status - Integer indicating the result of the comparison.  If the two field
%            files' contents are identical within tolerance, status is 0.  If
%            the two field files have different numbers of timesteps, though the
%            common contents are identical within tolerance, status is 1.  Otherwise,
%            status is 2.

if nargin < 5
    display_flag = [];
end

if nargin < 4
    verbosity_flag = [];
end

if nargin < 3
    tolerance = [];
end

if isempty( tolerance )
    tolerance = 1e-14;
end

if isempty( verbosity_flag )
    verbosity_flag = 1;
end

if isempty( display_flag )
    display_flag = 0;
end

% by default, assume that both field files contain the same thing.
status = 0;

data = smpm_read_fieldfile( field_file_name1 );
x1   = data.grid.x;    y1 = data.grid.y;    z1 = data.grid.z;
ux1  = data.field.ux; uy1 = data.field.uy; uz1 = data.field.uz; rho1 = data.field.rho;
data = smpm_read_fieldfile( field_file_name2 );
x2   = data.grid.x;    y2 = data.grid.y;    z2 = data.grid.z;
ux2  = data.field.ux; uy2 = data.field.uy; uz2 = data.field.uz; rho2 = data.field.rho;

number_timesteps1 = size( ux1, 4 );
number_timesteps2 = size( ux2, 4 );

% ensure that the field files are compatible for comparison.
if size( x1 ) ~= size( x2 )
    error( sprintf( 'X coordinates from ''%s'' have a different size than those from ''%s''.', ...
                    field_file_name1, field_file_name2 ) );
elseif size( y1 ) ~= size( y2 )
    error( sprintf( 'Y coordinates from ''%s'' have a different size than those from ''%s''.', ...
                    field_file_name1, field_file_name2 ) );
elseif size( z1 ) ~= size( z2 )
    error( sprintf( 'Z coordinates from ''%s'' have a different size than those from ''%s''.', ...
                    field_file_name1, field_file_name2 ) );
end

% Let the caller know things cannot be identical if one file can only be a
% potential subset of the other.
if number_timesteps1 ~= number_timesteps2
    warning( sprintf( '''%s'' has %d timesteps, while ''%s'' has %d.  Comparing the first %d.', ...
                      field_file_name1, number_timesteps1, field_file_name2, number_timesteps2, ...
                      min( number_timesteps1, number_timesteps2 ) ) );
end

% Walk through each of the common timesteps and note any differences.
%
% TODO: Indicate where in the grid the differences are.  Or at least one/the
%       first of them.

for timestep_index = 1:min( number_timesteps1, number_timesteps2 )
    % Keep a running count of the field variables that differ in this
    % timestep.  We use this to determine if we need to plot anything.
    local_differences = 0;

    % NOTE: The conditionals here are independent so that we can see exactly what is
    %       different within each timestep.
    if any( abs( (ux1(:, :, :, timestep_index) - ux2(:, :, :, timestep_index)) ./ ux1(:, :, :, timestep_index) ) >= tolerance )
        status            = 2;
        local_differences = local_differences + 1;

        if verbosity_flag
            warning( sprintf( 'Timestep #%d''s ux differs.', timestep_index ) );
        end
    end
    if any( abs( (uy1(:, :, :, timestep_index) - uy2(:, :, :, timestep_index)) ./ uy1(:, :, :, timestep_index) ) >= tolerance )
        status            = 2;
        local_differences = local_differences + 1;

        if verbosity_flag
            warning( sprintf( 'Timestep #%d''s uy differs.', timestep_index ) );
        end
    end
    if any( abs( (uz1(:, :, :, timestep_index) - uz2(:, :, :, timestep_index)) ./ uz1(:, :, :, timestep_index) ) >= tolerance )
        status            = 2;
        local_differences = local_differences + 1;

        if verbosity_flag
            warning( sprintf( 'Timestep #%d''s uz differs.', timestep_index ) );
        end
    end
    if any( abs( (rho1(:, :, :, timestep_index) - rho2(:, :, :, timestep_index)) ./ rho1(:, :, :, timestep_index) ) >= tolerance )
        status            = 2;
        local_differences = local_differences + 1;

        if verbosity_flag
            warning( sprintf( 'Timestep #%d''s rho differs.', timestep_index ) );
        end
    end

    % If requested, display a few plots per timestep showing the field
    % variables and their differences.
    if (local_differences > 0) && display_flag
        % Compute basic statistics about the differences.
        x_ratio   =  ux1(:, :, :, timestep_index) ./  ux2(:, :, :, timestep_index);
        y_ratio   =  uy1(:, :, :, timestep_index) ./  uy2(:, :, :, timestep_index);
        z_ratio   =  uz1(:, :, :, timestep_index) ./  uz2(:, :, :, timestep_index);
        rho_ratio = rho1(:, :, :, timestep_index) ./ rho2(:, :, :, timestep_index);

        % Deal with data that are zero (0/0 -> NaN) so we avoid warnings when
        % we plot.
        %
        % NOTE: This does not properly handle NaN's in solver output.
        x_ratio(isnan( x_ratio ))     = 0;
        y_ratio(isnan( y_ratio ))     = 0;
        z_ratio(isnan( z_ratio ))     = 0;
        rho_ratio(isnan( rho_ratio )) = 0;

        std_x   = max( std(   x_ratio(:) ), eps( 'single' ) );
        std_y   = min( std(   y_ratio(:) ), eps( 'single' ) );
        std_z   = min( std(   z_ratio(:) ), eps( 'single' ) );
        std_rho = min( std( rho_ratio(:) ), eps( 'single' ) );

        for y_index = 1:size( ux1, 3 )
            % Plot a transverse plane for each field variable.  We stack them
            % in rows since the SMPM solver targets high aspect ratio grids.
            figure( 'name', sprintf( 'Ratios for timestep #%d (Y=%d)', ...
                                     timestep_index, y_index ) );
            subplot( 4, 1, 1 ); imagesc(  x_ratio(:, :, y_index) ); axis equal;
            caxis( [-std_x, std_x] + 1 ); colorbar; axis xy; ylabel( 'Z' ); title( 'ux' ); ylabel( 'Z' );
            subplot( 4, 1, 2 ); imagesc(  y_ratio(:, :, y_index) ); axis equal;
            caxis( [-std_y, std_y] + 1); colorbar; axis xy; ylabel( 'Z' ); title( 'uy' )
            subplot( 4, 1, 3 ); imagesc(  z_ratio(:, :, y_index) ); axis equal;
            caxis( [-std_z, std_z] + 1 ); colorbar; axis xy; ylabel( 'Z' ); title( 'uz' )
            subplot( 4, 1, 4 ); imagesc( rho_ratio(:, :, y_index) ); axis equal;
            caxis( [-std_rho, std_rho] + 1 ); colorbar; axis xy; xlabel( 'X' ); ylabel( 'Z' ); title( 'rho' )
        end
    end
end

% handle the case where one field file was shorter than the other, but
% they contained the same data in the common portion.
if status == 0 && number_timesteps1 ~= number_timesteps2
    status = 1;
end

end
