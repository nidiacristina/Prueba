function [differences, fig_h] = smpm_quantify_field_differences( reference_file_name, compared_file_names, quantify_uy, quantify_rho, tolerance_range, figure_prefix, difference_range )
% [differences, fig_h] = smpm_quantify_field_differences( reference_file_name, compared_file_names, quantify_uy, quantify_rho, tolerance_range, figure_prefix, difference_range )
%
% Quantifies the difference in one or more field outputs against a reference
% field output.  This is intended to enable quantification of differences
% between different versions/configurations of the solver to aide in
% development and debugging.  Control over which grid variables are compared
% is available allowing for comparison of legacy SMPM outputs as well as field
% outputs that don't generate every field variable.
%
% Differences are plotted in a figure for easy comparison and analysis.
%
% Takes 7 arguments:
%
%  reference_file_name - Path to the SMPM field file to use as a reference.
%  compared_file_names - Cell array of paths to SMPM field files to compare
%                        against reference_file_name.  May also be specified
%                        as a single string if only a single field file needs
%                        to be compared.
%  quantify_uy         - Optional flag indicating whether velocity in the
%                        transverse direction should be compared.  If omitted,
%                        defaults to false.
%  quantify_rho        - Optional flag indicating whether density should be
%                        compared.  If omitted, defaults to false.
%  tolerance_range     - Optional vector specifying the relative error
%                        threshold to use for quantification.  This is
%                        given X for a relative error of 10^X.  If omitted,
%                        defaults to 0:-1:-10.
%  figure_prefix       - Optional string specifying a tag to prefix the
%                        generated figure's title with.  If omitted, defaults
%                        to an empty string which does not alter the title
%                        at all.
%  difference_range    - Optional two element vector specifying the range
%                        of differences to zoom in on.  This governs the
%                        figure's Y-axis limits and defaults to [0, 100] if
%                        omitted.
%
% Returns 2 values:
%
%  differences - Returns a matrix, FxT where F is the number of files in
%                compared_file_names and T the number of thresholds in
%                tolerance_range, containing the percentage of grid points
%                that differed.
%  fig_h       - Handle to the figure produced.
%

if nargin < 7
    difference_range = [];
end
if nargin < 6
    figure_prefix = [];
end
if nargin < 5
    tolerance_range = [];
end
if nargin < 4
    quantify_rho = [];
end
if nargin < 3
    quantify_uy = [];
end

if isempty( quantify_uy )
    quantify_uy = 1;
end
if isempty( quantify_rho )
    quantify_rho = 1;
end
if isempty( tolerance_range )
    tolerance_range = 0:-1:-10;
end
if isempty( figure_prefix )
    figure_prefix = '';
end
if isempty( difference_range )
    difference_range = [0, 5];
end

% Ensure that we have a list of files to compare against our reference, even
% if that list only has a single entry.
if ~iscell( compared_file_names )
    compared_file_names = { compared_file_names };
end

% Ensure that our flags are logicals.
quantify_uy  = logical( quantify_uy );
quantify_rho = logical( quantify_rho );

% Ensure our tolerances are ordered most stringent to most permissive so we
% can properly label our axes.
tolerance_range = sort( tolerance_range );

if ~isempty( figure_prefix )
    figure_prefix = sprintf( '[%s] - ', figure_prefix );
end

% Check that we're on a system that can run our tools.
if ~isunix
    error( 'The current system is not Unix-based.  Odds we can''t run.' )
end

[status, command_path] = system( 'which h5diff' );
if status ~= 0
    error( 'h5diff is not in the path.' )
end

% Ensure that every file specified can be compared.
%
% NOTE: This only verifies that the grid configuration is identical rather
%       than validating the parameterization.  We assume the caller has
%       intelligently selected files to compare.
reference_info   = smpm_fieldfile_info( reference_file_name );
reference_fields = fieldnames( reference_info );
for file_index = 1:length( compared_file_names )
    current_file_name = compared_file_names{file_index};
    current_info      = smpm_fieldfile_info( current_file_name );

    if reference_info.mx ~= current_info.mx
        error( '%s''s number of X domains does not match %s''s (%d != %d).', ...
               reference_file_name, current_file_name, reference_info.mx ~= current_info.mx );
    elseif reference_info.mz ~= current_info.mz
        error( '%s''s number of Z domains does not match %s''s (%d != %d).', ...
               reference_file_name, current_file_name, reference_info.mz ~= current_info.mz );
    elseif reference_info.n ~= current_info.n
        error( '%s''s number of collocation points does not match %s''s (%d != %d).', ...
               reference_file_name, current_file_name, reference_info.n ~= current_info.n );
    elseif ~all( reference_info.x(:) == current_info.x(:) )
        error( '%s''s X coordinates do not match %s''s.', ...
               reference_file_name, current_file_name );
    elseif ~all( reference_info.z(:) == current_info.z(:) )
        error( '%s''s Z coordinates do not match %s''s.', ...
               reference_file_name, current_file_name );
    end

    % We only validate the Y dimensions if we're quantifying differences
    % there.
    if quantify_uy && reference_info.my ~= current_info.my
        error( '%s''s number of Y domains does not match %s''s (%d != %d).', ...
               reference_file_name, current_file_name, reference_info.my ~= current_info.my );
    elseif ~all( reference_info.y(:) == current_info.y(:) )
        error( '%s''s Y coordinates do not match %s''s.', ...
               reference_file_name, current_file_name );
    end
end

mx           = reference_info.mx;
my           = reference_info.my;
mz           = reference_info.mz;
n            = reference_info.n;
number_steps = reference_info.number_steps;

% Compute the expected number of differences based on the grid configuration
% and the parameters we're requested to quantify.
%
% NOTE: The grid is always 3D even if we're not considering velocities in the
%       Y-dimension.
%
% NOTE: We ignore the first timestep's field as it *should* be the same
%       between two simulations started from the same initial conditions.
%       This ensures that extremely tight tolerances have the proper
%       normalization factor to show 100% different.
%
% NOTE: We have to cast the value back to double to ensure our arithmetic
%       does not become polluted by integral data types.  These can be
%       returned by smpm_fieldfile_info() depending on the implementation
%       of its HDF5 reader (e.g. MATLAB vs Octave).
number_grid_points = n.^2 * mx * my * mz * (2 + quantify_uy + quantify_rho);
number_grid_points = double( number_grid_points * (number_steps - 1));

% Pre-allocate our differences matrix.  Indexed by file compared and a
% tolerance.
differences = zeros( length( compared_file_names ), length( tolerance_range ) );

% Track the differences for each file.
for file_index = 1:length( compared_file_names )
    current_file_name = compared_file_names{file_index};

    for tolerance_index = 1:length( tolerance_range )
        tolerance       = tolerance_range(tolerance_index);

        % Build a command to compare the field outputs for these two files and
        % reduce the output to just the numeric differences found.  This
        % assumes that h5diff generates output of the form:
        %
        %   dataset: </field/step8/ux> and </field/step8/ux>
        %   272 differences found
        %   ...
        %
        % We don't care where the differences are, just that there are some
        % for the current tolerance level.
        %
        % NOTE: Care is taken to ignore the initial conditions to ensure our
        %       percentages are correct.  These *should* be identical though
        %       will result in differences if they aren't.
        tolerance_str   = sprintf( '%.*f', ...
                                   max( 1, abs( tolerance ) ), ...
                                   10 .^ tolerance );
        compare_command = sprintf( 'sh -c "(h5diff --exclude-path /field/step0 -p %s %s %s /field; echo 0 differences) | grep differences | cut -d\\" \\" -f1"', ...
                                   tolerance_str, ...
                                   reference_file_name, ...
                                   current_file_name );

        % Run the comparison and make sure that we didn't run into any
        % problems that produced non-numeric output.
        [status, comparison_output] = system( compare_command );
        current_differences = str2num( comparison_output );
        if isempty( current_differences )
            warning( 'Comparing %s to %s with tolerance of 1e%d failed!', ...
                     reference_file_name, current_file_name, tolerance );
            current_differences = 0;
        end

        % Compute the simplest statistic we can and sum the differences found.
        differences(file_index, tolerance_index) = sum( current_differences );
    end
end

% Convert absolute differences into percentage differences.
differences = differences ./ number_grid_points .* 100.0;

% Generate multiple line plots, one per compared file, showing the percentage
% of total grid points (across all timesteps) that are different when compared
% against a relative tolerance.
fig_h = figure;
plot( tolerance_range, differences' )
title( sprintf( '%sDifferences relative to ''%s''', ...
                figure_prefix, reference_file_name ), ...
       'Interpreter', 'None' )

% Zoom into a region of differences according to the user's request and do
% some contortions to get our tick labels set.  We have to tighten the
% X-axis to match our data and then identify which tolerances have a tick
% associated with them.  Thanks to MATLAB's lovely plotting mechanics we
% also have to deal with ticks that exist outside of the X-axis.
ylim( difference_range );
xlim( [tolerance_range(1), tolerance_range(end)] );
ticks  = get( gca, 'xtick' );
[~, tolerance_indices, tick_indices] = intersect( tolerance_range, ticks );
labels = cellfun( @(x) sprintf( '10^%d%%', x ), num2cell( tolerance_range(tolerance_indices) + 2 ), ...
                  'UniformOutput', false );
set( gca, ...
     'xtick',      ticks(tick_indices), ...
     'xticklabel', labels );

xlabel( 'Tolerance' )
ylabel( 'Difference (%)' )

% NOTE: We modify the legend's properties directly as Octave 4.x does not
%       handle specifying them during creation.
h_legend = legend( compared_file_names{:} );
set( h_legend, ...
     'Location',    'NorthEast', ...
     'Interpreter', 'None' )

% Suppress output when outputs were not explicitly requested.
if nargout < 1
    clear differences fig_h;
end

end
