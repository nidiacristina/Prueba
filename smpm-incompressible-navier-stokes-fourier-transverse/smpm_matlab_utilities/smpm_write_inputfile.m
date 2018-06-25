function smpm_write_inputfile( input_file_name, data )
% smpm_write_inputfile( input_file_name, data )
%
% Writes the contents of data to an input file for the SMPM code.
%
% Takes 2 arguments:
%
%   input_file_name - String indicating the field file to read from disk.
%   data            - Struct of variable name/values to write (c.f. as
%                     returned by smpm_read_inputfile()).
%
% Returns nothing.
%

% 11 Nov 2014
% Sumedh Joshi

% The SMPM solver's input file parsing and validation is not terribly
% sophisticated and silently truncates long lines.  Warn the user if
% we create an input file that has long lines so they can preemptively
% deal with the situation.
max_line_length = 132;

% Specify the order of fields to write, and whether they're mandatory for
% producing a minimal input file.  Any fields specified in data will
% be written with these first, and in the order below, with the rest in a
% random order.  All of the fields that are mandatory must be present in
% data.
fields_config = { ... % key file names in the simulation.
                  'fname_runname',                   true;  ...
                  'fname_restart',                   false; ...
                  'fname_init',                      true;  ...
                  'fname_setup',                     false; ...
                  ... % parameters governing bootstrapping a simulation.
                  'read_from_setupfile',             false;  ...
                  'write_to_setupfile',              false; ...
                  'apply_restart',                   false; ...
                  'setup_and_stop',                  false; ...
                  ... % grid parameters.
                   'n',                              true;  ...
                   'nsubx',                          true;  ...
                   'nsuby',                          true;  ...
                   'nsubz',                          true;  ...
                   ... % simulation time parameters.
                   'dt',                             true;  ...
                   'tend',                           true;  ...
                   'adaptive_timestep',              false; ...
                   ... % frequency of output, in various forms.
                  'timesteps_between_writes',        false; ...
                  'timesteps_between_restarts',      false; ...
                  'timesteps_between_logs',          false; ...
                   ... % physical (viscosity, density, etc) and simulation parameters.
                   'nu',                             true;  ...
                   'nu_d',                           true;  ...
                   'rho_0',                          true;  ...
                   'facrobin',                       true;  ...
                   'facrobin_ppe',                   true;  ...
                   'filter_order_xz',                true;  ...
                   'filter_order_y',                 true;  ...
                   ... % sub-problem iterative solution tolerances and thresholds.
                   'gmres_maxit_poisson',            false; ...
                   'gmres_maxit_viscous',            false; ...
                   'gmres_tol_poisson',              false; ...
                   'gmres_tol_viscous',              false; ...
                   'gmres_restart_viscous',          false; ...
                   'gmres_restart_poisson',          false; ...
                   ... % boundaries.
                   'sponge_layer_location',          false; ...
                   'bc_diffusion',                   true;  ...
                   'bc_viscous',                     true;  ...
                   'read_bcs_from_initfile',         false; ...
                   'apply_sponge_layer',             false; ...
                   'left_fraction',                  false; ...
                   'right_fraction',                 false; ...
                   'enforce_strong_velocity_bc',     false; ...
                   ... % knobs for developing the solver.
                   'check_null_error',               false; ...
                   'check_numerical_error',          false; ...
                   'do_interfacial_averaging',       false; ...
                   'use_capacitance_preconditioner', false; ...
                   'exact_nullspace_projection',     false; ...
                   'use_deflation',                  false; ...
                 };

% Get a list of fields of the data struct and determine the longest so we
% can left justify the output correctly.
F          = fieldnames( data );
pad_length = max( cellfun( @length, F ) );

% Compute the order in which we walk through our input structure's fields
% so that the output follows a logical ordering.  Any additional fields
% that aren't in fields_config will be handled separately afterwards.
ordered_indices = cellfun( @(x) find( strcmp( x, F ) ), { fields_config{:, 1} }, ...
                           'UniformOutput', false );
empty_indices   = cellfun( @isempty, ordered_indices );
ordered_indices = [ordered_indices{~empty_indices}];

% Ensure that the configuration we're writing contains at least all of the
% mandatory fields.
mandatory_indices = find( [fields_config{:, 2}] );
available_indices = find( ~empty_indices );

missing_indices = find( ~ismember( mandatory_indices, available_indices ) );
if length( missing_indices ) > 0
    missing_fields = '';

    % Build a list of fields that are missing to help the user out.
    for missing_index = missing_indices
        missing_fields = sprintf( '%s\n    %s', ...
                                  missing_fields, ...
                                  fields_config{mandatory_indices(missing_index), 1} );
    end
    error( 'Supplied input structure is missing the following mandatory fields:%s\n', ...
           missing_fields )
end

% Open a file stream.
fid = fopen( input_file_name, 'w+' );

% Write a time-stamp at the top of the file.
fprintf( fid, '#####################\n' );
fprintf( fid, '# Input file for SMPM code.\n# Generated: %s\n', datestr( now ) );
fprintf( fid, '#####################\n' );

% Get a cell array of the values of the data struct.
V = struct2cell( data );

% Loop over the data/value pairs, writing to the input file.
for ii = ordered_indices
    write_configuration( fid, F{ii}, V{ii}, pad_length, max_line_length );
end

% Remove all of the data/value pairs we just wrote so we can handle
% any that remain.
F(ordered_indices) = [];
V(ordered_indices) = [];

% Handle any configuration parameters that aren't in the ordered_fields
% list above.
for ii = 1:numel( F )
    write_configuration( fid, F{ii}, V{ii}, pad_length, max_line_length );
end

% Close the file.
fclose( fid );

return

function write_configuration( fid, field, value, pad_length, length_threshold )
% write_configuration( fid, field, value, pad_length, length_threshold )
%
% Writes a single configuration field to the file descriptor specified.  The
% field name will be padded with spaces (that is, left justified) before the
% equals sign and value are written.  Care is taken to identify the value's
% type and write out a value appropriate for reading within the solver.
%
% Takes 5 values:
%
%   fid              - File descriptor to write the configuration to.
%   field            - String specifying the name of the configuration.
%   value            - Value to be written.
%   pad_length       - Width of the field to write.  field will be left
%                      justified, padded with spaces, with pad_length.
%   length_threshold - Optional line length used to warn the caller of long
%                      lines when exceeded.  No length check is performed if
%                      specified as 0.  If omitted, defaults to 0.
%
% Returns nothing.

% Number of bytes this (field, value) pair took when written.
line_length = 0;

% NOTE: We check for logical values first since Octave treats logicals
%       as numeric 1's and 0's.  If this case wasn't considered before
%       a pure numeric, logicals would cause the SMPM solver to bark
%       when it parses the input file.
if islogical( value )
    % This is a logical - print true or false.
    if value
        line_length = fprintf( fid, '%-*s = .true.\n', pad_length, field );
    else
        line_length = fprintf( fid, '%-*s = .false.\n', pad_length, field );
    end
elseif isnumeric( value )
    % If this is just one number, print and move on.
    if length( value ) == 1
        line_length = fprintf( fid, '%-*s = %g\n', pad_length, field, value );
    elseif length( value ) == 2
        % This is an array of length two, which is probably
        % domain boundaries.
        line_length = fprintf( fid, '%-*s = %g,%g\n', ...
                               pad_length, field, value );
    else
        % This is an array, print it in Fortran format.
        %
        % NOTE: Currently this assumes that we're printing out the 4
        %       external boundaries.
        line_length = fprintf( fid, '%-*s = %g,%g,%g,%g\n', ...
                               pad_length, field, value );
    end
elseif ischar( value )
    % If it's a string, just print it.
    line_length = fprintf( fid, '%-*s = %s\n', pad_length, field, value );
end

% NOTE: We compare against one more than requested to account for the
%       newline.  These do not count against line length in the SMPM
%       input file parser.
if (length_threshold > 0) && (line_length > (length_threshold + 1))
    warning( ['Field ''%s'' requires a line of length %d which exceeds %d characters!\n' ...
              'This will likely be truncated and misparsed.  Please review.\n'], ...
             field, line_length, length_threshold )
end

return

