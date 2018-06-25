function [wall_time, start_time] = smpm_get_execution_statistics( field_file_name, print_flag )
% [wall_time, start_time] = smpm_get_execution_statistics( field_file_name, print_flag )
%
% Extracts the execution statistics from the specified SMPM field file and
% optionally prints a summary to standard output.
%
% Takes 2 arguments:
%
%   field_file_name - Path to the SMPM field file.
%   print_flag      - Optional flag specifying that the execution statistics
%                     should be printed to standard output.  If omitted,
%                     defaults to false unless no outputs are requested in
%                     which cause it defaults to true.
%
% Returns 2 values:
%
%   wall_time  - Structure containing execution statistics.
%   start_time - Start time of the simulation, in seconds since the epoch.
%

if nargin < 2
   print_flag = [];
end

% disable printing unless either explicitly requested or no outputs were
% provided.
if isempty( print_flag )
   if nargout == 0
      print_flag = true;
   else
      print_flag = false;
   end
end

wall_time = struct( );

% read a subset of the field output file and create a structure containing the
% important parts of the /execution/ sub-group.
start_time                  = h5read( field_file_name, '/execution/wall_time/start' );
wall_time.elapsed           = h5read( field_file_name, '/execution/wall_time/total' );

wall_time.setup.null_basis  = h5read( field_file_name, '/execution/wall_time/null_basis' );
wall_time.setup.null_error  = h5read( field_file_name, '/execution/wall_time/null_error' );
wall_time.setup.total       = h5read( field_file_name, '/execution/wall_time/setup' );

wall_time.steps             = h5read( field_file_name, '/execution/wall_time/steps' );
wall_time.io.field          = h5read( field_file_name, '/execution/wall_time/field_io' );
wall_time.io.restart        = h5read( field_file_name, '/execution/wall_time/restart_io' );

wall_time.compute.diffusion = h5read( field_file_name, '/execution/wall_time/solve_diffusion' );
wall_time.compute.poisson   = h5read( field_file_name, '/execution/wall_time/solve_poisson' );
wall_time.compute.viscous   = h5read( field_file_name, '/execution/wall_time/solve_viscous' );

% print the vitals to standard output if requested.
if print_flag
    timesteps_time         = sum( wall_time.steps );
    timesteps_io_time      = wall_time.io.field + wall_time.io.restart;
    timesteps_compute_time = timesteps_time - timesteps_io_time;

    fprintf( 1, 'Execution started at %s and ran for %f seconds.\n', ...
             datestr( datenum( 1970, 1, 1 ) + start_time / 86400, 31 ), ...
             wall_time.elapsed );
    fprintf( 1, '\n' )

    # setup.
    fprintf( 1, 'Setup time:\n' )
    fprintf( 1, '\n' )
    fprintf( 1, '   Null basis:    %.3f seconds\n', wall_time.setup.null_basis );
    fprintf( 1, '   Null error:    %.3f seconds\n', wall_time.setup.null_error );
    fprintf( 1, '\n' )
    fprintf( 1, '   Total:         %.3f seconds\n', wall_time.setup.total );
    fprintf( 1, '\n' )

    # timesteps.
    fprintf( 1, 'Solver time:\n' )
    fprintf( 1, '\n' )
    fprintf( 1, '   I/O:           %.3f seconds (%5.2f%%)\n', ...
             timesteps_io_time, timesteps_io_time / timesteps_time * 100 )
    fprintf( 1, '      Field:          %.3f seconds (%5.2f%%)\n', ...
             wall_time.io.field, wall_time.io.field / timesteps_time * 100 )
    fprintf( 1, '      Restart:        %.3f seconds (%5.2f%%)\n', ...
             wall_time.io.restart, wall_time.io.restart / timesteps_time * 100 )
    fprintf( 1, '\n' )
    fprintf( 1, '   Compute:       %.3f seconds (%5.2f%%)\n', ...
             timesteps_compute_time, timesteps_compute_time / timesteps_time * 100 )
    fprintf( 1, '      Diffusion:      %.3f seconds (%5.2f%%)\n', ...
             wall_time.compute.diffusion, wall_time.compute.diffusion / timesteps_time * 100 )
    fprintf( 1, '      Poisson:        %.3f seconds (%5.2f%%)\n', ...
             wall_time.compute.poisson, wall_time.compute.poisson / timesteps_time * 100 )
    fprintf( 1, '      Viscous:        %.3f seconds (%5.2f%%)\n', ...
             wall_time.compute.viscous, wall_time.compute.viscous / timesteps_time * 100 )
    fprintf( 1, '\n' )
end

% do not return outputs if they weren't explicitly requested.
if nargout == 0
   clear wall_time start_time;
end

end
