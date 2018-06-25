function [execute_flag, input_file_name] = smpm_prepare_test_case( test_path, np, force_flag, quiet_flag )
% [execute_flag, input_file_name] = smpm_prepare_test_case( test_path, np, force_flag, quiet_flag )
%
% Prepares a test case for execution by instantiating its inputs and
% determining whether it could be invoked on the current system with a
% specific number of processors.  Several indicators are provided allowing
% this function to be called programmatically either from within MATLAB or
% from another language that parses this function's output.
%
% Takes 4 arguments:
%
%   test_path  - Path to the SMPM input file generator associated with the test
%                case.
%   np         - Number of processors used during execution of the test case.
%   force_flag - Optional flag indicating that the specified test case should
%                be run regardless of memory availability.  If omitted, defaults
%                to false and does not force execution.
%   quiet_flag - Optional flag indicating that not output should be generated.
%                For use when invoking the function non-interactively from
%                within MATLAB.  If omitted, defaults to false and generates
%                output.
%
% Returns 2 values:
%
%   execute_flag    - Flag indicating whether the test case can be executed
%                     without exhausting available memory.
%   input_file_name - Path to the instantiated test cases' input file.
%

% handle our optional flags.
if nargin < 4
    quiet_flag = [];
end
if nargin < 3
    force_flag = [];
end

if isempty( force_flag )
    force_flag = false;
end
if isempty( quiet_flag )
    quiet_flag = false;
end

% figure out where the test is and move into that directory, so it is in our
% path.
[directory, file_name, file_extension] = fileparts( test_path );

% run the test to build the inputs.
%
% NOTE: we disable warnings to avoid MATLAB from helpfully telling us that the
%       function we're executing does not match the file name it resides in.
%       this can happen in some test cases that are symbolic links to examples.
%
original_pwd = cd( directory );
warning( 'off' )
[input_file_name, inputs] = feval( file_name );
warning( 'on' )
cd( original_pwd );

% run the memory checker to see if we fit.
[needed_memory, needed_per_node] = smpm_check_mem_requirements( inputs.n, ...
                                                                inputs.nsubx, ...
                                                                inputs.nsuby, ...
                                                                inputs.nsubz, ...
                                                                np, ...
                                                                inputs.gmres_maxit_viscous, ...
                                                                inputs.gmres_maxit_poisson, ...
                                                                true );
available_memory = get_available_memory();

% bail if we don't have enough memory and aren't being forced to execute.
%
% NOTE: we assume all of our ranks are local.
execute_flag = (needed_memory < available_memory) | force_flag;
if ~execute_flag
    return
end

% give the go/nogo by printing the input file name if we're not being quiet.
if ~quiet_flag
    fprintf( '%s\n', input_file_name );
end

end

function available = get_available_memory()
% available = get_available_memory()
%
% Determines the total amount of available memory on the system.  The value
% returned is total RAM not including swap and does not attempt to compute
% the amount of free memory.
%
% Takes no arguments.
%
% Returns 1 value:
%
%   available - The size of available memory on the current computer, in
%               bytes.  Is zero if the amount of memory available cannot be
%               acquired.
%

% MATLAB's memory command is basically worthless as it doesn't tell us useful
% information *AND* isn't available everywhere.  so we detect what type of
% system we're on and figure things out ourselves.
if isunix && ~ismac
    % Linux: free's memory capacity is listed in MiB.
    [status, output] = system( 'free | grep ^Mem | awk ''{ print $2 }'' ' );
    available = str2num( output ) * 1024;
else
    % if we don't know how to query the system's memory capacity, assume we
    % can't run anything.
    available = 0;
end

end
