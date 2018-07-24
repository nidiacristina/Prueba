% Example 01: Builds the inputs for the lid-driven cavity test case.  All
% files are written to the current directory.

% Add the SMPM Matlab directory to your path (just in case you haven't already
% done this).
addpath( genpath( '../../' ) );

[input_file_name, inputs] = create_vortex_dippole_collision_inputs();
