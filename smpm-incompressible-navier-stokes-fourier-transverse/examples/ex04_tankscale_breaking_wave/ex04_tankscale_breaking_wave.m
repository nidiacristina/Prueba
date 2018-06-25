% Example 04: Builds the input files for the tank-scale internal wave
% propagation problem with wave velocities artificially increased to induce
% breaking.

% Add the SMPM Matlab directory to your path (just in case you haven't already
% done this).
addpath( genpath( '../../' ) );

[input_file_name, inputs] = create_tankscale_breaking_wave_inputs();
