% Example 06: Builds the inputs for the linear advection-diffusion problem, 
% using a Gaussian bump that propagates in the periodic irection  All files are
% written to the current directory. 

% Add the SMPM Matlab directory to your path (just in case you haven't already
% done this).
addpath( genpath( '../../' ) );

[input_file_name, inputs] = create_linear_advection_diffusion_inputs();
init                      = smpm_read_initfile( inputs.fname_init );

