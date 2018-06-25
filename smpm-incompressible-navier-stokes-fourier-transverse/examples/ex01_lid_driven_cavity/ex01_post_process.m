%
% Example 01: post processing example.  This script assumes you've run the example
% shown in ex01_lid_driven_cavity.m.

% Add the SMPM Matlab directory to your path (just in case you haven't already done this).
addpath( genpath( '../../' ) );

% Set the time-index in the field file you want to post-process.
time_ndx = -1; % -1 means last index, -2 means second from last, etc.

% Set some file names of the run that was just completed.
field_file_name     = 'lid_driven_cavity_out.h5';
init_file_name      = 'lid_driven_cavity_init.h5';
input_file_name     = 'lid_driven_cavity_in';

% Set some output file names for the post-processing code.
post_init_file_name  = 'lid_driven_cavity_post_init.h5';
post_out_name        = 'lid_driven_cavity_post';
post_input_file_name = 'lid_driven_cavity_post_in';

%% Read in a bunch of data from the SMPM simulation.

   % Read in the last time-step of the SMPM simulation.
   data = smpm_read_fieldfile( field_file_name, time_ndx );

   % Read the initial conditions file for the SMPM simulation.
   init = smpm_read_initfile( init_file_name );

   % Read the input file for the SMPM simulation.
   input = smpm_read_inputfile( input_file_name );

%% Write input files for the SMPM post-processor.

   % Write an initial conditions file for the post-processor that has the last time-step of data.
   smpm_write_initfile( init.grid.n, init.grid.mx, init.grid.my, init.grid.mz, ...
       init.grid.x, init.grid.y, init.grid.z, ...
       data.field.rho, data.field.ux, data.field.uy, data.field.uz, ...
       init.ic.ubc, init.ic.dubcdz, init.ic.rho_bar, init.ic.rho_bar_z, post_init_file_name );

   % Write an input file for the post-processor with the correct file names.
   input.fname_init    = post_init_file_name;
   input.fname_runname = post_out_name;
   smpm_write_inputfile( post_input_file_name, input );

