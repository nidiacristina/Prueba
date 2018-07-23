% File that Creates the second part of the runs in excalibur
%
% Gustavo Rivera
% December2017

clear all
clc
close all
clear cmd

% Add path to matlab utlities
%addpath ~/RAID/Codes/smpm_navier_stokes/smpm_navier_stokes_3D/fourier-transverse-gfortran-mkl/smpm_matlab_utilities/
addpath ~/Prueba/smpm-incompressible-navier-stokes-fourier-transverse/smpm_matlab_utilities/

%Original folder
dir0='../examples/ex02_lock_exchange/';

% Specify path to restart file
%path_to_restart = (['../part_1/Shoaling_ISW_Dongsha_res_000020.h5']);
path_to_restart = (['../examples/ex02_lock_exchange/lock_exchange_restart_000047.h5']);

% Read the restart file
restart = smpm_read_restartfile(path_to_restart);

% Specify Path to input file
path_to_input = [dir0,'lock_exchange_in'];

% Read input file
inputs = smpm_read_inputfile(path_to_input);

% Create name for restart file
restart_name = 'restart_part_3.h5';

% Specify Setup File name
setup_file_name = 'lock_exchange_setup.h5';

% Give this run a name and specify its initial conditions file.
inputs.fname_setup                    = setup_file_name;
inputs.fname_restart                  = restart_name;

% Specify input file name
input_file_name                       = 'lock_exchange_in';

% Set the setupfile read flag
inputs.read_from_setupfile            = logical( 0 );
inputs.write_to_setupfile             = logical( 1);
inputs.apply_restart                  = logical( 1);


% Update simulation time (in seconds)
inputs.tend                           = 20; 

inputs.timesteps_between_writes       = round((inputs.tend/inputs.dt)/20,0);
inputs.timesteps_between_logs         = round((inputs.tend/inputs.dt)/20,0);
inputs.timesteps_between_restarts     = round((inputs.tend/inputs.dt)/20,0)*2;

% Create input file
smpm_write_inputfile( input_file_name, inputs );

% Create init file from restart field
smpm_create_file_from_checkpoint(restart, restart_name);

%Copy files created in directory 
copyfile input_file_name dir0
