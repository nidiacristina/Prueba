clc
clear
close all

%addpath ~/Prueba/smpm-incompressible-navier-stokes-fourier-transverse/smpm_matlab_utilities/

%Read the files with results inside the directory

fileList = dir('*lid_driven_cavity_out_00002*.h5');

for i=1:length(fileList)
    delete fileList(i).name
end