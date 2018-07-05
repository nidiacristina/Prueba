clear
clc
close all

addpath ~/Prueba/smpm-incompressible-navier-stokes-fourier-transverse/smpm_matlab_utilities/

%Read the files with results inside the directory

fileList = dir('*out*.h5');
%fileList = dir('*restart*.h5');

for i=1:length(fileList)
    data = smpm_read_fieldfile(fileList(i).name,'NEW');
    contourf(data.grid.x,data.grid.z,data.field.ux)
    pause
    %pause(0.2)
end

%figure(1)
