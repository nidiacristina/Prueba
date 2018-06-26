clear all
clc
close all

addpath ~/Prueba/smpm-incompressible-navier-stokes-fourier-transverse/smpm_matlab_utilities/

data = smpm_read_fieldfile('lock_exchange_out_000001.h5','NEW');

figure(1)
contourf(data.grid.x,data.grid.z,data.field.ux)