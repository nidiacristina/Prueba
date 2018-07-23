clear
clc
close all

addpath ~/Prueba/smpm-incompressible-navier-stokes-fourier-transverse/smpm_matlab_utilities/

%Read the files with results inside the directory

fileList = dir('*out*.h5');
ent=zeros(length(fileList),1);

for i=1:length(fileList)
%for i=1
    data = smpm_read_fieldfile(fileList(i).name,'NEW');
    
    %contourf(data.grid.x,data.grid.z,data.field.uz)
    
    x=data.grid.x;
    y=data.grid.z;
    u=data.field.ux;
    v=data.field.uz;
    t=data.field.time;
    
    ent(i)=sum(sum(0.5*(u.^2+v.^2).^0.5));
    
    %clf('reset')
    %h=subplot(2,2,3);
    %reset(h)
    clf
    plot_streamline(x,y,u,v,t)
%     j=subplot(2,2,4);
%     reset(j)
%     plot_vorticity(x,y,u,v,t)  
    %pause(0.001)
       
end


saveas(gcf,['NR_250_Stream_lines_',num2str(data.field.time,2),'.png'])

%figure(1)
