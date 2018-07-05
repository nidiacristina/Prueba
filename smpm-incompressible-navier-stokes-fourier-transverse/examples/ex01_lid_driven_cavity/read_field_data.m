clear
clc
close all

addpath ~/Prueba/smpm-incompressible-navier-stokes-fourier-transverse/smpm_matlab_utilities/

%Read the files with results inside the directory

fileList = dir('*out*.h5');

for i=1:length(fileList)
%for i=1
    data = smpm_read_fieldfile(fileList(i).name,'NEW');
    %contourf(data.grid.x,data.grid.z,data.field.uz)
    
    %%{
    x=data.grid.x;
    y=data.grid.z;
    u=data.field.ux;
    v=data.field.uz;
    
    dif=x(1,1:(end-1))-x(1,2:(end));
    [r,cols]=find(dif==0);
    
    dif2=y(1:(end-1),1)-y(2:(end),1);
    [rows,c]=find(dif2==0);
    
    x(rows,:)=[];
    y(rows,:)=[];
    u(rows,:)=[];
    v(rows,:)=[];    
    
    x(:,cols)=[];
    y(:,cols)=[];
    u(:,cols)=[];
    v(:,cols)=[]; 
    
    %%{
    %figure
    %quiver(x,y,u,v)

    %sy = 0:0.1:1;
    %sx = sy;
    sx = diag(x);
    sy = diag(y);   
    
    %[sx,sy] = meshgrid(0:1, 0:1);
    clf
    streamline(x,y,u,v,sx,sy)
    %streamline(x,y,u,v,0,1)
    
    pause(0.1)
    %pause
    %}
    %}
    %{
    X=x;
    Y=y;
    vx=u;
    vy=v;
    
    %figure
    pcolor(X,Y,hypot(vx,vy))
    shading interp

    sc = 1/10; 
    hold on
    quiver(imresize(X,sc),imresize(Y,sc),imresize(vx,sc),imresize(vy,sc),'r') 
    pause(0.2)
    %}

    
end

%figure(1)
