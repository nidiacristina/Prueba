function [] = plot_vorticity(data)
% plot_vorticity(x,y,u,v,t)
%
% Plots the vorticity of a 2D velocity field using MATLAB functions
%
% Takes 1 input argument:
%   
%   data - Field structure, as read by smpm_read_fieldfile, containing
%           the field data (i.e. grid variables, field variables, time, 
%           etc...) for the run of interest.
%
% Returns nothing.
%
%
% Jul 2018
% NCRG

%1) Prepare de data to be plotted
    %Remove de columns with the boundary nodes (...to improve this)
    %{
    
    x2=x;
    z2=z;
    ux2=ux;
    uz2=uz;
    dif=x(1,1:(end-1))-x(1,2:(end));
    [r,cols]=find(dif==0);
    dif2=z(1:(end-1),1)-y(2:(end),1);
    [rows,c]=find(dif2==0);
  
    x2(rows,:)=[];
    z2(rows,:)=[];
    ux2(rows,:)=[];
    uz2(rows,:)=[];    
    
    x2(:,cols)=[];
    z2(:,cols)=[];
    ux2(:,cols)=[];
    uz2(:,cols)=[]; 
    
    cav = curl(x2,z2,ux2,uz2);
    
    %}
    
%2) Create variables to be used by the MATLAB function 
     
    %fig=figure
    %clf(fig)  
    
    x   = data.grid.x;
    y   = data.grid.y;
    z   = data.grid.z;
    n   = data.grid.n;
    mx  = data.grid.mx;
    my  = data.grid.my;
    mz  = data.grid.mz;

    % Grab Field Data
    ux  = data.field.ux + data.ic.ubc;
    uz  = data.field.uz;
    
    % Differentiate the density field
    [~, ~, dux_z ] = smpm_compute_gradient( ux, n, mx, my, mz, x, y, z);
    [duz_x, ~, ~ ] = smpm_compute_gradient( uz, n, mx, my, mz, x, y, z);
    
    vorty=duz_x-dux_z;
    %[ phi_x, phi_y, phi_z ] = smpm_compute_gradient( phi, n, 1, 1, mz, x, zeros(size(x)), zeros(size(x)) )
    
    pcolor(x,z,vorty); 
    shading flat
    %hold on
    %contour(x,y,cav)
    %quiver(x,y,u,v,'k');
    %hold off
    colormap('copper');
    %colormap(hot(8));
    %hold on
    
    %t=data.field.time;
    title(['Time  ',num2str(data.field.time)])
    %pause(0.1)
    %hold on
end

