function [] = energy_calculations(x,y,u,v,t)
% plot_vorticity(x,y,u,v,t)
%
% Calculates kinetic energy
%
% Takes 5 arguments:
%
%   x              - Vector, of length n^2 * mx * mz * my, containing the
%                    x-coordinates of the mesh associated with the field
%                    variables.
%   y              - Vector, of length n^2 * mx * mz * my, containing the
%                    z-coordinates of the mesh associated with the field
%                    variables.
%   u            - Scalar, or vector, specifying the velocity in
%                    the x-direction of the field. Dimension n^2 * mx * mz * my.
%   v            - Scalar, or vector, specifying the velocity in
%                    the z-direction of the field. Dimension n^2 * mx * mz * my.
%   t            - Scalar specifying time.
%
% Returns nothing.
%
%
% Jul 2018
% NCRG

%1) Prepare de data to be plotted
    %Remove de columns with the boundary nodes (...to improve this)
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

%2) Calculate kinetic energy 
     
    %fig=figure
    %clf(fig)

    en=0.5*(u+v).^0.5;
    cav = curl(x,y,u,v);
    pcolor(x,y,cav); 
    shading interp
    %hold on
    %contour(x,y,cav)
    %quiver(x,y,u,v,'k');
    %hold off
    colormap('copper');
    %colormap(hot(8));
    %hold on
    
    %t=data.field.time;
    title(['Time  ',num2str(t)])
    %pause(0.1)
    %hold on
end

