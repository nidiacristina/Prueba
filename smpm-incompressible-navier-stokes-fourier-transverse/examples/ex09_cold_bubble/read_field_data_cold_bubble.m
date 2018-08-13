clear
clc
close all

addpath('../../smpm_matlab_utilities/')

%Read the files with results inside the directory

fileList = dir('*out*.h5');
ent=zeros(length(fileList),1);
ens=zeros(length(fileList),1);
pal=zeros(length(fileList),1);
ttime=zeros(length(fileList),1);

for i=1:length(fileList)
%for i=1
    data = smpm_read_fieldfile(fileList(i).name,'NEW');
    %contourf(data.grid.x,data.grid.z,data.field.uz)
   
    ttime(i)= data.field.time;
    ent(i)= compute_kinetic_energy(data);
    ens(i)= compute_enstrophy(data);
    pal(i) = compute_palinstrophy(data);
    
    clf
%    plot_streamline(x,y,u,v,t)

    plot_vorticity(data)  
    %pause(0.001)
    %pause
    saveas(gcf,['NR_5000_Voticity_',num2str(data.field.time*100,2),'_Mx_',num2str(data.grid.mx,2),'_Mz_',num2str(data.grid.mz,2),'_i_',num2str(i,2),'.png'])

end

%{
saveas(gcf,['NR_2500_Voticity_',num2str(data.field.time,2),'_Mx_',num2str(data.grid.mx,2),'_Mz_',num2str(data.grid.mz,2),'.png'])

figure
plot(ttime,ent,'b')
xlabel('t')
ylabel('E (t)')
saveas(gcf,['NR_2500_KineticEnergy_N_',num2str(data.grid.n,2),'_Mx_',num2str(data.grid.mx,2),'_Mz_',num2str(data.grid.mz,2),'.png'])

figure
plot(ttime,ens,'g')
xlabel('t')
ylabel('\Omega (t)')
saveas(gcf,['NR_2500_Enstrophy_N_',num2str(data.grid.n,2),'_Mx_',num2str(data.grid.mx,2),'_Mz_',num2str(data.grid.mz,2),'.png'])

%}
