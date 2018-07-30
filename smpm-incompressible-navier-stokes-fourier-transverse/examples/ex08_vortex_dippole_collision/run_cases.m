
clear
clc
close all

%addpath ~/Prueba/smpm-incompressible-navier-stokes-fourier-transverse/smpm_matlab_utilities/
%addpath ../smpm_matlab_utilities/


n     = 10;     % The number of GLL points per direction per element.
nsubx = 24;    % The number of x elements.
nsuby = 1;     % The number of y elements.
nsubz = 24;    % The number of z elements.

p_n=[8,10];
p_nsubx=[6,8]*8;

for k=1:length(p_n)  
    n=p_n(k);       
    for j=1:length(p_nsubx)

        nsubx=p_nsubx(j);
        nsubz=nsubx;
        %Remove files resulting from the previous simulations
        delete *vortex_dippole_collision_out_0000
        delete *vortex_dippole_collision_restart_0000
        %Builds the inputs for the case

        %%%Crear input
        [input_file_name, inputs] = create_vortex_dippole_collision_inputs(n, nsubx,nsuby,nsubz);

        %%%%Correr codigo
        command = './run_ex08_vortex_dippole_collision.sh';
        [status,cmdout] = system(command,'-echo')
        status

        try
            %%%%Leer resultados

            fileList = dir('*out*.h5');
            ent=zeros(length(fileList),1);
            ens=zeros(length(fileList),1);
            pal=zeros(length(fileList),1);
            ttime=zeros(length(fileList),1);

            for i=1:length(fileList)
                data = smpm_read_fieldfile(fileList(i).name,'NEW');
                ttime(i)= data.field.time;
                ent(i)= compute_kinetic_energy(data);
                ens(i)= compute_enstrophy(data);
                pal(i)= compute_palinstrophy(data);
            end

            save(['NR_2500_Variables_N_',num2str(n,2),'_Mx_',num2str(nsubx,2),'_Mz_',num2str(nsubz,2),'.mat'],'ent','ens','ttime')

            figure
            plot(ttime,ent,'b')
            xlabel('t')
            ylabel('E (t)')
            saveas(gcf,['NR_2500_KineticEnergy_N_',num2str(n,2),'_Mx_',num2str(nsubx,2),'_Mz_',num2str(nsubz,2),'.png'])

            figure
            plot(ttime,ens,'g')
            xlabel('t')
            ylabel('\Omega (t)')
            saveas(gcf,['NR_2500_Enstrophy_N_',num2str(n,2),'_Mx_',num2str(nsubx,2),'_Mz_',num2str(nsubz,2),'.png'])

            figure
            plot(ttime,pal,'r')
            xlabel('t')
            ylabel('\P (t)')
            saveas(gcf,['NR_2500_Palinstrophy_N_',num2str(n,2),'_Mx_',num2str(nsubx,2),'_Mz_',num2str(nsubz,2),'.png'])
            close all
        
        catch
            
        end
    end
 end